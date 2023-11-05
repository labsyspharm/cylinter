import os
import sys
import re
import glob
import pickle
import logging
from dataclasses import dataclass
from typing import Dict

import numpy as np
import pandas as pd

import math

import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib import colors

from sklearn.preprocessing import MinMaxScaler

import zarr
import dask.array as da
import tifffile
import skimage

## these imports are for classical artifact detection
import napari
from skimage.filters.rank import minimum as rank_min
from skimage.filters.rank import maximum as rank_max
from skimage.filters.rank import gradient, mean
from skimage.measure import label 

from skimage.morphology import disk, erosion, dilation, h_maxima, flood, local_maxima
from skimage.exposure import rescale_intensity
############

from napari.utils.notifications import notification_manager, Notification, NotificationSeverity

logger = logging.getLogger(__name__)

SUPPORTED_EXTENSIONS = ['.csv']


def log_banner(log_function, msg):
    """Call log_function with a blank line, msg, and an underline."""
    log_function("")
    log_function(msg)
    log_function("-" * len(msg))


def log_multiline(log_function, msg):
    """Call log_function once for each line of msg."""
    for line in msg.split("\n"):
        log_function(line)


def input_check(self):

    # check for redundant sampleMetadata keys
    if len(set(self.sampleNames.keys())) != len(self.sampleNames.keys()):
        logger.info('Aborting; sampleMetadata contains redundant keys.')
        sys.exit()
    
    contents = os.listdir(self.inDir)

    # check whether input directory contains expected files and folders:
    if any(element not in contents for element in
           ['config.yml', 'markers.csv', 'csv', 'tif', 'seg', 'mask']):

        ##########################################################################################
        # if not, check for mcmicro input directory
        
        # test for TMA data
        if 'dearray' in contents:
            try:
                markers = pd.read_csv(os.path.join(self.inDir, 'markers.csv'))
            except FileNotFoundError:
                logger.info('Aborting; markers.csv file not found.')
                sys.exit()
            
            for key in self.sampleNames.keys():

                sample_name = key.split('--')[0]
                segmentation_method = key.split('--')[1].split('_')[0]
                segmentation_object = key.split('--')[1].split('_')[1]

                try:
                    glob.glob(
                        os.path.join(self.inDir, 'quantification', f'{key}*.csv'))[0]
                except IndexError:
                    logger.info(
                        f'sampleMetadata key {sample_name} in config.yml does '
                        'not match a CSV filename.'
                    )
                    sys.exit()

                try:
                    glob.glob(
                        os.path.join(
                            self.inDir, 'dearray', f'{sample_name}*.tif'))[0] 
                except IndexError:
                    logger.info(f'Aborting; OME-TIFF for sample {sample_name} not found.')
                    sys.exit()

                try:
                    glob.glob(
                        os.path.join(self.inDir, 'qc/s3seg',
                                     f"{segmentation_method}-{sample_name}", 
                                     f"{segmentation_object}*.tif"))[0]
                except IndexError:
                    logger.info(f'Aborting; segmentation outlines for sample {sample_name} '
                                'not found.'
                                )
                    sys.exit()
                
                try:
                    glob.glob(
                        os.path.join(self.inDir, 'segmentation',
                                     f"{segmentation_method}-{sample_name}",
                                     f"{segmentation_object}*.tif"))[0]
                except IndexError:
                    logger.info(f'Aborting; segmentation mask for sample {sample_name} '
                                'not found.'
                                )
                    sys.exit()

            markers_filepath = os.path.join(self.inDir, 'markers.csv')
            return 'mcmicro_TMA', markers_filepath

        else:
            # test for WSI data
            # check that samples specified in config.yml each have a csv, tif, seg, and mask file
            markers_list = []
            for key in self.sampleNames.keys():
                
                sample_name = key.split('--')[0]
                segmentation_method = key.split('--')[1].split('_')[0]
                segmentation_object = key.split('--')[1].split('_')[1]

                try:
                    markers = pd.read_csv(os.path.join(self.inDir, sample_name, 'markers.csv'))
                    markers_list.append(markers)
                except FileNotFoundError:
                    logger.info(f'Aborting; markers.csv file for sample {sample_name} not found.')
                    sys.exit()
                
                try:
                    glob.glob(
                        os.path.join(self.inDir, sample_name, 'quantification', f'{key}*.csv'))[0]
                except IndexError:
                    logger.info(
                        f'sampleMetadata key {sample_name} in config.yml does '
                        'not match a CSV filename.'
                    )
                    sys.exit()
                
                try:
                    glob.glob(
                        os.path.join(
                            self.inDir, sample_name, 'registration', f'{sample_name}*.tif'))[0] 
                except IndexError:
                    logger.info(f'Aborting; OME-TIFF for sample {sample_name} not found.')
                    sys.exit()

                try:
                    glob.glob(
                        os.path.join(self.inDir, sample_name, 'qc/s3seg',
                                     f"{segmentation_method}-{sample_name}", 
                                     f"{segmentation_object}*.tif"))[0]
                except IndexError:
                    logger.info(f'Aborting; segmentation outlines for sample {sample_name} '
                                'not found.'
                                )
                    sys.exit()
                
                try:
                    glob.glob(
                        os.path.join(self.inDir, sample_name, 'segmentation',
                                     f"{segmentation_method}-{sample_name}",
                                     f"{segmentation_object}*.tif"))[0]
                except IndexError:
                    logger.info(f'Aborting; segmentation mask for sample {sample_name} '
                                'not found.')
                    sys.exit()

            # check that all markers.csv files are identical (if not, which is one is correct?)
            if not all(markers.equals(markers_list[0]) for markers in markers_list):
                logger.info('Aborting; markers.csv files differ between samples.')
                sys.exit()
            
            markers_filepath = os.path.join(
                self.inDir, list(self.sampleNames.keys())[0].split('--')[0], 'markers.csv'
            )
            return 'mcmicro_WSI', markers_filepath

        ##########################################################################################

    # next, check that csv, tif, seg, and mask subdirectories each contain files for all samples
    csv_names = set(
        [os.path.basename(path).split('.')[0] for path
         in glob.glob(os.path.join(self.inDir, 'csv', '*.csv'))]
    )

    tif_names = set(
        [os.path.basename(path).split('.')[0] for path
         in glob.glob(os.path.join(self.inDir, 'tif', '*.tif'))]
    )
    
    seg_names = set(
        [os.path.basename(path).split('.')[0] for path
         in glob.glob(os.path.join(self.inDir, 'seg', '*.tif'))]
    )
    
    mask_names = set(
        [os.path.basename(path).split('.')[0] for path
         in glob.glob(os.path.join(self.inDir, 'mask', '*.tif'))]
    )
    if not all(s == csv_names for s in [csv_names, tif_names, seg_names, mask_names]):

        return False
    
    # check that file names specified in config.yml are contained in input directory
    if not set(self.sampleNames.keys()).issubset(csv_names):
        logger.info(
            'Aborting; at least 1 sampleMetadata key in config.yml does not match a CSV filename.'
        )
        sys.exit()

    markers_filepath = os.path.join(self.inDir, 'markers.csv')
    return 'standard', markers_filepath


def get_filepath(self, check, sample, file_type):

    sampleMetadata_key = next(
        (key for key, val in self.sampleNames.items() if val == sample), None
    )

    if check == 'standard':
        if file_type == 'CSV':
            file_path = os.path.join(self.inDir, 'csv', f"{sampleMetadata_key}.csv")
        if file_type == 'TIF':
            file_path = glob.glob(
                os.path.join(self.inDir, 'tif', f"{sampleMetadata_key}.*tif"))[0]
        if file_type == 'SEG':
            file_path = glob.glob(
                os.path.join(self.inDir, 'seg', f"{sampleMetadata_key}.*tif"))[0]
        if file_type == 'MASK':
            file_path = glob.glob(
                os.path.join(self.inDir, 'mask', f"{sampleMetadata_key}.*tif"))[0]

    elif check == 'mcmicro_TMA':
        sample_name = sampleMetadata_key.split('--')[0]
        segmentation_method = sampleMetadata_key.split('--')[1].split('_')[0]
        segmentation_object = sampleMetadata_key.split('--')[1].split('_')[1]

        if file_type == 'CSV':
            file_path = os.path.join(
                self.inDir, 'quantification', f'{sampleMetadata_key}.csv'
            )
        if file_type == 'TIF':
            file_path = glob.glob(
                os.path.join(self.inDir, 'dearray', f"{sample_name}.*tif"))[0]
        if file_type == 'SEG':
            file_path = glob.glob(
                os.path.join(self.inDir, 'qc/s3seg',
                             f"{segmentation_method}-{sample_name}", 
                             f"{segmentation_object}*.tif"))[0]
        if file_type == 'MASK':
            file_path = glob.glob(
                os.path.join(self.inDir, 'segmentation',
                             f"{segmentation_method}-{sample_name}", 
                             f"{segmentation_object}*.tif"))[0]

    elif check == 'mcmicro_WSI':
        sample_name = sampleMetadata_key.split('--')[0]
        segmentation_method = sampleMetadata_key.split('--')[1].split('_')[0]
        segmentation_object = sampleMetadata_key.split('--')[1].split('_')[1]
        
        if file_type == 'CSV':
            file_path = os.path.join(
                self.inDir, sample_name, 'quantification', f'{sampleMetadata_key}.csv'
            )
        if file_type == 'TIF':
            file_path = glob.glob(
                os.path.join(self.inDir, sample_name, 'registration', f"{sample_name}.*tif"))[0]
        if file_type == 'SEG':
            file_path = glob.glob(
                os.path.join(self.inDir, sample_name, 'qc/s3seg',
                             f"{segmentation_method}-{sample_name}", 
                             f"{segmentation_object}*.tif"))[0]
        if file_type == 'MASK':
            file_path = glob.glob(
                os.path.join(self.inDir, sample_name, 'segmentation',
                             f"{segmentation_method}-{sample_name}", 
                             f"{segmentation_object}*.tif"))[0]

    return file_path


def napari_notification(msg):
    notification_ = Notification(msg, severity=NotificationSeverity.INFO)
    notification_manager.dispatch(notification_)


def read_markers(markers_filepath, markers_to_exclude, data):

    markers = pd.read_csv(markers_filepath, dtype={0: 'int16', 1: 'int16', 2: 'str'}, comment='#')
    if data is None:
        markers_to_include = [i for i in markers['marker_name'] if i not in markers_to_exclude]
    else:
        markers_to_include = [
            i for i in markers['marker_name'] if i not in markers_to_exclude if i in data.columns
        ]

    markers = markers[markers['marker_name'].isin(markers_to_include)]

    dna1 = markers['marker_name'][markers['channel_number'] == markers['channel_number'].min()][0]
    dna_moniker = str(re.search(r'[^\W\d]+', dna1).group())

    # abx channels
    abx_channels = [i for i in markers['marker_name'] if dna_moniker not in i]

    return markers, dna1, dna_moniker, abx_channels


def marker_channel_number(markers, marker_name):

    channel_number = markers.index[markers['marker_name'] == marker_name].item()

    return channel_number


def reorganize_dfcolumns(data, markers, cluster_dim):

    first_cols = (
        ['CellID'] +
        [i for i in markers['marker_name'] if
         i in data.columns] +
        [f'{i}_bool' for i in markers['marker_name'] if
         f'{i}_bool' in data.columns] +
        ['X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength',
         'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent',
         'Orientation', 'Sample', 'Condition', 'Replicate']
    )

    # (for BAF project)
    # first_cols = (
    #     ['CellID'] + [i for i in markers['marker_name'] if
    #                   i in data.columns] +
    #     ['X_centroid', 'Y_centroid', 'Area', 'CytArea', 'Solidity',
    #      'AreaSubstruct', 'MeanInsideSubstruct', 'Corenum', 'CoreCoord',
    #      'CoreFlag', 'Sample', 'Condition', 'Replicate']
    # )

    last_cols = [col for col in data.columns if col not in first_cols]

    data = data[first_cols + last_cols]

    return data


def single_channel_pyramid(tiff_path, channel):

    tiff = tifffile.TiffFile(tiff_path)

    if 'Faas' not in tiff.pages[0].software:

        if len(tiff.series[0].levels) > 1:

            pyramid = [zarr.open(s[channel].aszarr()) for s in tiff.series[0].levels]

            pyramid = [da.from_zarr(z) for z in pyramid]

            min_val = pyramid[-1].min()
            max_val = pyramid[-1].max()
            vmin, vmax = da.compute(min_val, max_val)

        else:

            img = tiff.pages[channel].asarray()

            pyramid = [img[::4**i, ::4**i] for i in range(4)]

            pyramid = [da.from_array(z) for z in pyramid]

            min_val = pyramid[-1].min()
            max_val = pyramid[-1].max()
            vmin, vmax = da.compute(min_val, max_val)

        return pyramid, vmin, vmax

    else:  # support legacy OME-TIFF format

        if len(tiff.series) > 1:

            pyramid = [zarr.open(s[channel].aszarr()) for s in tiff.series]

            pyramid = [da.from_zarr(z) for z in pyramid]

            min_val = pyramid[-1].min()
            max_val = pyramid[-1].max()
            vmin, vmax = da.compute(min_val, max_val)

        else:
            img = tiff.pages[channel].asarray()

            pyramid = [img[::4**i, ::4**i] for i in range(4)]

            pyramid = [da.from_array(z) for z in pyramid]

            min_val = pyramid[-1].min()
            max_val = pyramid[-1].max()
            vmin, vmax = da.compute(min_val, max_val)

        return pyramid, vmin, vmax


def matplotlib_warnings(fig):

    # suppress pyplot warning:
    # MatplotlibDeprecationWarning: Toggling axes navigation
    # from the keyboard is deprecated since 3.3 and will be
    # removed two minor releases later.
    fig.canvas.mpl_disconnect(
        fig.canvas.manager.key_press_handler_id)


def napari_warnings():

    # suppress warning:
    # vispy_camera.py:109: RuntimeWarning: divide by
    # zero encountered in true_divide
    # zoom = np.min(canvas_size / scale)
    np.seterr(divide='ignore')


# scatter point selection tool
class SelectFromCollection(object):
    """Select indices from a matplotlib collection using `LassoSelector`.

    Selected indices are saved in the `ind` attribute. This tool fades out the
    points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently
    alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        Axes to interact with.

    collection : :class:`matplotlib.collections.Collection` subclass
        Collection you want to select from.

    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to `alpha_other`.
    """

    def __init__(self, ax, collection, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other

        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))

        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


def read_dataframe(modules_list, start_module, outDir):

    # target_module_idx = modules_list.index(start_module) - 1
    # target_module_name = modules_list[target_module_idx]

    fileList = glob.glob(os.path.join(outDir, 'dataframe_archive/*.parquet'))
    parquet_name = fileList[0].split('/')[-1]
    logger.info(f'Reading: {parquet_name}')

    for filePath in fileList:
        df = pd.read_parquet(filePath)

    return df


def save_dataframe(df, outDir, moduleName):

    fileList = glob.glob(os.path.join(outDir, 'dataframe_archive/*.parquet'))

    for filePath in fileList:
        os.remove(filePath)

    df.to_parquet(
        os.path.join(outDir, f'dataframe_archive/{moduleName}.parquet'), index=True
    )


def categorical_cmap(numUniqueSamples, numCatagories, cmap='tab10', continuous=False):

    numSubcatagories = math.ceil(numUniqueSamples / numCatagories)

    if numCatagories > plt.get_cmap(cmap).N:
        raise ValueError('Too many categories for colormap.')
    if continuous:
        ccolors = plt.get_cmap(cmap)(np.linspace(0, 1, numCatagories))
    else:
        ccolors = plt.get_cmap(cmap)(np.arange(numCatagories, dtype=int))
        # rearrange hue order to taste
        cd = {
            'B': 0, 'O': 1, 'G': 2, 'R': 3, 'Pu': 4,
            'Br': 5, 'Pi': 6, 'Gr': 7, 'Y': 8, 'Cy': 9,
        }
        myorder = [
            cd['B'], cd['O'], cd['G'], cd['Pu'], cd['Y'],
            cd['R'], cd['Cy'], cd['Br'], cd['Gr'], cd['Pi']
        ]
        ccolors = [ccolors[i] for i in myorder]

        # use Okabe and Ito color-safe palette for first 6 colors
        # ccolors[0] = np.array([0.91, 0.29, 0.235]) #E84A3C
        # ccolors[1] = np.array([0.18, 0.16, 0.15]) #2E2926
        ccolors[0] = np.array([0.0, 0.447, 0.698, 1.0])  # blue
        ccolors[1] = np.array([0.902, 0.624, 0.0, 1.0])  # orange
        ccolors[2] = np.array([0.0, 0.620, 0.451, 1.0])  # bluish green
        ccolors[3] = np.array([0.8, 0.475, 0.655, 1.0])  # reddish purple
        ccolors[4] = np.array([0.941, 0.894, 0.259, 1.0])  # yellow
        ccolors[5] = np.array([0.835, 0.369, 0.0, 1.0])  # vermillion

    cols = np.zeros((numCatagories * numSubcatagories, 3))
    for i, c in enumerate(ccolors):
        chsv = colors.rgb_to_hsv(c[:3])
        arhsv = np.tile(chsv, numSubcatagories).reshape(numSubcatagories, 3)
        arhsv[:, 1] = np.linspace(chsv[1], 0.25, numSubcatagories)
        arhsv[:, 2] = np.linspace(chsv[2], 1, numSubcatagories)
        rgb = colors.hsv_to_rgb(arhsv)
        cols[i * numSubcatagories:(i + 1) * numSubcatagories, :] = rgb
    cmap = colors.ListedColormap(cols)

    # trim colors if necessary
    if len(cmap.colors) > numUniqueSamples:
        trim = len(cmap.colors) - numUniqueSamples
        cmap_colors = cmap.colors[:-trim]
        cmap = colors.ListedColormap(cmap_colors, name='from_list', N=None)

    return cmap


def cluster_expression(df, markers, cluster, num_proteins, clus_dim, norm_ax):

    if cluster != -1:

        df = df[df[f'cluster_{clus_dim}d'] != -1]

        cluster_means = (
            df[markers + [f'cluster_{clus_dim}d']]
            .groupby(f'cluster_{clus_dim}d').mean()
        )

        scaler = MinMaxScaler(feature_range=(0, 1), copy=True)

        if norm_ax == 'clusters':

            # rescale across clusters
            vectors = cluster_means.values.T
            scaled_rows = scaler.fit_transform(vectors.T)
            scaled_cluster_means = pd.DataFrame(
                data=scaled_rows, index=cluster_means.index,
                columns=cluster_means.columns
            )

            # isolate markers with highest expression values across clusters
            hi_markers = (
                scaled_cluster_means
                .loc[cluster]
                .sort_values(ascending=False)[:num_proteins]
                .index.tolist()
            )

            return norm_ax, hi_markers

        elif norm_ax == 'channels':

            # rescale across channels
            vectors = cluster_means.values
            scaled_rows = scaler.fit_transform(vectors.T).T
            scaled_cluster_means = pd.DataFrame(
                data=scaled_rows, index=cluster_means.index,
                columns=cluster_means.columns
            )

            # isolate markers with highest expression values across channels
            hi_markers = (
                scaled_cluster_means
                .loc[cluster]
                .sort_values(ascending=False)[:num_proteins]
                .index.tolist()
            )

            return norm_ax, hi_markers

    else:
        hi_markers = 'unclustered cells'

        return norm_ax, hi_markers


def gate_expression(pop, gate_dir):

    if os.path.exists(os.path.join(gate_dir, 'signatures.pkl')):
        f = open(os.path.join(gate_dir, 'signatures.pkl'), 'rb')
        signatures = pickle.load(f)
        hi_markers = [str(i) for i in signatures[pop] if i.negated is False]
        return hi_markers
    else:
        logger.info('Aborting; Cell classification dictionary does not exit in gating output.')
        sys.exit()


def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):

    # from statsmodels.stats.multitest

    '''
    pvalue correction for false discovery rate.

    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests.

    Parameters
    ----------
    pvals : array_like, 1d
        Set of p-values of the individual tests.
    alpha : float, optional
        Family-wise error rate. Defaults to ``0.05``.
    method : {'i', 'indep', 'p', 'poscorr', 'n', 'negcorr'}, optional
        Which method to use for FDR correction.
        ``{'i', 'indep', 'p', 'poscorr'}`` all refer to ``fdr_bh``
        (Benjamini/Hochberg for independent or positively
        correlated tests). ``{'n', 'negcorr'}`` both refer to ``fdr_by``
        (Benjamini/Yekutieli for general or negatively correlated tests).
        Defaults to ``'indep'``.
    is_sorted : bool, optional
        If False (default), the p_values will be sorted, but the corrected
        pvalues are in the original order. If True, then it assumed that the
        pvalues are already sorted in ascending order.

    Returns
    -------
    rejected : ndarray, bool
        True if a hypothesis is rejected, False if not
    pvalue-corrected : ndarray
        pvalues adjusted for multiple hypothesis testing to limit FDR

    Notes
    -----
    If there is prior information on the fraction of true hypothesis,
    then alpha should be set to ``alpha * m/m_0`` where m is the number of
    tests, given by the p-values, and m_0 is an estimate of the true
    hypothesis. (see Benjamini, Krieger and Yekuteli)

    The two-step method of Benjamini, Krieger and Yekutiel that estimates the
    number of false hypotheses will be available (soon).

    Both methods exposed via this function
    (Benjamini/Hochberg, Benjamini/Yekutieli) are also available in the
    function ``multipletests``, as ``method="fdr_bh"`` and
    ``method="fdr_by"``, respectively.

    See also
    --------
    multipletests

    '''

    def _ecdf(x):

        '''no frills empirical cdf used in fdrcorrection
        '''
        nobs = len(x)
        return np.arange(1, nobs + 1) / float(nobs)

    pvals = np.asarray(pvals)
    assert pvals.ndim == 1, "pvals must be 1-dimensional, of shape (n,)"

    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    else:
        pvals_sorted = pvals  # alias

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))  # corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
    # elif method in ['n', 'negcorr']:
    #     cm = np.sum(np.arange(len(pvals)))
    #     ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and negcorr implemented')
    reject = pvals_sorted <= ecdffactor * alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected
        reject_ = np.empty_like(reject)
        reject_[pvals_sortind] = reject
        return reject_, pvals_corrected_
    else:
        return reject, pvals_corrected


def triangulate_ellipse(corners, num_segments=100):

    # from Napari's GitHub page: napari/napari/layers/shapes/_shapes_utils.py

    """Determines the triangulation of a path. The resulting `offsets` can
    multiplied by a `width` scalar and be added to the resulting `centers`
    to generate the vertices of the triangles for the triangulation, i.e.
    `vertices = centers + width*offsets`. Using the `centers` and `offsets`
    representation thus allows for the computed triangulation to be
    independent of the line width.
    Parameters
    ----------
    corners : np.ndarray
        4xD array of four bounding corners of the ellipse. The ellipse will
        still be computed properly even if the rectangle determined by the
        corners is not axis aligned
    num_segments : int
        Integer determining the number of segments to use when triangulating
        the ellipse
    Returns
    -------
    vertices : np.ndarray
        Mx2 array coordinates of vertices for triangulating an ellipse.
        Includes the center vertex of the ellipse, followed by `num_segments`
        vertices around the boundary of the ellipse
    triangles : np.ndarray
        Px2 array of the indices of the vertices for the triangles of the
        triangulation. Has length given by `num_segments`
    """
    if not corners.shape[0] == 4:
        raise ValueError(
            trans._(
                "Data shape does not match expected `[4, D]` shape"
                "specifying corners for the ellipse",
                deferred=True,
            )
        )

    center = corners.mean(axis=0)
    adjusted = corners - center

    vec = adjusted[1] - adjusted[0]
    len_vec = np.linalg.norm(vec)
    if len_vec > 0:
        # rotate to be axis aligned
        norm_vec = vec / len_vec
        if corners.shape[1] == 2:
            transform = np.array(
                [[norm_vec[0], -norm_vec[1]], [norm_vec[1], norm_vec[0]]]
            )
        else:
            transform = np.array(
                [
                    [0, 0],
                    [norm_vec[0], -norm_vec[1]],
                    [norm_vec[1], norm_vec[0]],
                ]
            )
        adjusted = np.matmul(adjusted, transform)
    else:
        transform = np.eye(corners.shape[1])

    radii = abs(adjusted[0])
    vertices = np.zeros((num_segments + 1, 2), dtype=np.float32)
    theta = np.linspace(0, np.deg2rad(360), num_segments)
    vertices[1:, 0] = radii[0] * np.cos(theta)
    vertices[1:, 1] = radii[1] * np.sin(theta)

    if len_vec > 0:
        # rotate back
        vertices = np.matmul(vertices, transform.T)

    # Shift back to center
    vertices = vertices + center

    triangles = np.array([[0, i + 1, i + 2] for i in range(num_segments)])
    triangles[-1, 2] = 1

    return vertices, triangles


def artifact_detector_v3(pyramid, 
                      downscale=2, 
                      erosion_kernel_size=5,
                      lower_quantile_cutoff=0.2,
                      max_contrast=20,
                      h=None,
                      debug=False):
    if debug:
        fig, axes = plt.subplots(1, 6, figsize=(24, 4))
    kernel = disk(erosion_kernel_size)
    kernel_large = disk(erosion_kernel_size*4)
    im = rescale_intensity(pyramid[downscale].compute(), out_range=np.uint8)
    im_eroded = rank_min(im, kernel)
    lower_quantile = np.quantile(im_eroded.ravel(), lower_quantile_cutoff)
    im_eroded[im_eroded < lower_quantile] = lower_quantile
    im_transformed = rank_max(mean(im_eroded, kernel), kernel)
    local_contrast = gradient(im_transformed, kernel_large)
    
    if h is None:
        avg_ccmp_areas = []
        for i in range(max_contrast):
            num = np.sum(local_contrast > i)
            denom = np.max(label(local_contrast > i))
            if denom==0:
                break
            else: 
                avg_ccmp_areas.append(num/denom)
        h = np.argmin(avg_ccmp_areas)
    local_maxima_ = h_maxima(im_transformed, h=h, footprint=kernel_large) if h > 0 else\
                    local_maxima(im_transformed)
    local_maxima_labeled = label(local_maxima_)
    artifact_mask = np.zeros_like(im, dtype=np.int16)
    max_contrast = np.max(im_transformed) - np.min(im_transformed)
    
    num_local_maxima = np.max(local_maxima_labeled)
    seeds = np.empty((num_local_maxima, 2), dtype=int)
    optimal_tols = np.empty((num_local_maxima,), dtype=int)
    
    for seg_id in range(1, num_local_maxima+1):
        seed = np.argwhere(local_maxima_labeled==seg_id)[0]
        num_filled_pixels = [np.sum(flood(im_transformed, seed_point=tuple(seed), 
                                                tolerance=tol)) for tol in range(max_contrast+5)]
        delta_num = np.diff(num_filled_pixels)
        try: 
            optimal_tol = argrelextrema(delta_num, np.greater)[0][-2] 
            if optimal_tol >= 20:
                optimal_tol -= 2 #kind of a fudging factor here
        except:
            optimal_tol = max(0, np.argmax(delta_num > np.mean(delta_num))-1) #-1 is a conservative fudging factor here
        
        seeds[seg_id-1, :] = seed
        optimal_tols[seg_id-1] = optimal_tol
        artifact_mask += flood(im_transformed, seed_point=tuple(seed), 
                                                tolerance=optimal_tol)        

    if debug:
        axes[0].imshow(im, cmap='gray')
        axes[1].imshow(im_eroded, cmap='gray')
        axes[2].imshow(im_transformed, cmap='gray')
        axes[3].imshow(local_contrast, cmap='gray')
        axes[4].imshow(local_maxima_labeled, cmap='gray')
        axes[5].imshow(artifact_mask, cmap='gray')
    return artifact_mask, im_transformed, seeds, optimal_tols


def upscale(raw_im, target_im):
    return skimage.transform.resize(raw_im, target_im.shape,
                                     order=0, preserve_range=True, 
                                     anti_aliasing=False)

@dataclass
class ArtifactInfo():
    params: Dict
    mask: np.ndarray
    transformed: np.ndarray
    seeds: np.ndarray
    tols: np.ndarray
    features: pd.DataFrame = None
    artifact_layer: napari.layers.Image = None
    seed_layer: napari.layers.Points = None
        
    def update_mask(self, new_mask):
        self.mask = new_mask
        self.artifact_layer.data = upscale(self.mask, 
                                           self.artifact_layer.data)
        self.artifact_layer.refresh()

    def render_seeds(self, viewer, loaded_ims, layer_name, abx_channel):
        seeds = np.vstack(list(self.seeds.values()))
        self.seed_layer = viewer.add_points(seeds*(2**self.params['downscale']), 
                        name=layer_name[abx_channel+'_seeds'],
                        face_color=[1,0,0,1],
                        edge_color=[0,0,0,0],
                        size=int(max(*self.mask.shape) * (2**self.params['downscale']) / 100),
                        features={
                            'id': list(self.seeds.keys()),
                            'tol': self.tols
                        }, visible=False)
        self.features = self.seed_layer.features
        self.seed_layer.metadata['abx_channel'] = abx_channel
        self.seed_layer.metadata['downscale'] = self.params['downscale']

    def render_mask(self, viewer, loaded_ims, layer_name, abx_channel):
        grayscale = upscale(self.mask > 0, loaded_ims[abx_channel][0])
        self.artifact_layer = viewer.add_image(grayscale,
                                        name=layer_name[abx_channel+'_mask'], 
                                        opacity=0.5, visible=False,
                                        blending='additive')
        self.artifact_layer.metadata['abx_channel'] = abx_channel

    def render(self, viewer, loaded_ims, layer_name, abx_channel):
        self.render_mask(viewer, loaded_ims, layer_name, abx_channel)
        self.render_seeds(viewer, loaded_ims, layer_name, abx_channel)
        