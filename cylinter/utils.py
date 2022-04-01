import os
import sys
import subprocess
import re
import gc
import csv
import glob
import math
import yaml
import zarr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib import colors
from sklearn.preprocessing import MinMaxScaler
from hurry.filesize import size
import psutil
import zarr
import dask.array as da
import tifffile

SUPPORTED_EXTENSIONS = ['.csv']


def dataset_files(root):
    """
    Return a list of all supported extension
    files in the specified directory.
    """
    total = os.listdir(root)
    pruned = [i for i in total if i.endswith(tuple(SUPPORTED_EXTENSIONS))]
    return pruned


def log_banner(log_function, msg):
    """Call log_function with a blank line, msg, and an underline."""
    log_function("")
    log_function(msg)
    log_function("-" * len(msg))


def log_multiline(log_function, msg):
    """Call log_function once for each line of msg."""
    for line in msg.split("\n"):
        log_function(line)


def dict_to_csv(dict, path):
    file = open(path, 'w')
    writer = csv.writer(file)
    for key, value in dict.items():
        writer.writerow([key, repr(value)])
    file.close()


def csv_to_dict(path):
    with open(path, 'r') as inp:
        reader = csv.reader(inp)
        dict = {rows[0]: rows[1] for rows in reader}
    return dict


def reorganize_dfcolumns(data, markers, cluster_dim):
    cols = (
        ['CellID'] + [i for i in markers['marker_name'] if
                      i in data.columns] +
        ['X_centroid', 'Y_centroid', 'Area', 'MajorAxisLength',
         'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent',
         'Orientation', 'Sample', 'Condition', 'Replicate']
         )

    if f'cluster_{cluster_dim}d' in data:
        cols = cols + [f'cluster_{cluster_dim}d']

    data = data[cols]
    return data


def single_channel_pyramid(tiff_path, channel):

    target_filepath = tiff_path
    tiff = tifffile.TiffFile(target_filepath, is_ome=False)

    pyramid = [
        zarr.open(s[channel].aszarr())
        for s in tiff.series[0].levels
        ]

    pyramid = [
        da.from_zarr(z)
        for z in pyramid
        ]

    return pyramid


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
    print(f'Reading: {parquet_name}')

    for filePath in fileList:
        df = pd.read_parquet(filePath)

    return df


def save_dataframe(df, outDir, moduleName):

    fileList = glob.glob(os.path.join(outDir, 'dataframe_archive/*.parquet'))

    for filePath in fileList:
        os.remove(filePath)

    df.to_parquet(
        os.path.join(outDir, f'dataframe_archive/{moduleName}.parquet'),
        index=True
        )


def read_markers(markers_filepath, markers_to_exclude, data):

    markers = pd.read_csv(
        markers_filepath,
        dtype={0: 'int16', 1: 'int16', 2: 'str'},
        comment='#'
        )
    if data is None:
        markers_to_include = [
            i for i in markers['marker_name']
            if i not in markers_to_exclude
            ]
    else:
        markers_to_include = [
            i for i in markers['marker_name']
            if i not in markers_to_exclude if i in data.columns
            ]

    markers = markers[markers['marker_name'].isin(markers_to_include)]

    dna1 = markers['marker_name'][markers['channel_number'] == 1][0]
    dna_moniker = str(re.search(r'[^\W\d]+', dna1).group())

    # abx channels
    abx_channels = [
        i for i in markers['marker_name'] if
        dna_moniker not in i
        ]

    return markers, dna1, dna_moniker, abx_channels


def marker_channel_number(markers, marker_name):

    channel_number = markers['channel_number'][
                markers['marker_name']
                == marker_name]

    return channel_number


def categorical_cmap(numUniqueSamples, numCatagories, cmap='tab10', continuous=False):

    numSubcatagories = math.ceil(numUniqueSamples/numCatagories)

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
    If there is prior information on the fraction of true hypothesis, then alpha
    should be set to ``alpha * m/m_0`` where m is the number of tests,
    given by the p-values, and m_0 is an estimate of the true hypothesis.
    (see Benjamini, Krieger and Yekuteli)

    The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
    of false hypotheses will be available (soon).

    Both methods exposed via this function (Benjamini/Hochberg, Benjamini/Yekutieli)
    are also available in the function ``multipletests``, as ``method="fdr_bh"`` and
    ``method="fdr_by"``, respectively.

    See also
    --------
    multipletests

    '''

    def _ecdf(x):

        '''no frills empirical cdf used in fdrcorrection
        '''
        nobs = len(x)
        return np.arange(1, nobs + 1)/float(nobs)

    pvals = np.asarray(pvals)
    assert pvals.ndim == 1, "pvals must be 1-dimensional, that is of shape (n,)"

    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    else:
        pvals_sorted = pvals  # alias

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
##    elif method in ['n', 'negcorr']:
##        cm = np.sum(np.arange(len(pvals)))
##        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and negcorr implemented')
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected
        reject_ = np.empty_like(reject)
        reject_[pvals_sortind] = reject
        return reject_, pvals_corrected_
    else:
        return reject, pvals_corrected


def clearRAM(print_usage=False):

    gc.collect()

    if print_usage:
        process = psutil.Process(os.getpid())

        print(size(process.memory_info().rss))


def open_file(filename):
    if sys.platform == "win32":
        os.startfile(filename)
    else:
        opener = "open" if sys.platform == "darwin" else "xdg-open"
        subprocess.run([opener, filename])
