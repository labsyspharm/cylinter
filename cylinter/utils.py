import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib import colors
from sklearn.preprocessing import MinMaxScaler
import os
import math
import subprocess
import yaml
import zarr
from hurry.filesize import size
import psutil
import gc

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


def read_markers(markers_filepath, mask_object, markers_to_exclude):

    markers = pd.read_csv(
        markers_filepath,
        dtype={0: 'int16', 1: 'int16', 2: 'str'},
        comment='#'
        )

    markers_to_include = [
        i for i in markers['marker_name']
        if i not in markers_to_exclude
        ]
    markers = markers[markers['marker_name'].isin(markers_to_include)]

    dna1 = markers['marker_name'][markers['channel_number'] == 1][0]
    dna_moniker = str(re.search(r'[^\W\d]+', dna1).group())

    # abx channels
    abx_channels = [
        f'{i}_{mask_object}' for i in markers['marker_name'] if
        dna_moniker not in i
        ]

    return markers, dna1, dna_moniker, abx_channels


def categorical_cmap(numUniqueSamples, numCatagories, cmap='tab10', continuous=False):

    numSubcatagories = math.ceil(numUniqueSamples/numCatagories)

    if numCatagories > plt.get_cmap(cmap).N:
        raise ValueError('Too many categories for colormap.')
    if continuous:
        ccolors = plt.get_cmap(cmap)(np.linspace(0, 1, numCatagories))
    else:
        ccolors = plt.get_cmap(cmap)(np.arange(numCatagories, dtype=int))
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


def cluster_expression(df, markers, cluster, num_proteins, across_or_within='across'):

    if cluster > -1:

        df = df[df['cluster'] > -1]

        cluster_means = df[markers + ['cluster']].groupby('cluster').mean()

        # rescale across clusters first
        # (this will make the clustermap results agree with the top
        # expressed markers shown for each cluster)
        min_max_scaler = MinMaxScaler()
        rescaled_vals = min_max_scaler.fit_transform(cluster_means.values)
        cluster_means_rescaled = pd.DataFrame(
            index=cluster_means.index,
            columns=cluster_means.columns,
            data=rescaled_vals
            )

        if across_or_within == 'across':
            max_cluster_ids = [
                cluster_means[i].idxmax() for
                i in cluster_means
                ]

            max_cluster_vals = [
                cluster_means[i].max() for
                i in cluster_means
                ]

            id_val_tuples = list(tuple(zip(max_cluster_ids, max_cluster_vals)))

            means_dict = dict(
                zip(markers, id_val_tuples)
                )

            means_dict = {
                k: v for k, v in
                means_dict.items() if v[0] == cluster
                }

            # sort means_dict by second element of 2-tuple
            # (i.e. mean expression value)
            total_markers = [
                k for k, v in sorted(
                    means_dict.items(),
                    key=lambda item: item[1][1])
                ]

            total_markers.reverse()

            hi_markers = total_markers[:num_proteins]

            return hi_markers

        elif across_or_within == 'within':

            cluster_means_rescaled = cluster_means_rescaled.loc[
                cluster].sort_values(ascending=False)

            hi_markers = list(cluster_means_rescaled.index[:num_proteins])

            return hi_markers

    else:
        hi_markers = []
        return hi_markers


def loadZarrs(df, outDir, markers_filepath, mask_object):

    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath,
        mask_object=mask_object,
        )

    zarrs_dir = os.path.join(outDir, 'zarrs')

    zs = {}
    for sample_name in df['Sample'].unique():

        if 'unmicst-' in sample_name:
            img_name = sample_name.split('unmicst-')[1]

        else:
            img_name = sample_name

        zs[f'{img_name}_{dna1}'] = (
            zarr.open(f'{zarrs_dir}/{img_name}_{dna1}.zarr', mode='r')
            )

        abx_channels = [
            i for i in markers['marker_name'] if
            dna_moniker not in i
            ]
        for ab in abx_channels:
            zs[f'{img_name}_{ab}'] = zarr.open(
                f'{zarrs_dir}/{img_name}_{ab}.zarr', mode='r'
                )

    return zs


def clearRAM(print_usage=False):

    gc.collect()

    if print_usage:
        process = psutil.Process(os.getpid())

        print(size(process.memory_info().rss))
