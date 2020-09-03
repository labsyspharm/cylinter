import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
from matplotlib import colors
import os
import math
import subprocess
import yaml
import zarr


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


def save_dataframe(df, outDir, moduleName):
    print('Saving dataframe...')
    df.to_csv(os.path.join(outDir, f'dataframe_archive/{moduleName}.csv'))


def read_dataframe(outDir):
    df = pd.read_csv(
        os.path.join(outDir, 'current_dataframe.csv'), index_col=0
        )
    return df


def read_markers(markers_filepath):
    markers = pd.read_csv(
        markers_filepath,
        dtype={0: 'int16', 1: 'int16', 2: 'str'},
        comment='#'
        )
    dna1 = markers['marker_name'][markers['channel_number'] == 1][0]
    dna_moniker = str(re.search(r'[^\W\d]+', dna1).group())
    return markers, dna1, dna_moniker


def categorical_cmap(numUniqueSamples, numCatagories, cmap='tab10',
                     continuous=False):

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


def cluster_expression(df, markers, cluster, num_proteins):

    cluster_means = df[
        markers +
        ['cluster']].groupby('cluster').mean()

    cluster_means = cluster_means[
        cluster_means.index != -1
        ]

    marker_names = [i for i in cluster_means]

    max_cluster_ids = [
        cluster_means[i].idxmax() for
        i in cluster_means
        ]

    max_cluster_vals = [
        cluster_means[i].max() for
        i in cluster_means
        ]

    id_val_tuples = list(
        tuple(
            zip(max_cluster_ids, max_cluster_vals))
        )

    means_dict = dict(
        zip(marker_names, id_val_tuples)
        )

    means_dict = {
        k: v for k, v in
        means_dict.items() if v[0] == cluster
        }

    total_markers = [
        k for k, v in sorted(
            means_dict.items(),
            key=lambda item: item[1][1])
        ]

    total_markers.reverse()

    hi_markers = total_markers[:num_proteins]

    return hi_markers


def loadZarrs(df, outDir, markers_filepath):

    markers, dna1, dna_prefix = read_markers(
        markers_filepath=markers_filepath
        )

    zarrs_dir = os.path.join(outDir, 'zarrs')

    zs = {}
    for sample_name in df['Sample'].unique():

        if '-' in sample_name:
            img_name = sample_name.split('-')[1]

        else:
            img_name = sample_name

        zs[f'{img_name}_{dna1}'] = (
            zarr.open(f'{zarrs_dir}/{img_name}_{dna1}.zarr', mode='r')
            )

        abx_channels = [
            i for i in markers['marker_name'] if
            dna_prefix not in i
            ]
        for ab in abx_channels:
            zs[f'{img_name}_{ab}'] = zarr.open(
                f'{zarrs_dir}/{img_name}_{ab}.zarr', mode='r'
                )

    return zs
