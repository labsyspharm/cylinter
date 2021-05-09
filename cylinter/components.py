import logging
import functools

import os
import re
import yaml
import math
import pickle
import subprocess

import gc
import hdbscan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from skimage.color import gray2rgb
from skimage.filters import gaussian
from skimage.util.dtype import img_as_float
from skimage.util.dtype import img_as_uint

from matplotlib.lines import Line2D
from matplotlib.widgets import Slider, Button
from matplotlib.widgets import TextBox
from matplotlib.colors import ListedColormap
from matplotlib.path import Path
import matplotlib.patheffects as path_effects

import napari
from tifffile import imread
import zarr

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

from sklearn.preprocessing import MinMaxScaler
from umap import UMAP
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize as norm

from natsort import natsorted, order_by_index, index_natsorted
from numcodecs import Blosc
from datetime import datetime
from joblib import Memory
from scipy.stats import ttest_ind
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from decimal import Decimal
from bridson import poisson_disc_samples

from .utils import (
    dataset_files, log_banner, log_multiline,
    SelectFromCollection, read_dataframe, save_dataframe, read_markers,
    categorical_cmap, cluster_expression, loadZarrs, clearRAM, unmicst_version
    )

logger = logging.getLogger(__name__)

# map matplotlib color codes to the default seaborn palette
sns.set()
sns.set_color_codes()
_ = plt.plot([0, 1], color='r')
sns.set_color_codes()
_ = plt.plot([0, 2], color='b')
sns.set_color_codes()
_ = plt.plot([0, 3], color='g')
sns.set_color_codes()
_ = plt.plot([0, 4], color='m')
sns.set_color_codes()
_ = plt.plot([0, 5], color='y')
plt.close('all')

# Pipeline module order, to be filled in by the @module decorator.
pipeline_modules = []
pipeline_module_names = []


def module(func):
    """
    Annotation for pipeline module functions.

    This function adds the given function to the registry list. It also wraps
    the given function to log a pre/post-call banner.

    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger.info("=" * 70)
        logger.info("RUNNING MODULE: %s", func.__name__)
        result = func(*args, **kwargs)
        logger.info("=" * 70)
        logger.info("")
        return result
    pipeline_modules.append(wrapper)
    pipeline_module_names.append(wrapper.__name__)
    return wrapper


class QC(object):
    def __init__(self,

                 # config.yaml —
                 in_dir=None,
                 out_dir=None,
                 random_sample_size=None,
                 mask_object=None,
                 start_module=None,
                 sample_conditions=None,
                 sample_abbreviations=None,
                 sample_statuses=None,
                 sample_replicates=None,
                 samples_to_exclude=None,
                 markers_to_exclude=None,

                 # setContrast -
                 view_sample=None,

                 # selectROIs -
                 delint_mode=None,
                 show_ab_channels=None,

                 # crossCycleCorrelation -
                 cutoffAxis=None,
                 log_ratio_rnge=None,  # Can be None or (float, float)

                 # pruneOutliers -
                 hexbins=None,
                 hexbin_grid_size=None,

                 # PCA module —
                 channelExclusionsPCA=None,
                 samplesToRemovePCA=None,
                 dimensionPCA=None,
                 pointSize=None,
                 normalize=None,
                 labelPoints=None,
                 distanceCutoff=None,
                 samplesToSilhouette=None,

                 # clustering module —
                 embeddingAlgorithm1=None,
                 embeddingAlgorithm2=None,
                 channelExclusionsClustering1=None,
                 channelExclusionsClustering2=None,
                 samplesToRemoveClustering1=None,
                 samplesToRemoveClustering2=None,
                 normalizeTissueCounts1=None,
                 normalizeTissueCounts2=None,
                 fracForEmbedding1=None,
                 fracForEmbedding2=None,
                 dimensionEmbedding1=None,
                 dimensionEmbedding2=None,

                 perplexity1=None,
                 perplexity2=None,
                 earlyExaggeration1=None,
                 earlyExaggeration2=None,
                 learningRateTSNE1=None,
                 learningRateTSNE2=None,
                 metric1=None,
                 metric2=None,
                 random_state1=None,
                 random_state2=None,

                 nNeighbors1=None,
                 nNeighbors2=None,
                 learningRateUMAP1=None,
                 learningRateUMAP2=None,
                 minDist1=None,
                 minDist2=None,
                 repulsionStrength1=None,
                 repulsionStrength2=None,

                 # frequencyStats —
                 controlGroups=None,
                 denominatorCluster1=None,
                 denominatorCluster2=None,
                 FDRCorrection=None,

                 # curateThumbnails —
                 numThumbnails=None,
                 squareWindowDimension=None,

                 # dropClusters
                 clustersToDrop=None,

                 # clusterBoxplots —
                 bonferroniCorrection=None,

                 # spatialAnalysis —
                 cropDict=None,
                 spatialDict1=None,
                 spatialDict2=None,
                 radiusRange=None,
                 ):

        # assert(SOMETHING)  # placeholder for now

        self.in_dir = in_dir
        self.out_dir = out_dir
        self.random_sample_size = random_sample_size
        self.mask_object = mask_object
        self.start_module = start_module
        self.sample_conditions = sample_conditions
        self.sample_abbreviations = sample_abbreviations
        self.sample_statuses = sample_statuses
        self.sample_replicates = sample_replicates
        self.samples_to_exclude = samples_to_exclude
        self.markers_to_exclude = markers_to_exclude

        self.view_sample = view_sample

        self.delint_mode = delint_mode
        self.show_ab_channels = show_ab_channels

        self.cutoffAxis = cutoffAxis
        self.log_ratio_rnge = log_ratio_rnge

        self.hexbins = hexbins
        self.hexbin_grid_size = hexbin_grid_size

        self.channelExclusionsPCA = channelExclusionsPCA
        self.samplesToRemovePCA = samplesToRemovePCA
        self.dimensionPCA = dimensionPCA
        self.pointSize = pointSize
        self.normalize = normalize
        self.labelPoints = labelPoints
        self.distanceCutoff = distanceCutoff
        self.samplesToSilhouette = samplesToSilhouette

        self.embeddingAlgorithm1 = embeddingAlgorithm1
        self.embeddingAlgorithm2 = embeddingAlgorithm2
        self.channelExclusionsClustering1 = channelExclusionsClustering1
        self.channelExclusionsClustering2 = channelExclusionsClustering2
        self.samplesToRemoveClustering1 = samplesToRemoveClustering1
        self.samplesToRemoveClustering2 = samplesToRemoveClustering2
        self.normalizeTissueCounts1 = normalizeTissueCounts1
        self.normalizeTissueCounts2 = normalizeTissueCounts2
        self.fracForEmbedding1 = fracForEmbedding1
        self.fracForEmbedding2 = fracForEmbedding2
        self.dimensionEmbedding1 = dimensionEmbedding1
        self.dimensionEmbedding2 = dimensionEmbedding2

        self.perplexity1 = perplexity1
        self.perplexity2 = perplexity2
        self.earlyExaggeration1 = earlyExaggeration1
        self.earlyExaggeration2 = earlyExaggeration2
        self.learningRateTSNE1 = learningRateTSNE1
        self.learningRateTSNE2 = learningRateTSNE2
        self.metric1 = metric1
        self.metric2 = metric2
        self.random_state1 = random_state1
        self.random_state2 = random_state2

        self.nNeighbors1 = nNeighbors1
        self.nNeighbors2 = nNeighbors2
        self.learningRateUMAP1 = learningRateUMAP1
        self.learningRateUMAP2 = learningRateUMAP2
        self.minDist1 = minDist1
        self.minDist2 = minDist2
        self.repulsionStrength1 = repulsionStrength1
        self.repulsionStrength2 = repulsionStrength2

        self.controlGroups = controlGroups
        self.denominatorCluster1 = denominatorCluster1
        self.denominatorCluster2 = denominatorCluster2
        self.FDRCorrection = FDRCorrection

        self.numThumbnails = numThumbnails
        self.squareWindowDimension = squareWindowDimension

        self.clustersToDrop = clustersToDrop

        self.bonferroniCorrection = bonferroniCorrection

        self.cropDict = cropDict
        self.spatialDict1 = spatialDict1
        self.spatialDict2 = spatialDict2
        self.radiusRange = radiusRange

    @module
    def getSingleCellData(data, self, args):

        files = dataset_files(f'{self.in_dir}/csv')

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        df_list = []
        raw_sample_names_dict = {}
        channel_setlist = []
        for file in files:
            if not file.startswith('.'):
                raw_sample_name = file.split('.', -1)[0]

                if raw_sample_name.startswith('unmicst-'):
                    sample_name = raw_sample_name.split('-', 1)[1]
                    raw_sample_names_dict[sample_name] = raw_sample_name
                else:
                    sample_name = raw_sample_name
                    raw_sample_names_dict[sample_name] = raw_sample_name

                # disregard samples specified in
                # "samples_to_exclude" config param
                if raw_sample_name not in self.samples_to_exclude:

                    print(
                        f'Importing single-cell data for sample {sample_name}.'
                        )

                    csv = pd.read_csv(
                        os.path.join(f'{self.in_dir}/csv', file)
                        )

                    # drop markers specified in
                    # "markers_to_exclude" config param
                    csv.drop(
                        columns=[
                            f'{i}_{self.mask_object}' for i
                            in self.markers_to_exclude],
                        axis=1,
                        inplace=True
                        )

                    csv['Sample'] = sample_name

                    channel_setlist.append(set(csv.columns))

                    # append dataframe to list
                    df_list.append(csv)
                else:
                    print(
                        f'Censoring single-cell data for sample {sample_name}.'
                        )
        print()

        # stack dataframes row-wise
        df = pd.concat(df_list, axis=0)
        del df_list

        # add condition column
        df['Condition'] = [
            self.sample_conditions[s] for s in
            [raw_sample_names_dict[i] for i in df['Sample']]
            ]

        # add replicate column
        df['Replicate'] = [
            self.sample_replicates[s] for s in
            [raw_sample_names_dict[i] for i in df['Sample']]
            ]

        # organize columns
        cols = (
            ['CellID', 'Sample', 'Condition', 'Replicate', 'X_centroid',
             'Y_centroid', 'column_centroid', 'row_centroid', 'Area',
             'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity',
             'Extent', 'Orientation'] +
            [f'{i}_{self.mask_object}' for i in markers['marker_name']]
             )
        df = df[cols]

        # ONLY SELECT CHANNELS COMMON BETWEEN ALL SAMPLES
        channels_set = list(set.intersection(*channel_setlist))
        channels_set.extend(['Condition', 'Replicate'])

        print(f'{len(df.columns)} total features.')
        print(f'{len(channels_set)} features shared among all samples.')
        print()

        before = set(df.columns)
        after = set(channels_set)

        if len(before.difference(after)) == 0:
            pass
        else:
            markers_to_drop = list(before.difference(after))
            print(
                f'Features {markers_to_drop} were not probed for in all' +
                ' samples and will be dropped from the analysis.'
                )
        df = df[channels_set].copy()

        # handle data subsetting
        df = df.sample(frac=self.random_sample_size, random_state=1)
        df.sort_values(by=['Sample', 'CellID'], inplace=True)

        # assign global index
        df.reset_index(drop=True, inplace=True)
        print()

        return df

    @module
    def selectROIs(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        selection_dir = os.path.join(self.out_dir, 'ROIs')
        if not os.path.exists(selection_dir):
            os.makedirs(selection_dir)

        if os.path.exists(os.path.join(selection_dir, 'polygon_dict.pkl')):
            f = open(os.path.join(selection_dir, 'polygon_dict.pkl'), 'rb')
            polygon_dict = pickle.load(f)
            completed_samples = set(polygon_dict.keys())
            total_samples = set(data['Sample'].unique())
            samples_to_draw = total_samples.difference(completed_samples)
            print(f'Samples to draw: {len(samples_to_draw)}')

        else:
            samples_to_draw = data['Sample'].unique()
            print(f'Samples requiring ROI selection: {len(samples_to_draw)}')
            print()
            polygon_dict = {}

        if (
          (len(samples_to_draw) > 0) or
          (len([name for name in os.listdir(selection_dir) if
           name.endswith('.txt')]) < len(data['Sample'].unique()))):
            for sample_name in natsorted(samples_to_draw):

                dna = imread(
                    f'{self.in_dir}/tif/{sample_name}.*tif', key=0
                    )

                polygons = []
                with napari.gui_qt():
                    viewer = napari.view_image(
                        dna, rgb=False, blending='additive',
                        colormap='gray', visible=True,
                        name=f'{dna1}:{sample_name}'
                        )

                    if self.show_ab_channels:
                        for ch in abx_channels:
                            ch = ch.split(f'_{self.mask_object}')[0]
                            channel_number = markers['channel_number'][
                                        markers['marker_name'] == ch]

                            # read antibody image
                            img = imread(
                                f'{self.in_dir}/tif/' +
                                f'{sample_name}.*tif',
                                key=(channel_number-1)
                                )

                            viewer.add_image(
                                img, rgb=False, blending='additive',
                                colormap='green', visible=False,
                                name=ch
                                )

                    if len(polygons) == 0:
                        selection_layer = viewer.add_shapes(
                            shape_type='polygon',
                            ndim=2,
                            face_color=[1.0, 1.0, 1.0, 0.2],
                            edge_color=[0.0, 0.66, 1.0, 1.0],
                            edge_width=10.0,
                            name='ROI(s)'
                            )
                    else:
                        selection_layer = viewer.add_shapes(
                            polygons,
                            shape_type='polygon',
                            ndim=2,
                            face_color=[1.0, 1.0, 1.0, 0.2],
                            edge_color=[0.0, 0.66, 1.0, 1.0],
                            edge_width=10.0,
                            name='ROI(s)'
                            )

                    # @viewer.mouse_drag_callbacks.append
                    # def get_cell_indices(viewer, event):
                    #
                    #     # on mouse press
                    #     yield
                    #
                    #     # on mouse move
                    #     while event.type == 'mouse_move':
                    #         yield
                    #
                    #     # on mouse release
                    #     selection_layer = viewer.layers['ROI_1']
                    #     yield

                for roi in selection_layer.data:
                    # only process new polygon points
                    if not next(
                        (True for elem in polygons if
                         np.array_equal(elem, roi)), False):

                        polygons.append(roi)

                polygon_dict[sample_name] = polygons

                os.chdir(selection_dir)
                f = open(os.path.join(selection_dir, 'polygon_dict.pkl'), 'wb')
                pickle.dump(polygon_dict, f)
                f.close()

            os.chdir(selection_dir)
            samples_for_cell_selection = set(
                data['Sample'].unique()).difference(
                    set([i.split('.txt')[0] for i in os.listdir()
                         if i.endswith('.txt')]))

            for sample_name, group in natsorted(data.groupby('Sample')):
                if sample_name in samples_for_cell_selection:

                    print(f'Selecting cells for sample {sample_name}')

                    sample_data = group[
                        ['X_centroid', 'Y_centroid', 'CellID']].astype(int)
                    sample_data['tuple'] = list(
                        zip(sample_data['Y_centroid'],
                            sample_data['X_centroid'])
                        )

                    dna = imread(
                        f'{self.in_dir}/tif/{sample_name}.*tif', key=0
                        )

                    columns, rows = np.meshgrid(
                        np.arange(dna.shape[1]),
                        np.arange(dna.shape[0])
                        )
                    columns, rows = columns.flatten(), rows.flatten()
                    pixel_coords = np.vstack((rows, columns)).T

                    cell_coords = set(
                        [tuple(i) for i in np.array(
                            sample_data[['Y_centroid', 'X_centroid']])]
                        )

                    clearRAM(print_usage=False)
                    del columns, rows
                    clearRAM(print_usage=False)

                    cell_ids = set()
                    mask_coords = set()
                    if polygon_dict[sample_name]:
                        for poly in polygon_dict[sample_name]:
                            selection_verts = np.round(poly).astype(int)
                            polygon = Path(selection_verts)
                            # print('C')
                            grid = polygon.contains_points(pixel_coords)
                            mask = grid.reshape(
                                dna.shape[0], dna.shape[1])
                            # print('D')
                            mask_coords.update(
                                [tuple(i) for i in np.argwhere(mask)]
                                )
                            # print('E')
                        clearRAM(print_usage=False)
                        del grid, mask, dna, pixel_coords
                        clearRAM(print_usage=False)

                    inter = mask_coords.intersection(cell_coords)

                    if self.delint_mode:
                        cell_ids.update(
                            [i[1]['CellID'] for i in sample_data.iterrows() if
                             i[1]['tuple'] in inter]
                             )
                    else:
                        if polygon_dict[sample_name]:
                            cell_ids.update(
                                [i[1]['CellID'] for i in sample_data.iterrows()
                                 if i[1]['tuple'] not in inter]
                                 )  # take all cells if no ROIs drawn

                    clearRAM(print_usage=False)
                    del sample_data, inter, cell_coords
                    clearRAM(print_usage=False)

                    os.chdir(selection_dir)
                    f = open(f'{sample_name}.txt', 'wb')
                    pickle.dump(cell_ids, f)
                    f.close()
            print()

            # drop cell IDs from full dataframe
            os.chdir(selection_dir)
            for file in os.listdir(selection_dir):
                if file.endswith('.txt'):
                    file_name = file.split('.txt')[0]
                    print(f'Dropping cells from sample {file_name}')
                    f = open(file, 'rb')
                    cell_ids = pickle.load(f)
                    cell_ids = set(cell_ids)
                    global_idxs = data[
                        (data['Sample'] == file_name)
                        & (data['CellID'].isin(cell_ids))
                        ].index
                    data.drop(global_idxs, inplace=True)
            print()

            return data

        else:

            # drop cell IDs from full dataframe
            os.chdir(selection_dir)
            for file in os.listdir(selection_dir):
                if file.endswith('.txt'):
                    file_name = file.split('.txt')[0]
                    print(f'Dropping cells from {file_name}')
                    f = open(file, 'rb')
                    cell_ids = pickle.load(f)
                    cell_ids = set(cell_ids)
                    global_idxs = data[
                        (data['Sample'] == file_name)
                        & (data['CellID'].isin(cell_ids))
                        ].index
                    data.drop(global_idxs, inplace=True)
            print()

            return data

    @module
    def dnaIntensityCutoff(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        bins = 100

        histtype = 'stepfilled'

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        n, bins, patches = plt.hist(
            data[f'{dna1}_{self.mask_object}'], bins=bins,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype=histtype,
            range=None, label='before'
            )

        plt.title('mean DNA intensity')
        plt.ylabel('count')

        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        rnge = [bins.min(), bins.max()]

        sLower = Slider(
            axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sLower.label.set_color('b')

        sUpper = Slider(
            axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sUpper.label.set_color('r')

        def update(val):
            [i.remove() for i in ax.get_lines()]
            lowerCutoff = sLower.val
            upperCutoff = sUpper.val
            blueLine = ax.axvline(
                x=lowerCutoff, c='b', linewidth=2.5)
            redLine = ax.axvline(
                x=upperCutoff, c='r', linewidth=2.5)
            return lowerCutoff, upperCutoff

        sLower.on_changed(update)
        sUpper.on_changed(update)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sLower.reset()
            sUpper.reset()
        button.on_clicked(reset)

        def submit(text):
            lowerCutoff, upperCutoff = update(val=None)

            # apply lower and upper cutoffs
            df_test = data[
                (data[f'{dna1}_{self.mask_object}'] > lowerCutoff) &
                (data[f'{dna1}_{self.mask_object}'] < upperCutoff)
                ]

            if text in data['Sample'].unique():

                dna = imread(
                    f'{self.in_dir}/tif/{text}.*tif',
                    key=0
                    )

                centroids = df_test[
                    ['Y_centroid', 'X_centroid']][df_test['Sample'] == text]

                dna_intensity = np.array(df_test[
                    [f'{dna1}_{self.mask_object}']][df_test['Sample'] == text])

                point_properties = {
                    'dna_intensity': dna_intensity
                    }

                with napari.gui_qt():
                    viewer = napari.view_image(
                        dna, rgb=False, name=dna1
                        )

                    viewer.add_points(
                        centroids, name='DNA intensity',
                        properties=point_properties,
                        face_color='dna_intensity', face_colormap='viridis',
                        edge_color='k', edge_width=0.0, size=4.0
                        )

        axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
        text_box = TextBox(
            axbox, 'evaluation sample name', initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.on_submit(submit)
        plt.show(block=True)

        lowerCutoff, upperCutoff = update(val=None)

        fig, ax = plt.subplots()
        # plot DNA intensity histogram BEFORE filtering
        plt.hist(
            data[f'{dna1}_{self.mask_object}'], bins=bins,
            density=False, color='b', ec='none',
            alpha=0.5, histtype=histtype,
            range=None, label='before'
            )

        if lowerCutoff == upperCutoff:
            lowerCutoff = data[f'{dna1}_{self.mask_object}'].min()
            upperCutoff = data[f'{dna1}_{self.mask_object}'].max()

        # apply lower and upper cutoffs
        df = data[
            (data[f'{dna1}_{self.mask_object}'] > lowerCutoff) &
            (data[f'{dna1}_{self.mask_object}'] < upperCutoff)
            ]

        # plot DNA intensity histogram AFTER filtering
        plt.hist(
            df[f'{dna1}_{self.mask_object}'], bins=bins, color='r', ec='none',
            alpha=0.5, histtype=histtype, range=None, label='after')
        plt.xlabel('mean DNA intensity')
        plt.ylabel('count')

        legend_elements = []
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='excluded data',
                   markerfacecolor='b', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='included data',
                   markerfacecolor='r', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        plt.legend(
            handles=legend_elements, prop={'size': 10},
            loc='best'
            )

        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, 'histogram_dna.pdf'))
        plt.close()
        print()

        return df

    @module
    def dnaAreaCutoff(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        bins = 100

        histtype = 'stepfilled'

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        n, bins, patches = plt.hist(
            data['Area'], bins=bins,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype=histtype,
            range=None, label='before'
            )

        plt.title('mean DNA area')
        plt.ylabel('count')

        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        rnge = [bins.min(), bins.max()]

        sLower = Slider(
            axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sLower.label.set_color('b')

        sUpper = Slider(
            axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sUpper.label.set_color('r')

        def update(val):
            [i.remove() for i in ax.get_lines()]
            lowerCutoff = sLower.val
            upperCutoff = sUpper.val
            blueLine = ax.axvline(
                x=lowerCutoff, c='b', linewidth=2.5)
            redLine = ax.axvline(
                x=upperCutoff, c='r', linewidth=2.5)
            return lowerCutoff, upperCutoff

        sLower.on_changed(update)
        sUpper.on_changed(update)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sLower.reset()
            sUpper.reset()
        button.on_clicked(reset)

        def submit(text):
            lowerCutoff, upperCutoff = update(val=None)

            # apply lower and upper cutoffs
            df_test = data[
                (data['Area'] > lowerCutoff) &
                (data['Area'] < upperCutoff)
                ]

            if text in data['Sample'].unique():

                dna = imread(
                    f'{self.in_dir}/tif/{text}.*tif',
                    key=0
                    )

                # read antibody image
                seg = imread(
                    f'{self.in_dir}/seg/{text}.*tif',
                    key=0
                    )

                centroids = df_test[
                    ['Y_centroid', 'X_centroid']][df_test['Sample'] == text]

                dna_area = np.array(
                    df_test[['Area']][df_test['Sample'] == text]
                    )

                point_properties = {
                    'dna_area': dna_area
                    }

                with napari.gui_qt():
                    viewer = napari.view_image(
                        dna, rgb=False, name=dna1
                        )

                    viewer.add_image(
                        seg, rgb=False, blending='additive',
                        colormap='green', visible=False,
                        name='segmentation'
                        )

                    viewer.add_points(
                        centroids, name='DNA area',
                        properties=point_properties,
                        face_color='dna_area', face_colormap='viridis',
                        edge_color='k', edge_width=0.0, size=4.0
                        )

        axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
        text_box = TextBox(
            axbox, 'evaluation sample name', initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.on_submit(submit)
        plt.show(block=True)

        lowerCutoff, upperCutoff = update(val=None)

        fig, ax = plt.subplots()
        # plot DNA area histogram BEFORE filtering
        plt.hist(
            data['Area'], bins=bins,
            density=False, color='b', ec='none',
            alpha=0.5, histtype=histtype,
            range=None, label='before'
            )

        if lowerCutoff == upperCutoff:
            lowerCutoff = data['Area'].min()
            upperCutoff = data['Area'].max()

        # apply lower and upper cutoffs
        df = data[
            (data['Area'] > lowerCutoff) &
            (data['Area'] < upperCutoff)
            ]

        # plot DNA area histogram AFTER filtering
        plt.hist(
            df['Area'], bins=bins, color='r', ec='none',
            alpha=0.5, histtype=histtype, range=None, label='after')
        plt.xlabel('mean DNA area')
        plt.ylabel('count')

        legend_elements = []
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='excluded data',
                   markerfacecolor='b', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='included data',
                   markerfacecolor='r', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        plt.legend(
            handles=legend_elements, prop={'size': 10},
            loc='best'
            )

        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, 'histogram_area.pdf'))
        plt.close()

        # save images
        lasso_dir = os.path.join(self.out_dir, 'lassos')
        if not os.path.exists(lasso_dir):
            os.mkdir(lasso_dir)

        for sample_name in natsorted(df['Sample'].unique()):

            print(f'Plotting ROI selections for {sample_name}.')

            dna = imread(f'{self.in_dir}/tif/{sample_name}.*tif', key=0)

            fig, ax = plt.subplots()
            ax.imshow(dna, cmap='gray')
            ax.grid(False)
            coords = df[['X_centroid', 'Y_centroid', 'Area']][
                df['Sample'] == sample_name]
            sp = ax.scatter(
                coords['X_centroid'], coords['Y_centroid'],
                s=0.35, lw=0.0,
                c=coords['Area'], cmap='viridis'
                )
            plt.title(
                f'Sample {sample_name}. ' +
                'Selected cells colored by nuclear area')
            plt.colorbar(sp)
            plt.savefig(
                os.path.join(
                    lasso_dir, f'{sample_name}.png'), dpi=1000)
            plt.close('all')
        print()

        return df

    @module
    def crossCycleCorrelation(data, self, args):

        cycles_dir = os.path.join(self.out_dir, 'cycles')
        if not os.path.exists(cycles_dir):
            os.makedirs(cycles_dir)

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        dna_cycles = natsorted(
            data.columns[data.columns.str.contains(dna_moniker)]
            )

        # log(cycle 1/n) ratios
        ratios = pd.DataFrame(
            [np.log10(
                (data[f'{dna1}_{self.mask_object}'] + 0.00001) /
                (data[i] + 0.00001)) for i in dna_cycles]).T

        list1 = [i for i in ratios.columns if i.startswith('Unnamed')]
        list2 = [
            f'{dna1}/{i}' for i in
            [j.split(f'_{self.mask_object}')[0] for j in dna_cycles[1:]]
            ]

        ratio_columns = dict(zip(list1, list2))
        ratio_columns[f'{dna1}_{self.mask_object}'] = f'{dna1}/{dna1}'
        ratios.rename(columns=ratio_columns, inplace=True)
        ratios['sample'] = data['Sample']

        # log(cycle n/n+1) ratios
        # ratios = pd.DataFrame(
        #     [np.log10(
        #         (data[i] + 0.00001) /
        #         (data[dna_moniker
        #          + str(int([re.findall(r'(\w+?)(\d+)', i)[0]][0][1]) + 1)
        #          + '_' + self.mask_object] + 0.00001)
        #          )
        #         for i in
        #         natsorted(
        #             data.columns[data.columns.str.contains(dna_moniker)])[0:-1]]).T
        # list1 = [i for i in ratios.columns]
        # list2 = [f'{i}/{i+1}' for i in range(1, len(list1)+1)]
        # ratio_columns = dict(zip(list1, list2))
        # ratios.rename(columns=ratio_columns, inplace=True)
        # ratios['sample'] = data['Sample']

        ratios_melt = (
            ratios
            .reset_index()
            .melt(id_vars=['sample', 'index'], var_name='cycle',
                  value_name='log10(ratio)')
            )

        ratios_melt['sample'] = pd.Categorical(
            ratios_melt['sample'], ordered=True,
            categories=natsorted(
                ratios_melt['sample'].unique()))

        ratios_melt['cycle'] = pd.Categorical(
            ratios_melt['cycle'], ordered=True,
            categories=natsorted(
                ratios_melt['cycle'].unique()))

        ratios_melt = (
            ratios_melt.sort_values(['sample', 'cycle', 'index'])
            )

        if not os.path.exists(
          os.path.join(cycles_dir, 'cycle_correlation(logRatio).pdf')):

            sns.set(font_scale=0.5)
            sns.set_style('whitegrid')

            g = sns.FacetGrid(
                ratios_melt, row='sample',
                col='cycle', sharey=False
                )

            g = g.map(
                plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
                ec='none', range=self.log_ratio_rnge, bins=200, density=True
                )

            plt.savefig(
                os.path.join(cycles_dir, 'cycle_correlation(logRatio).pdf'))
            plt.close()

            subprocess.call(
                ['open', '-a', 'Preview', os.path.join(
                    cycles_dir, 'cycle_correlation(logRatio).pdf')])

        else:
            subprocess.call(
                ['open', '-a', 'Preview', os.path.join(
                    cycles_dir, 'cycle_correlation(logRatio).pdf')])

        if self.cutoffAxis == 'y':
            def submit(text):

                count_cutoff = float(text.split(', ')[0])
                os.chdir(cycles_dir)
                f = open('count_cutoff.pkl', 'wb')
                pickle.dump(count_cutoff, f)
                f.close()

                sample_to_inspect = str(text.split(', ')[1])

                if sample_to_inspect in data['Sample'].unique():
                    sample_df = data[data['Sample'] == sample_to_inspect]
                    sample_centroids = sample_df[['Y_centroid', 'X_centroid']]

                    with napari.gui_qt():

                        # add dna images
                        for e, i in enumerate(sorted(markers.index[
                          markers['marker_name'].str.contains(dna_moniker)],
                          reverse=True)):

                            dna = imread(
                                f'{self.in_dir}/tif/' +
                                f'{sample_to_inspect}.*tif',
                                key=i
                                )

                            name = markers.loc[i]['marker_name']
                            cycle_num = markers.loc[i]['cycle_number']

                            if name == dna1:
                                visible = True
                            else:
                                visible = False

                            if e == 0:
                                viewer = napari.view_image(
                                    dna, rgb=False, blending='opaque',
                                    visible=visible,
                                    name=name
                                    )
                            else:
                                viewer.add_image(
                                    dna, rgb=False, blending='opaque',
                                    visible=visible,
                                    name=name
                                    )

                            # log(cycle n/n+1) ratios
                            if cycle_num < markers['cycle_number'].max():
                                sample_ratios = np.log10(
                                    (sample_df[
                                        f'{name}_{self.mask_object}']
                                        + 0.00001)
                                    /
                                    (sample_df[dna_moniker
                                     + str(int(
                                        [re.findall(r'(\w+?)(\d+)', name)[0]
                                         ][0][1]
                                        ) + 1) + '_' + self.mask_object]
                                        + 0.00001)
                                    )

                            # log(cycle 1/n) ratios
                            # if cycle_num != 1:
                            #     sample_ratios = np.log10(
                            #         (sample_df[
                            #             f'{dna1}_{self.mask_object}']
                            #             + 0.00001)
                            #         /
                            #         (sample_df[
                            #             f'{name}_{self.mask_object}']
                            #             + 0.00001)
                            #         )

                                sample_ratios = np.clip(
                                    sample_ratios,
                                    a_min=np.percentile(sample_ratios, 1.0),
                                    a_max=np.percentile(sample_ratios, 99.0)
                                    )

                                point_properties = {
                                    'face_color': sample_ratios
                                    }

                                viewer.add_points(
                                    sample_centroids,
                                    name=f'log({cycle_num}/{cycle_num + 1})',
                                    visible=False,
                                    properties=point_properties,
                                    face_color='face_color',
                                    face_colormap='PiYG',
                                    edge_color='k', edge_width=0.0, size=7.0
                                    )

                        # rearrange viwer layers (DNA images first)
                        layer_names = [str(i) for i in viewer.layers]

                        current_order = tuple(
                            [layer_names.index(i) for i in layer_names]
                            )

                        dna_idxs = [
                            layer_names.index(i) for i in layer_names
                            if dna_moniker in i
                            ]
                        log_idxs = [
                            layer_names.index(i) for i in layer_names
                            if 'log(' in i
                            ]

                        target_order = tuple(log_idxs + dna_idxs)

                        viewer.layers[current_order] = viewer.layers[
                            target_order
                            ]

                sns.set(font_scale=0.5)
                sns.set_style('whitegrid')

                g = sns.FacetGrid(
                    ratios_melt, row='sample',
                    col='cycle', sharey=False
                    )
                g = g.map(
                    plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
                    ec='none', range=self.log_ratio_rnge,
                    bins=200, density=True
                    )

                for ax in g.axes.ravel():
                    ax.axhline(y=count_cutoff, c='k', linewidth=0.5)

                plt.savefig(
                    os.path.join(
                        cycles_dir, 'cycle_correlation(logRatio).pdf')
                        )

                subprocess.call(
                    ['open', '-a', 'Preview',
                     os.path.join(
                        cycles_dir, 'cycle_correlation(logRatio).pdf')]
                     )

                plt.show(block=False)
                plt.close()

            plt.rcParams['figure.figsize'] = (6, 2)
            axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
            text_box = TextBox(
                axbox, 'countCutoff, sampleName', initial='',
                color='0.95',
                hovercolor='1.0',
                label_pad=0.05
                )
            text_box.label.set_size(12)

            text_box.on_submit(submit)

            plt.show(block=True)

            if os.path.exists(
              os.path.join(cycles_dir, 'count_cutoff.pkl')):
                os.chdir(cycles_dir)
                pickle_in = open('count_cutoff.pkl', 'rb')
                count_cutoff = pickle.load(pickle_in)
            else:
                count_cutoff = 0.0

                sns.set(font_scale=0.5)
                sns.set_style('whitegrid')

                g = sns.FacetGrid(
                    ratios_melt, row='sample',
                    col='cycle', sharey=False
                    )
                g = g.map(
                    plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
                    ec='none', range=self.log_ratio_rnge,
                    bins=200, density=True
                    )

                for ax in g.axes.ravel():
                    ax.axhline(y=count_cutoff, c='k', linewidth=0.5)

                plt.savefig(
                    os.path.join(
                        cycles_dir, 'cycle_correlation(logRatio).pdf')
                        )

            # initialize a list to append indices to drop
            indices_to_drop = set()

            for name, group in ratios_melt.groupby(['sample']):
                for cycle_ratio in group['cycle'].unique():
                    if not (
                      cycle_ratio.split('/')[0]
                      == cycle_ratio.split('/')[1]):

                        # isolate ratio data
                        cycle_data = group[group['cycle'] == cycle_ratio]

                        # get histogram elements
                        sns.set_style('whitegrid')
                        fig, ax = plt.subplots()
                        counts, bins, patches = plt.hist(
                            cycle_data['log10(ratio)'],
                            color='r', histtype='stepfilled',
                            ec='none', range=None,
                            bins=200, density=True
                            )
                        plt.close('all')

                        # plot histogram of log(ratios)
                        # with a horizontal line at cutoff point
                        # ax.axhline(y=count_cutoff, c='k', linewidth=0.5)
                        # plt.title(sample + '_' + col_name)
                        # plt.show(block=True)

                        # get bin values (i.e. ratios) where cell
                        # counts are > than count_cutoff
                        count_indices = np.where(counts > count_cutoff)
                        bin_values = [bins[i] for i in count_indices[0]]

                        if len(bin_values) > 1:
                            min_bin_val = min(bin_values)
                            max_bin_val = max(bin_values)

                            # get indices in log(ratio) series outside
                            # min_bin_val and max_bin_val
                            idxs = list(
                                cycle_data['index'][
                                    (cycle_data['log10(ratio)']
                                     < min_bin_val) |
                                    (cycle_data['log10(ratio)']
                                     > max_bin_val)]
                                    )

                            # append indices of uncorrelated
                            # log(ratios) to idx_list
                            indices_to_drop.update(set(idxs))

            # filter dataframe by selecting indices NOT in the
            # indices_to_drop list
            df = data.loc[~data.index.isin(indices_to_drop)]

        elif self.cutoffAxis == 'x':

            # pick up where samples loop left off
            if os.path.exists(
              os.path.join(cycles_dir, 'sample_drop_idxs.pkl')):
                f = open(
                    os.path.join(cycles_dir, 'sample_drop_idxs.pkl'), 'rb'
                    )

                sample_drop_idxs = pickle.load(f)

                samples_to_threshold = (
                    len(data['Sample'].unique())
                    - len(sample_drop_idxs.keys())
                    )

                print(f'Samples to threshold: {samples_to_threshold}')
            else:
                # initialize a dictionary to append indices to drop
                stt = len(data['Sample'].unique())
                print(f'Samples to threshold: {stt}')
                sample_drop_idxs = {}

            # drop samples previously run
            ratios_melt = ratios_melt[
                ~ratios_melt['sample'].isin(sample_drop_idxs.keys())
                ]

            # indices_to_drop list
            indices_to_drop = []

            for name, group in ratios_melt.groupby(['sample']):
                for cycle_ratio in group['cycle'].unique():
                    if not (
                      cycle_ratio.split('/')[0]
                      == cycle_ratio.split('/')[1]):

                        # isolate ratio data
                        cycle_data = group[group['cycle'] == cycle_ratio]

                        cycle_num = cycle_ratio.split('/')[1]

                        num_bins = 300
                        histtype = 'stepfilled'
                        sns.set_style('whitegrid')

                        if not cycle_data.empty:
                            cycle_data.to_parquet(
                                os.path.join(cycles_dir, 'cycle_data.parquet')
                                )

                            fig, ax = plt.subplots()
                            plt.subplots_adjust(left=0.25, bottom=0.25)
                            counts, bins, patches = plt.hist(
                                cycle_data['log10(ratio)'], bins=num_bins,
                                density=False, color='grey', ec='none',
                                alpha=0.75, histtype=histtype,
                                range=None, label='before'
                                )

                            plt.title(
                                f'Sample={name}   log10({dna1}/{cycle_num})',
                                size=10
                                )
                            plt.ylabel('count', size=10)

                            axcolor = 'lightgoldenrodyellow'
                            axLowerCutoff = plt.axes(
                                [0.25, 0.15, 0.65, 0.03], facecolor=axcolor
                                )
                            axUpperCutoff = plt.axes(
                                [0.25, 0.1, 0.65, 0.03], facecolor=axcolor
                                )

                            rnge = [bins.min(), bins.max()]

                            sLower = Slider(
                                axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
                                valinit=0.00, valstep=(rnge[1]/100000)
                                )
                            sLower.label.set_color('b')

                            sUpper = Slider(
                                axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
                                valinit=0.00, valstep=(rnge[1]/100000)
                                )
                            sUpper.label.set_color('r')

                            def update(val):
                                [i.remove() for i in ax.get_lines()]
                                lowerCutoff = sLower.val
                                upperCutoff = sUpper.val
                                blueLine = ax.axvline(
                                    x=lowerCutoff, c='b', linewidth=2.5)
                                redLine = ax.axvline(
                                    x=upperCutoff, c='r', linewidth=2.5)
                                return lowerCutoff, upperCutoff

                            sLower.on_changed(update)
                            sUpper.on_changed(update)

                            resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
                            button = Button(
                                resetax, 'Reset', color=axcolor,
                                hovercolor='0.975'
                                )

                            def reset(event):
                                sLower.reset()
                                sUpper.reset()
                            button.on_clicked(reset)

                            def submit(text):
                                lowerCutoff, upperCutoff = update(val=None)

                                # read group
                                cycle_data = pd.read_parquet(
                                    os.path.join(
                                        cycles_dir, 'cycle_data.parquet')
                                    )

                                # get indices in log(ratio) series outside
                                # lower and upper cutoffs
                                idxs = list(
                                    cycle_data['index'][
                                        (cycle_data['log10(ratio)']
                                         < lowerCutoff) |
                                        (cycle_data['log10(ratio)']
                                         > upperCutoff)]
                                        )

                                if text == name:

                                    channel_number = markers['channel_number'][
                                                markers['marker_name']
                                                == f'{cycle_num}']

                                    dna_first = imread(
                                        f'{self.in_dir}/tif/{text}.*tif',
                                        key=0
                                        )

                                    dna_last = imread(
                                        f'{self.in_dir}/tif/{text}.*tif',
                                        key=channel_number - 1
                                        )

                                    # filter group data by selecting
                                    # indices NOT in idxs
                                    sample_data = data[data['Sample'] == name]
                                    drop_df = sample_data.index.isin(idxs)
                                    sample_centroids = sample_data[
                                        ['Y_centroid', 'X_centroid']][~drop_df]

                                    with napari.gui_qt():
                                        viewer = napari.view_image(
                                            dna_last, rgb=False,
                                            blending='additive',
                                            colormap='magenta',
                                            name=f'{cycle_num}'
                                            )
                                        viewer.add_image(
                                            dna_first, rgb=False,
                                            blending='additive',
                                            colormap='green', name=f'{dna1}'
                                            )
                                        viewer.add_points(
                                            sample_centroids,
                                            name='Stable Nuclei',
                                            properties=None,
                                            face_color='yellow',
                                            edge_color='k',
                                            edge_width=0.0, size=4.0
                                            )
                                else:
                                    print(
                                        'Must enter name of current' +
                                        f' sample: {name}.'
                                        )

                            axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
                            text_box = TextBox(
                                axbox, 'evaluation sample name', initial='',
                                color='0.95',
                                hovercolor='1.0',
                                label_pad=0.05
                                )
                            text_box.on_submit(submit)
                            plt.show(block=True)

                            # update sample_drop_idxs dictionary
                            lowerCutoff, upperCutoff = update(val=None)

                            cycle_data = pd.read_parquet(
                                os.path.join(cycles_dir, 'cycle_data.parquet')
                                )

                            # take all data if sliders not moved
                            if lowerCutoff == upperCutoff:
                                lowerCutoff = cycle_data['log10(ratio)'].min()
                                upperCutoff = cycle_data['log10(ratio)'].max()

                            idxs = list(
                                cycle_data['index'][
                                    (cycle_data['log10(ratio)'] < lowerCutoff)
                                    |
                                    (cycle_data['log10(ratio)'] > upperCutoff)]
                                    )
                            sample_drop_idxs[name] = idxs

                            # update drop indices pickle
                            os.chdir(cycles_dir)
                            f = open('sample_drop_idxs.pkl', 'wb')
                            pickle.dump(sample_drop_idxs, f)
                            f.close()

                            # filter dataframe by selecting indices NOT in the
                            for k, v in sample_drop_idxs.items():
                                for j in v:
                                    indices_to_drop.append(j)

                        else:
                            # filter dataframe by selecting indices NOT in the
                            indices_to_drop = []
                            for k, v in sample_drop_idxs.items():
                                for j in v:
                                    indices_to_drop.append(j)

            df = data.loc[~data.index.isin(indices_to_drop)]

        # grab dna and sample columns
        facet_input = df.loc[
            :, df.columns.str.contains(f'{dna_moniker}|Sample')].copy()

        facet_per_cycle_melt = (
            facet_input
            .sample(frac=1.0)
            .reset_index()
            .melt(id_vars=['Sample', 'index'], var_name='cycle',
                  )
            )

        facet_per_cycle_melt['Sample'] = pd.Categorical(
            facet_per_cycle_melt['Sample'], ordered=True,
            categories=natsorted(
                facet_per_cycle_melt['Sample'].unique()))

        facet_per_cycle_melt['cycle'] = pd.Categorical(
            facet_per_cycle_melt['cycle'], ordered=True,
            categories=natsorted(
                facet_per_cycle_melt['cycle'].unique()))

        facet_per_cycle_melt = (
            facet_per_cycle_melt.sort_values(['Sample', 'cycle', 'index'])
            )

        # plot dna intensity correlation per cycle
        fig, ax = plt.subplots(figsize=(5, 5))
        g = sns.FacetGrid(
            facet_per_cycle_melt, col='cycle', col_wrap=5,
            sharex=True, sharey=False
            )

        g.map(
            lambda y, color: plt.scatter(
                facet_per_cycle_melt['value'].loc[
                    facet_per_cycle_melt['cycle']
                    == f'{dna1}_{self.mask_object}'], y,
                s=0.05, alpha=0.1, linewidth=None,
                marker='o', c='r'), 'value')
        plt.savefig(
            os.path.join(
                cycles_dir, 'cycle_correlation(perCycle).png'), dpi=800
                )
        plt.close('all')

        # plot dna intensity correlation per cycle (color by sample)
        fig, ax = plt.subplots(figsize=(5, 5))

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(
                facet_per_cycle_melt['Sample'].unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False
            )

        sample_color_dict = dict(
            zip(
                natsorted(facet_per_cycle_melt['Sample'].unique()),
                cmap.colors)
                )

        g = sns.FacetGrid(
            facet_per_cycle_melt, col='cycle', hue='Sample',
            col_wrap=5, sharex=True, sharey=True
            )

        g.map(
            lambda sam, y, color, **kwargs: plt.scatter(
                facet_per_cycle_melt.loc[
                    (facet_per_cycle_melt['Sample'] ==
                     sam.unique()[0])
                    & (facet_per_cycle_melt['cycle'] ==
                       f'{dna1}_{self.mask_object}'),
                    'value'], y,
                c=np.reshape(sample_color_dict[sam.unique()[0]], (-1, 3)),
                s=0.05, linewidth=None, marker='o', **kwargs),
            'Sample', 'value'
            )

        plt.legend(markerscale=10, bbox_to_anchor=(1.1, 1.05))

        plt.savefig(
            os.path.join(
                cycles_dir, 'cycle_correlation(perSample).png'), dpi=800,
            bbox_inches='tight'
            )
        plt.close('all')
        print()

        # remove last sample groupby dataframe
        if os.path.exists(os.path.join(cycles_dir, 'cycle_data.parquet')):
            os.remove(os.path.join(cycles_dir, 'cycle_data.parquet'))

        return df

    @module
    def log10transform(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        # rescale abx intensities between 0 and 1
        # min_max_scaler = MinMaxScaler()
        # data[abx_channels] = min_max_scaler.fit_transform(data[abx_channels])

        # log10 transform
        data[abx_channels] += 0.00000000001
        data[abx_channels] = np.log10(data[abx_channels])
        print()

        return data

    @module
    def pruneOutliers(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        # save images
        pruning_dir = os.path.join(self.out_dir, 'pruning')
        if not os.path.exists(pruning_dir):
            os.mkdir(pruning_dir)

        # store upper and lower percentile cutoffs for antibody channels
        # pick up where left off
        if os.path.exists(os.path.join(pruning_dir, 'pruning_dict.pkl')):
            f = open(os.path.join(pruning_dir, 'pruning_dict.pkl'), 'rb')
            pruning_dict = pickle.load(f)
            total_markers = set(abx_channels)
            completed_markers = set(pruning_dict.keys())
            scrambled_markers_to_prune = total_markers.difference(
                completed_markers
                )
            # revert scrambled markers back to marker.csv order
            markers_dict = dict(zip(abx_channels, range(0, len(abx_channels))))
            sorted_marker_idxs = sorted(
                [markers_dict[i] for i in scrambled_markers_to_prune]
            )
            markers_to_prune = [
                list(markers_dict.keys())[list(markers_dict.values()).index(i)]
                for i in sorted_marker_idxs
                ]
            print(f'Immunomarker channels to prune: {len(markers_to_prune)}')
            print()
            if os.path.exists(os.path.join(pruning_dir, 'data_copy1.parquet')):
                os.chdir(pruning_dir)
                data_copy1 = pd.read_parquet('data_copy1.parquet')
        else:
            markers_to_prune = abx_channels
            print(f'Immunomarker channels to prune: {len(markers_to_prune)}')
            print()
            pruning_dict = {}
            data_copy1 = data.copy()

        # plot raw signal intensity distributions for inspection
        for ab in sorted(markers_to_prune):

            print(ab.split(f'_{self.mask_object}')[0])

            hist_facet = (
                data_copy1[['Sample', 'Condition', 'Area'] + [ab]]
                .sample(frac=1.0)
                .melt(
                    id_vars=['Sample', 'Condition', 'Area'],
                    var_name='channel', value_name='signal')
                )

            # sort and create labels column for plotting
            hist_facet.sort_values(by=['Condition', 'Sample'], inplace=True)
            hist_facet['for_plot'] = (
                hist_facet['Condition'] + ', ' +
                hist_facet['Sample']
                )

            sns.set_style('white')
            g = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=15,
                height=3.0, aspect=1.0, sharex=True, sharey=False,
                )

            if self.hexbins:
                g.map(
                    lambda x, y, color: plt.hexbin(
                        x, y, gridsize=self.hexbin_grid_size,
                        linewidths=0.02, color='dimgrey'),
                    'signal', 'Area')
            else:
                g.map(
                    lambda x, y, color: plt.scatter(
                        x, y, s=1.0, linewidths=0.0, color='k'),
                    'signal', 'Area')

            g.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=9.0, pad=2.0
                )

            g.fig.suptitle(
                ab.split(f'_{self.mask_object}')[0], y=1.1, size=20.0
                )

            for ax in g.axes.flatten():
                ax.tick_params(
                    axis='both', which='major',
                    labelsize=5.0, pad=-2
                    )
                ax.xaxis.label.set_size(10.0)
                ax.yaxis.label.set_size(10.0)

                if self.hexbins:
                    ax.spines['left'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                else:
                    ax.spines['left'].set_linewidth(0.1)
                    ax.spines['bottom'].set_linewidth(0.1)

            plt.subplots_adjust(
                left=0.01, bottom=0.01, right=0.99,
                top=0.90, hspace=0.4, wspace=0.4
                )

            plt.savefig(
                os.path.join(
                    pruning_dir,
                    f'{ab.split(f"_{self.mask_object}")[0]}_raw.pdf'),
                bbox_inches='tight')
            plt.close()

            subprocess.call(
                ['open', '-a', 'Preview', os.path.join(
                    pruning_dir,
                    f'{ab.split(f"_{self.mask_object}")[0]}_raw.pdf')]
                    )

            def submit(text):

                lowerPercentileCutoff = float(text.split(',')[0])
                upperPercentileCutoff = float(text.split(',')[1])

                pruning_dict[ab] = (
                    lowerPercentileCutoff,
                    upperPercentileCutoff
                    )

                os.chdir(pruning_dir)
                f = open(os.path.join(pruning_dir, 'pruning_dict.pkl'), 'wb')
                pickle.dump(pruning_dict, f)
                f.close()

                data_copy2 = (
                    data_copy1[['Sample', 'Condition', 'Area'] + [ab]].copy()
                    )

                channel_data = data_copy2[ab]

                indices_to_drop = []

                # add row index to list if column value < lower bound
                indices_to_drop.extend(
                    channel_data.index[
                        channel_data < np.percentile(
                            channel_data, lowerPercentileCutoff)]
                        )

                # add row index to list if column value > upper bound
                indices_to_drop.extend(
                    channel_data.index[
                        channel_data > np.percentile(
                            channel_data, upperPercentileCutoff)])

                data_copy2.drop(
                    labels=set(indices_to_drop), axis=0,
                    inplace=True, errors='raise'
                    )

                pruned_data = data_copy2[ab]

                # rescale channel signal intensities
                scaler = (
                    MinMaxScaler(feature_range=(0, 1), copy=True)
                    .fit(pruned_data.values.reshape(-1, 1))
                        )
                rescaled_data = scaler.transform(
                    pruned_data.values.reshape(-1, 1)
                    )

                rescaled_data = pd.DataFrame(
                    data=rescaled_data,
                    index=pruned_data.index,
                    ).rename(columns={0: ab})

                data_copy2.update(rescaled_data)

                # plot pruned and rescaled signal intensity histrograms
                hist_facet = (
                    data_copy2
                    .sample(frac=1.0)
                    .melt(
                        id_vars=['Sample', 'Condition', 'Area'],
                        var_name='channel', value_name='signal')
                    )

                # sort and create labels column for plotting
                hist_facet.sort_values(
                    by=['Condition', 'Sample'], inplace=True
                    )
                hist_facet['for_plot'] = (
                    hist_facet['Condition'] + ', ' +
                    hist_facet['Sample']
                    )

                g = sns.FacetGrid(
                    hist_facet, col='for_plot', col_wrap=15,
                    height=3.0, aspect=1.0, sharex=True, sharey=False
                    )

                if self.hexbins:
                    g.map(
                        lambda x, y, color: plt.hexbin(
                            x, y, gridsize=self.hexbin_grid_size,
                            linewidths=0.02, color='dimgrey'),
                        'signal', 'Area')
                else:
                    g.map(
                        lambda x, y, color: plt.scatter(
                            x, y, s=1.0, linewidths=0.0, color='k'),
                        'signal', 'Area')

                g.set_titles(
                    col_template="{col_name}", fontweight='bold',
                    size=9.0, pad=2.0
                    )

                g.fig.suptitle(
                    ab.split(f'_{self.mask_object}')[0], y=1.1, size=20.0
                    )

                for ax in g.axes.flatten():
                    ax.tick_params(
                        axis='both', which='major',
                        labelsize=5.0, pad=-2
                        )
                    ax.xaxis.label.set_size(10.0)
                    ax.yaxis.label.set_size(10.0)

                    if self.hexbins:
                        ax.spines['left'].set_visible(False)
                        ax.spines['bottom'].set_visible(False)
                    else:
                        ax.spines['left'].set_linewidth(0.1)
                        ax.spines['bottom'].set_linewidth(0.1)

                plt.subplots_adjust(
                    left=0.01, bottom=0.01, right=0.99,
                    top=0.90, hspace=0.4, wspace=0.4
                    )

                plt.savefig(
                    os.path.join(
                        pruning_dir,
                        f'{ab.split(f"_{self.mask_object}")[0]}' +
                        '_pruned_rescaled.pdf'),
                    bbox_inches='tight'
                        )
                plt.close()

                subprocess.call(
                    ['open', '-a', 'Preview', os.path.join(
                        pruning_dir,
                        f'{ab.split(f"_{self.mask_object}")[0]}' +
                        '_pruned_rescaled.pdf')]
                        )
                ##################

            plt.rcParams['figure.figsize'] = (6, 2)
            axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
            text_box = TextBox(
                axbox,
                'lowerCutoff, upperCutoff',
                initial='',
                color='0.95',
                hovercolor='1.0',
                label_pad=0.05
                )
            text_box.label.set_size(12)

            text_box.on_submit(submit)

            plt.show(block=True)

            if ab not in pruning_dict.keys():
                # default crop percentiles if not explicitly set
                pruning_dict[ab] = (0.0, 100.0)

                os.chdir(pruning_dir)
                f = open(os.path.join(pruning_dir, 'pruning_dict.pkl'), 'wb')
                pickle.dump(pruning_dict, f)
                f.close()

            channel_data = data_copy1[ab]

            indices_to_drop = []

            # add row index to list if column value < lower bound
            indices_to_drop.extend(
                channel_data.index[
                    channel_data < np.percentile(
                        channel_data,
                        pruning_dict[ab][0])]
                    )

            # add row index to list if column value > upper bound
            indices_to_drop.extend(
                channel_data.index[
                    channel_data > np.percentile(
                        channel_data,
                        pruning_dict[ab][1])])

            data_copy1.drop(
                labels=set(indices_to_drop), axis=0,
                inplace=True, errors='raise'
                )

            pruned_data = data_copy1[ab]

            # rescale channel signal intensities
            scaler = (
                MinMaxScaler(feature_range=(0, 1), copy=True)
                .fit(pruned_data.values.reshape(-1, 1))
                    )
            rescaled_data = scaler.transform(
                pruned_data.values.reshape(-1, 1)
                )

            rescaled_data = pd.DataFrame(
                data=rescaled_data,
                index=pruned_data.index,
                ).rename(columns={0: ab})

            data_copy1.update(rescaled_data)

            # plot pruned and rescaled signal intensity histrograms
            hist_facet = (
                data_copy1[['Sample', 'Condition', 'Area'] + [ab]]
                .sample(frac=1.0)
                .melt(
                    id_vars=['Sample', 'Condition', 'Area'],
                    var_name='channel', value_name='signal')
                )

            # sort and create labels column for plotting
            hist_facet.sort_values(
                by=['Condition', 'Sample'], inplace=True
                )
            hist_facet['for_plot'] = (
                hist_facet['Condition'] + ', ' +
                hist_facet['Sample']
                )

            g = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=15,
                height=3.0, aspect=1.0, sharex=True, sharey=False
                )

            if self.hexbins:
                g.map(
                    lambda x, y, color: plt.hexbin(
                        x, y, gridsize=self.hexbin_grid_size,
                        linewidths=0.02, color='dimgrey'),
                    'signal', 'Area')
            else:
                g.map(
                    lambda x, y, color: plt.scatter(
                        x, y, s=1.0, linewidths=0.0, color='k'),
                    'signal', 'Area')

            g.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=9.0, pad=2.0
                )

            g.fig.suptitle(
                ab.split(f'_{self.mask_object}')[0], y=1.1, size=20.0
                )

            for ax in g.axes.flatten():
                ax.tick_params(
                    axis='both', which='major',
                    labelsize=5.0, pad=-2
                    )
                ax.xaxis.label.set_size(10.0)
                ax.yaxis.label.set_size(10.0)

                if self.hexbins:
                    ax.spines['left'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                else:
                    ax.spines['left'].set_linewidth(0.1)
                    ax.spines['bottom'].set_linewidth(0.1)

            plt.subplots_adjust(
                left=0.01, bottom=0.01, right=0.99,
                top=0.90, hspace=0.4, wspace=0.4
                )

            plt.savefig(
                os.path.join(
                    pruning_dir,
                    f'{ab.split(f"_{self.mask_object}")[0]}' +
                    '_pruned_rescaled.pdf'),
                bbox_inches='tight'
                    )
            plt.close()

            os.chdir(pruning_dir)
            data_copy1.to_parquet('data_copy1.parquet')
        print()

        # apply percentile cutoffs to their respective channels
        # in the same order as originally curated
        for k, v in sorted(pruning_dict.items()):
            print(
                'Applying percentile cutoffs to the ' +
                f'{k.split(f"_{self.mask_object}")[0]} channel.'
                )
            channel_data = data[k]

            indices_to_drop = []

            # add row index to list if column value < lower bound
            indices_to_drop.extend(
                channel_data.index[
                    channel_data < np.percentile(
                        channel_data, v[0])]
                    )

            # add row index to list if column value > upper bound
            indices_to_drop.extend(
                channel_data.index[
                    channel_data > np.percentile(
                        channel_data, v[1])])

            data.drop(
                labels=set(indices_to_drop), axis=0,
                inplace=True, errors='raise'
                )

            pruned_data = data[k]

            # rescale channel signal intensities
            scaler = (
                MinMaxScaler(feature_range=(0, 1), copy=True)
                .fit(pruned_data.values.reshape(-1, 1))
                    )
            rescaled_data = scaler.transform(
                pruned_data.values.reshape(-1, 1)
                )

            rescaled_data = pd.DataFrame(
                data=rescaled_data,
                index=pruned_data.index,
                ).rename(columns={0: k})

            data.update(rescaled_data)

        print()

        # remove data_copy1 file (equivalent to data returned by this module)
        # note: data_copy1 is updated after each channel is
        # pruned/rescaled for dynamic restarts
        if os.path.exists(os.path.join(pruning_dir, 'data_copy1.parquet')):
            os.remove(os.path.join(pruning_dir, 'data_copy1.parquet'))

        return data

    @module
    def PCA(data, self, args):

        if len(data['Sample'].unique()) > 1:

            markers, dna1, dna_moniker, abx_channels = read_markers(
                markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
                mask_object=self.mask_object,
                markers_to_exclude=self.markers_to_exclude,
                )

            abx_channels = [
                i for i in abx_channels if i not in self.channelExclusionsPCA
                ]

            medians = (
                data
                .groupby(['Sample'])
                .median()[abx_channels]
                )
            medians = medians.reindex(natsorted(medians.index))

            # cancer_cores = [
            #     k.split('-')[1] for (k, v) in self.sample_statuses.items()
            #     if v.split('-')[1] == 'TRUE'
            #     ]
            # medians = medians[medians.index.isin(cancer_cores)]
            medians = medians[~medians.index.isin(self.samplesToRemovePCA)]

            # specify PCA parameters
            pca = PCA(self.dimensionPCA, random_state=1)

            idx = medians.index

            # normalize signal intensities across samples (axis=0)
            if self.normalize is True:
                medians = norm(
                    medians, norm='l2', axis=0, copy=True, return_norm=False
                    )
            else:
                medians = medians.values

            # apply PCA parameters to data
            projected = pca.fit_transform(medians)

            # generate dataframe for plot input
            scatter_input = pd.DataFrame(data=projected, index=idx)
            scatter_input.rename(columns={0: 'PC1', 1: 'PC2'}, inplace=True)

            # get sample metadata keys from configuration file
            metadata_keys = [
                i for i in self.sample_conditions.keys()
                if i.split('-', 1)[1] in scatter_input.index
                ]

            # set sample abbreviations as row index
            scatter_input.index = [
                self.sample_abbreviations[i]
                for i in metadata_keys
                ]
            scatter_input['condition'] = [
                self.sample_conditions[i]
                for i in metadata_keys
                ]

            ordered_samples = (
                natsorted(
                    set(scatter_input.index.unique())
                    .difference(set(self.samplesToSilhouette)))
                    )

            # build cmap
            cmap = categorical_cmap(
                numUniqueSamples=len(ordered_samples),
                numCatagories=10,
                cmap='tab10',
                continuous=False
                )

            sample_color_dict = dict(
                zip(
                    ordered_samples,
                    cmap.colors)
                    )

            for s in self.samplesToSilhouette:
                sample_color_dict[s] = [0.863, 0.863, 0.863]

            # loop over plot input data
            for e, sample in enumerate(scatter_input.iterrows()):
                if sample[0] in self.samplesToSilhouette:
                    sns.set_style('whitegrid')
                    data_point = pd.DataFrame(scatter_input.iloc[e]).T
                    g = sns.scatterplot(
                        data=data_point, x='PC1', y='PC2',
                        hue=data_point.index,
                        palette=sample_color_dict,
                        edgecolor='gainsboro', linewidth=0.2, zorder=2,
                        s=self.pointSize, alpha=1.0, legend=False
                        )
                else:
                    # plot scores plot for first 2 PCs
                    sns.set_style('whitegrid')
                    data_point = pd.DataFrame(scatter_input.iloc[e]).T
                    g = sns.scatterplot(
                        data=data_point, x='PC1', y='PC2',
                        hue=data_point.index,
                        palette=sample_color_dict,
                        edgecolor='k', linewidth=0.2, zorder=3,
                        s=self.pointSize, alpha=1.0, legend=False
                        )

            g.grid(color='gray', linewidth=0.05, linestyle='-', alpha=1.0)
            plt.setp(g.spines.values(), color='k', lw=0.5)

            scatter_input = scatter_input.reset_index().rename(
                columns={'index': 'abbreviation'}
                )

            # annotate data points
            if self.labelPoints is True:

                # generate squareform distance matrix
                sq = squareform(
                    pdist(scatter_input[['PC1', 'PC2']], metric='euclidean')
                    )
                # squareform matrix with numerical labels
                d = pd.DataFrame(
                    sq, index=scatter_input.index,
                    columns=scatter_input.index
                    )
                # upper triangle
                d1 = d.where(np.triu(np.ones(d.shape)).astype(np.bool))
                # apply distance cutoff
                d2 = d1[d1 < self.distanceCutoff]
                # flatten, set multi-index as columns (core1, core2)
                d3 = (
                    d2
                    .stack()
                    .reset_index()
                    .rename(
                        columns={
                            'level_0': 'core1_id', 'level_1': 'core2_id',
                            0: 'dist'
                            }
                        )
                    )
                d3['core1_tissue'] = [
                    scatter_input.loc[i, 'abbreviation']
                    for i in d3['core1_id']
                    ]
                d3['core2_tissue'] = [
                    scatter_input.loc[i, 'abbreviation']
                    for i in d3['core2_id']
                    ]
                # drop diagonal values
                d4 = d3[d3['dist'] != 0.0].dropna()
                d5 = pd.DataFrame(columns=d4.columns)
                for i in d4.iterrows():
                    if i[1]['core1_tissue'] == i[1]['core2_tissue']:
                        d5 = d5.append(i[1])

                # list of proximal cores
                proximal_labels = set(
                    list(d5['core1_id']) + list(d5['core2_id'])
                    )

                unique_labels = set()
                neighbors_set = set()
                for e, (label, x, y) in enumerate(zip(
                  scatter_input['abbreviation'],
                  scatter_input['PC1'], scatter_input['PC2'])):

                    if e in proximal_labels:
                        if e not in neighbors_set:
                            d6 = (
                                d5[
                                    (d5['core1_id'] == e) |
                                    (d5['core2_id'] == e)
                                    ]
                                [['core1_id', 'core2_id']]
                                )
                            neighbors = set(
                                list(d6['core1_id']) + list(d6['core2_id'])
                                )
                            neighbors_set = neighbors_set.union(neighbors)

                            neighbors_df = scatter_input.loc[neighbors]

                            pc1 = neighbors_df['PC1']
                            pc2 = neighbors_df['PC2']
                            centroid = (
                                sum(pc1) / len(neighbors_df),
                                sum(pc2) / len(neighbors_df)
                                )

                            if label not in self.samplesToSilhouette:

                                text = plt.annotate(
                                    label,
                                    xy=(centroid[0], centroid[1]),
                                    xytext=(0, 0), size=4.75,
                                    fontweight='bold',
                                    color=sample_color_dict[label],
                                    textcoords='offset points', ha='center',
                                    va='center',
                                    bbox=dict(
                                        boxstyle='round,pad=0.1',
                                        fc='yellow', alpha=0.0)
                                        )
                                unique_labels.add(label)

                                text.set_path_effects(
                                    [path_effects.Stroke(
                                        linewidth=0.75, foreground='k'),
                                     path_effects.Normal()]
                                    )
                    else:
                        if label not in self.samplesToSilhouette:

                            text = plt.annotate(
                                label, xy=(x, y),
                                xytext=(0, 0), size=4.75, fontweight='bold',
                                color=sample_color_dict[label],
                                textcoords='offset points', ha='center',
                                va='center', bbox=dict(
                                    boxstyle='round,pad=0.1',
                                    fc='yellow', alpha=0.0)
                                    )
                            unique_labels.add(label)

                            text.set_path_effects(
                                [path_effects.Stroke(
                                    linewidth=0.75, foreground='k'),
                                 path_effects.Normal()]
                                )

            # get n per tissue type according to full sample names
            n_per_tissue_type = (
                scatter_input
                .groupby(['condition'])
                .count()
                .reindex(scatter_input['condition'])['PC1'].values
                )
            scatter_input['n'] = n_per_tissue_type

            legend_data = (
                scatter_input
                .drop(
                    ['PC1', 'PC2'], axis=1)
                .drop_duplicates()
                )

            # natural sort by abbreviated sample names
            legend_data['abbreviation'] = pd.Categorical(
                legend_data['abbreviation'],
                ordered=True,
                categories=natsorted(legend_data['abbreviation'].unique())
                )
            legend_data.sort_values('abbreviation', inplace=True)

            legend_handles = []
            for i in legend_data.iterrows():
                cond = i[1]['condition']
                abbr = i[1]['abbreviation']
                n = i[1]['n']

                if abbr in self.samplesToSilhouette:
                    markerfacecolor = 'gainsboro'
                    markeredgecolor = 'gainsboro'
                else:
                    markerfacecolor = sample_color_dict[abbr]
                    markeredgecolor = 'k'

                legend_handles.append(
                    Line2D([0], [0], marker='o', color='none',
                           label=f'{abbr} ({cond}, n={n})',
                           markerfacecolor=markerfacecolor,
                           markeredgecolor=markeredgecolor,
                           markeredgewidth=0.2,
                           markersize=5.0)
                           )

            g.legend(
                handles=legend_handles,
                prop={'size': 5.0},
                bbox_to_anchor=[1.03, 1.0]
                )

            plt.xlabel(
                f'PC1 ({round((pca.explained_variance_ratio_[0] * 100), 2)}'
                '% of variance)', fontsize=10, labelpad=7.0)
            plt.ylabel(
                f'PC2 ({round((pca.explained_variance_ratio_[1] * 100), 2)}'
                '% of variance)', fontsize=10, labelpad=4.0)
            plt.tick_params(axis='both', which='major', labelsize=7.0)
            plt.savefig(
                os.path.join(self.out_dir, 'pcaScoresPlot.pdf'),
                bbox_inches='tight'
                )
            plt.close('all')
        print()

        return data

    @module
    def clustering1(data, self, args):

        clustering_dir = os.path.join(self.out_dir, 'clustering/clustering1')
        if not os.path.exists(clustering_dir):
            os.makedirs(clustering_dir)

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering1
            ]

        if os.path.exists(os.path.join(clustering_dir, 'embedding.npy')):

            # recapitulate df index at the point of embedding

            data = data[~data['Sample'].isin(self.samplesToRemoveClustering1)]

            if self.normalizeTissueCounts1:

                # calculate per tissue cell-count weighted random sample
                groups = data.groupby('Sample')
                sample_weights = pd.DataFrame({
                    'weights': 1 / (groups.size() * len(groups))
                })
                weights = pd.merge(
                    data[['Sample']], sample_weights,
                    left_on='Sample', right_index=True
                    )

                df = data.sample(
                    frac=self.fracForEmbedding1, replace=False,
                    weights=weights['weights'], random_state=5, axis=0
                    )
                print('Tissue counts normalized')

            else:

                df = data.sample(frac=self.fracForEmbedding1, random_state=5)

            df.reset_index(drop=True, inplace=True)
            print(df[abx_channels])

            embedding = np.load(os.path.join(clustering_dir, 'embedding.npy'))
            df['emb1'] = embedding[:, 0]
            df['emb2'] = embedding[:, 1]

        else:
            startTime = datetime.now()

            data = data[~data['Sample'].isin(self.samplesToRemoveClustering1)]

            if self.normalizeTissueCounts1:

                # calculate per tissue cell-count weighted random sample
                groups = data.groupby('Sample')
                sample_weights = pd.DataFrame({
                    'weights': 1 / (groups.size() * len(groups))
                })
                weights = pd.merge(
                    data[['Sample']], sample_weights,
                    left_on='Sample', right_index=True
                    )

                df = data.sample(
                    frac=self.fracForEmbedding1, replace=False,
                    weights=weights['weights'], random_state=5, axis=0
                    )
                print('Tissue counts normalized')

            else:

                df = data.sample(frac=self.fracForEmbedding1, random_state=5)

            df.reset_index(drop=True, inplace=True)
            print(df[abx_channels])

            if self.embeddingAlgorithm1 == 'TSNE':
                print('Computing TSNE embedding.')
                embedding = TSNE(
                    n_components=self.dimensionEmbedding1,
                    perplexity=self.perplexity1,
                    early_exaggeration=self.earlyExaggeration1,
                    learning_rate=self.learningRateTSNE1,
                    metric=self.metric1,
                    random_state=self.random_state1,
                    init='pca', n_jobs=-1).fit_transform(df[abx_channels])

            elif self.embeddingAlgorithm1 == 'UMAP':
                print('Computing UMAP embedding.')
                embedding = UMAP(
                    n_components=self.dimensionEmbedding1,
                    n_neighbors=self.nNeighbors1,
                    learning_rate=self.learningRateUMAP1,
                    output_metric=self.metric1,
                    min_dist=self.minDist1,
                    repulsion_strength=self.repulsionStrength1,
                    random_state=3,
                    n_epochs=1000,
                    init='spectral',
                    metric='euclidean',
                    metric_kwds=None,
                    output_metric_kwds=None,
                    n_jobs=-1,
                    low_memory=False,
                    spread=1.0,
                    local_connectivity=1.0,
                    set_op_mix_ratio=1.0,
                    negative_sample_rate=5,
                    transform_queue_size=4.0,
                    a=None,
                    b=None,
                    angular_rp_forest=False,
                    target_n_neighbors=-1,
                    target_metric='categorical',
                    target_metric_kwds=None,
                    target_weight=0.5,
                    transform_seed=42,
                    transform_mode='embedding',
                    force_approximation_algorithm=False,
                    verbose=False,
                    unique=False,
                    densmap=False,
                    dens_lambda=2.0,
                    dens_frac=0.6,
                    dens_var_shift=0.1,
                    disconnection_distance=None,
                    output_dens=False).fit_transform(df[abx_channels])

            print('Embedding completed in ' + str(datetime.now() - startTime))

            np.save(os.path.join(clustering_dir, 'embedding'), embedding)
            df['emb1'] = embedding[:, 0]
            df['emb2'] = embedding[:, 1]

        sns.set_style('white')

        def submit(text):

            markers, dna1, dna_moniker, abx_channels = read_markers(
                markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
                mask_object=self.mask_object,
                markers_to_exclude=self.markers_to_exclude,
                )
            abx_channels = [
                i for i in abx_channels
                if i not in self.channelExclusionsClustering1
                ]

            numerical_input = text.split('.')[0].strip()
            tup = tuple(map(int, numerical_input.split('-')))

            if len(tup) == 1:

                mylist = [tup[0]]

                for i in mylist:

                    min_cluster_size = i

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=min_cluster_size, min_samples=None,
                        metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,
                        memory=Memory(
                            location=None),
                        approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False).fit(
                            df[['emb1', 'emb2']]
                            )
                    df['cluster'] = clustering.labels_

                    print(
                        f'min_cluster_size={i}', np.unique(clustering.labels_)
                        )

                    plt.rcParams['figure.figsize'] = (15, 7)
                    fig, (ax1, ax2) = plt.subplots(1, 2)

                    plt.subplots_adjust(
                        wspace=0.7,
                        left=0.04
                        )

                    # PLOT TSNE
                    for color_by in ['cluster', 'Condition']:

                        highlight = 'none'

                        if color_by == 'cluster':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(df[color_by].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            # make black the first color to specify
                            # cluster outliers (i.e. cluster -1 cells)
                            # colors_list.insert(0, (0.0, 0.0, 0.0))
                            cmap = ListedColormap(
                                np.insert(
                                    arr=cmap.colors, obj=0,
                                    values=[0, 0, 0], axis=0)
                                    )

                            # trim qualitative cmap to number of unique samples
                            trim = (
                                len(cmap.colors) - len(df[color_by].unique())
                                )
                            cmap = ListedColormap(
                                cmap.colors[:-trim]
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df[color_by].unique()),
                                    list(range(len(df[color_by].unique()))))
                                    )

                            c = [sample_dict[i] for i in df[color_by]]

                            ax1.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=105000/len(df),
                                ec=[
                                    'k' if i == highlight else 'none' for
                                    i in df[color_by]
                                    ],
                                linewidth=0.1
                                )

                            ax1.axis('equal')
                            ax1.tick_params(labelsize=5)
                            ax1.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(df[color_by].unique())
                              ):

                                hi_markers = cluster_expression(
                                    df=df, markers=abx_channels,
                                    cluster=i, num_proteins=3,
                                    across_or_within='within',
                                    )

                                legend_elements.append(
                                    Line2D([0], [0], marker='o', color='none',
                                           label=f'Cluster: {i} {hi_markers}',
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor='none', lw=0.001,
                                           markersize=4)
                                           )

                            cluster_lgd = ax1.legend(
                                handles=legend_elements,
                                prop={'size': 7},
                                bbox_to_anchor=[1.02, 1.0]
                                )

                        elif color_by == 'Condition':

                            # label PCA plot points by condition and replicate
                            df['label'] = (
                                df[color_by] + '_' + df['Replicate'].map(str)
                                )

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(df['label'].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df['label'].unique()),
                                    list(range(len(df['label'].unique()))))
                                    )

                            c = [sample_dict[i] for i in df['label']]

                            ax2.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=105000/len(df),
                                ec=[
                                    'k' if i == highlight else 'none' for
                                    i in df['label']
                                    ],
                                linewidth=0.1
                                )

                            ax2.axis('equal')
                            ax2.tick_params(labelsize=5)
                            ax2.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(df['label'].unique())
                              ):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                uv = unmicst_version(self.sample_conditions)

                                sample_to_map = df['Sample'][
                                    df['label'] == i].unique()[0]
                                abbr = self.sample_abbreviations[
                                    f"{uv}-{sample_to_map }"
                                    ]

                                legend_elements.append(
                                    Line2D([0], [0], marker='o', color='none',
                                           label=f'Sample: {abbr} ({i})',
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001,
                                           markersize=4)
                                           )

                            sample_lgd = ax2.legend(
                                handles=legend_elements,
                                prop={'size': 7},
                                bbox_to_anchor=[1.27, 1.0]
                                )

                    plt.tight_layout()

                    if '.save' in text:

                        plt.savefig(
                            os.path.join(
                                clustering_dir,
                                f'{self.embeddingAlgorithm1}_'
                                f'{min_cluster_size}.png'),
                            bbox_extra_artists=(cluster_lgd, sample_lgd),
                            bbox_inches='tight', dpi=1000
                            )
                        plt.close('all')

                    plt.show(block=False)

            else:

                # df = df.sample(frac=0.01, random_state=22)

                mylist = list(range(tup[0], tup[1] + 1, 1))
                mylist.reverse()  # run higher sizes first for plot order

                for i in mylist:

                    min_cluster_size = i

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=min_cluster_size, min_samples=None,
                        metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,
                        memory=Memory(
                            location=None),
                        approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False).fit(
                            df[['emb1', 'emb2']]
                            )
                    df['cluster'] = clustering.labels_

                    print(
                        f'min_cluster_size={i}', np.unique(clustering.labels_)
                        )

        plt.rcParams['figure.figsize'] = (6, 2)
        axbox = plt.axes([0.6, 0.525, 0.35, 0.15])
        text_box = TextBox(
            axbox,
            'min_cluster_size (single # or range #-# .save)',
            initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.label.set_size(12)
        text_box.on_submit(submit)
        plt.show(block=True)
        print()

        df.drop(columns='label', inplace=True)

        return df

    @module
    def getClustermap1(data, self, args):

        clustering_dir = os.path.join(self.out_dir, 'clustering/clustering1')
        if not os.path.exists(clustering_dir):
            os.makedirs(clustering_dir)

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        clustermap_input = data[data['cluster'] != -1]

        cluster_heatmap_input = clustermap_input[
            abx_channels + ['cluster']].groupby('cluster').mean()

        sns.set(font_scale=0.8)
        g = sns.clustermap(
            cluster_heatmap_input, cmap='viridis', standard_scale=1,
            square=False, yticklabels=1, linewidth=0.1, cbar=True
            )

        plt.gcf().set_size_inches(8.0, 8.0)

        plt.savefig(
            os.path.join(
                clustering_dir, 'clustermap.pdf'), bbox_inches='tight')

        plt.show(block=True)
        print()

        return data

    @module
    def lassoClusters(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        subplot_kw = dict(
            xlim=(data['emb1'].min(), data['emb1'].max()),
            ylim=(data['emb2'].min(), data['emb2'].max()),
            autoscale_on=False)

        subplot_kw = dict()

        plt.rcParams['figure.figsize'] = (8, 8)

        fig, ax = plt.subplots(subplot_kw=subplot_kw)

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(data['cluster'].unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False
            )
        # add black as first element to represent HDBSCAN outliers
        cmap.colors = np.append([[0.0, 0.0, 0.0]], cmap.colors, axis=0)
        cmap.colors = np.delete(cmap.colors, -1, axis=0)

        pts = ax.scatter(
            data['emb1'],
            data['emb2'],
            c=data['cluster'],
            cmap=cmap,
            s=2.0,
            ec='none'
            )
        # pts = ax.scatter(
        #     lasso_X['emb1'], lasso_X['emb2'], c='lime', s=0.0005
        #     )
        selector = SelectFromCollection(ax, pts)

        def accept(event):
            if event.key == "enter":
                # print("Selected points:")
                # print(selector.xys[selector.ind])
                selector.disconnect()
                ax.set_title("")
                fig.canvas.draw()

        fig.canvas.mpl_connect("key_press_event", accept)
        ax.set_title(
            f'{self.embeddingAlgorithm1} embedding. ' +
            'Press enter to accept selected points.')
        ax.set_aspect('equal')
        plt.show(block=True)

        # Isolate indices from "selector object" that WERE NOT lassoed.
        # selector_idxs_to_drop = (
        #     set(np.array(list(range(selector.xys.data.shape[0]))))
        #     - set(selector.ind))

        # Use iloc to isolate the indices of the original dataframe
        # corresponding to those in the "selector object".
        # (This step relies on the fact that the "selector object"
        # maintains the original dataframe row order.)
        df_idxs_lassoed = data.iloc[selector.ind].index
        # df_idxs_to_drop = df.iloc[list(selector_idxs_to_drop)].index

        # show highest expression channels
        markers = data.copy()
        markers.loc[markers.index.isin(df_idxs_lassoed), 'cluster'] = 1000
        # markers.loc[~markers.index.isin(df_idxs_to_drop), 'cluster'] = 1000
        if len(markers[markers['cluster'] == 1000]) > 0:
            hi_markers = cluster_expression(
                df=markers, markers=abx_channels, cluster=1000,
                num_proteins=3, across_or_within='within',
                )
            print(hi_markers)
        print()

        return data

    @module
    def frequencyStats1(data, self, args):

        stats_input = data[['Sample', 'Replicate', 'cluster']][
            data['cluster'] >= 0]

        for i in range(
          len(list(self.sample_statuses.values())[0].split(', '))):

            comparison = set(
                [j.split(', ')[i] for j in self.sample_statuses.values()
                 if '-UNK' not in j.split(', ')[i]]
                )

            test = [i for i in comparison if i not in self.controlGroups][0]
            control = [i for i in comparison if i in self.controlGroups][0]

            frequency_dir = os.path.join(
                self.out_dir, 'clustering/clustering1/frequency_stats',
                f"{test}_v_{control}"
                )
            if not os.path.exists(frequency_dir):
                os.makedirs(frequency_dir)

            # create single-column dataFrame containing all sample names
            # to pad counts tables with zeros if a celltype is not in a tissue
            pad = pd.DataFrame(
                natsorted(stats_input['Sample'].unique())).rename(
                    columns={0: 'Sample'}
                    )

            uv = unmicst_version(self.sample_conditions)

            cluster_list = []
            ratio_list = []
            dif_list = []
            pval_list = []

            # intialize a dataframe to collect catplot data
            catplot_input = pd.DataFrame()

            # loop over clusters
            for w, group in stats_input.groupby('cluster'):

                print(
                    f'Calculating log2({test}/{control})'
                    f' of mean cell density for cluster {str(w)}.'
                    )

                group = (
                    group.groupby(['Sample', 'Replicate', 'cluster'])
                    .size()
                    .reset_index(drop=False)
                    .rename(columns={0: 'count'})
                    )

                group = (
                    group
                    .merge(pad, how='right', on='Sample')
                    .sort_values(by='count', ascending=False)
                    )

                # guard against NaNs induced by the absence
                # of a given cluster in one or
                # more of the tissue samples
                group['count'] = [
                    0 if np.isnan(i) else int(i) for i in group['count']]

                group['status'] = [
                    self.sample_statuses[f"{uv}-{j}"]
                    .split(', ')[i] for j in group['Sample']
                    ]

                group['Replicate'] = [
                    self.sample_replicates[f"{uv}-{i}"]
                    for i in group['Sample']
                    ]

                group['cluster'] = w

                group = group[~group['status'].str.contains('-UNK')]

                group.reset_index(drop=True, inplace=True)

                # get denominator cell count for each sample
                if self.denominatorCluster1 is None:
                    group['tissue_count'] = [
                        len(stats_input[stats_input['Sample'] == i]) for
                        i in group['Sample']]
                else:
                    group['tissue_count'] = [
                        len(stats_input[(stats_input['Sample'] == i) &
                            (stats_input['cluster'] ==
                             self.denominatorCluster1)])
                        for i in group['Sample']
                        ]

                # compute density of cells per sample
                group['density'] = group['count']/group['tissue_count']

                # append group data to catplot_input
                catplot_input = catplot_input.append(group)

                cnd1_values = group['density'][group['status'] == test]
                cnd2_values = group['density'][group['status'] == control]

                # Welch's t-test (equal_var=False)
                stat, pval = ttest_ind(
                    cnd1_values, cnd2_values,
                    axis=0, equal_var=False, nan_policy='propagate')

                stat = round(stat, 3)
                pval = round(pval, 3)

                cnd1_mean = np.mean(cnd1_values)
                cnd2_mean = np.mean(cnd2_values)

                ratio = np.log2((cnd1_mean + 0.000001)/(cnd2_mean + 0.000001))

                dif = cnd1_mean-cnd2_mean

                cluster_list.append(w)
                ratio_list.append(ratio)
                dif_list.append(dif)
                pval_list.append(pval)

            statistics = pd.DataFrame(
                list(zip(cluster_list, ratio_list, dif_list, pval_list)),
                columns=['cluster', 'ratio', 'dif', 'pval']).sort_values(
                    by='cluster')

            statistics.to_csv(
                os.path.join(
                    frequency_dir, 'stats_total.csv'), index=False)

            stats = importr('stats')

            p_adjust = stats.p_adjust(
                FloatVector(statistics['pval'].tolist()),
                method='BH')

            statistics['qval'] = p_adjust

            if self.FDRCorrection:
                stat = 'qval'

            else:
                stat = 'pval'

            significant = statistics[
                statistics[stat] <= 0.05].sort_values(by=stat)

            significant.to_csv(
                os.path.join(
                    frequency_dir, 'stats_sig.csv'), index=False)

            sns.set_style('whitegrid')
            fig, ax = plt.subplots()
            plt.scatter(abs(significant['dif']), significant['ratio'])

            for label, qval, x, y in zip(
              significant['cluster'], significant[stat],
              abs(significant['dif']), significant['ratio']):

                plt.annotate(
                    (label, f'{stat[0]}=' + str(qval)), size=3,
                    xy=(x, y), xytext=(10, 10),
                    textcoords='offset points', ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.1', fc='yellow',
                              alpha=0.0))

            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(5)

            plt.title(
                f'{test} vs. {control} ({stat[0]}<0.05)',
                fontsize=12
                )
            plt.xlabel(
                f'abs({test} - {control})', fontsize=10
                )
            plt.ylabel(
                f'log2({test} / {control})',
                fontsize=10
                )
            plt.savefig(os.path.join(frequency_dir, 'plot.pdf'))
            plt.close()

            catplot_input.reset_index(drop=True, inplace=True)

            catplot_input[stat] = [
                 'ns' if i not in significant['cluster'].unique() else
                 significant[stat][significant['cluster'] == i].values[0]
                 for i in catplot_input['cluster']]

            # build cmap
            cmap = categorical_cmap(
                numUniqueSamples=len(catplot_input['Sample'].unique()),
                numCatagories=10,
                cmap='tab10',
                continuous=False
                )

            sample_color_dict = dict(
                zip(
                    natsorted(catplot_input['Sample'].unique()),
                    cmap.colors)
                    )

            catplot_input.sort_values(
                by=['cluster', 'status', 'density'],
                ascending=[True, False, True], inplace=True
                )

            catplot_input['cluster'] = (
                catplot_input['cluster'].astype(str) + f'; {stat} = ' +
                catplot_input[stat].astype(str)
                )
            catplot_input['cluster'] = [
                i.split(f'; {stat} = ns')[0] for i in catplot_input['cluster']
                ]

            sns.set(font_scale=0.4)
            g = sns.catplot(
                x='status', y='density',
                hue=catplot_input['Sample'], col='cluster', col_wrap=6,
                data=catplot_input, kind='bar', palette=sample_color_dict,
                height=2, aspect=0.8, sharex=True, sharey=False,
                edgecolor='k', linewidth=0.1,
                legend=False
                )

            g.set(ylim=(0.0, None))

            sample_conds = [
                self.sample_conditions[f'{uv}-{i}']
                for i in natsorted(catplot_input['Sample'].unique())
                ]

            sample_abbrs = [
                self.sample_abbreviations[f'{uv}-{i}']
                for i in natsorted(catplot_input['Sample'].unique())
                ]

            cond_abbr = [
                f'{i}-{j}' for i, j in zip(sample_conds, sample_abbrs)
                ]

            handles_dict = dict(
                zip(natsorted(catplot_input['Sample'].unique()),
                    cond_abbr)
                    )

            legend_handles = []
            for k, v in handles_dict.items():
                legend_handles.append(
                    Line2D([0], [0], marker='o', color='none',
                           label=v, markerfacecolor=sample_color_dict[k],
                           markeredgecolor='k', markeredgewidth=0.2,
                           markersize=5.0)
                           )

            plt.legend(
                handles=legend_handles,
                prop={'size': 5.0},
                bbox_to_anchor=[1.03, 1.0]
                )

            plt.savefig(
                os.path.join(frequency_dir, 'catplot.pdf'),
                bbox_inches='tight'
                )
            plt.close('all')
            print()

        return data

    @module
    def setContrast(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        dna = imread(
            f'{self.in_dir}/tif/{self.view_sample}.*tif', key=0
            )

        if os.path.exists(
          os.path.join(self.out_dir, 'contrast_limits.yml')):

            contrast_limits = yaml.safe_load(
                open(f'{self.out_dir}/contrast_limits.yml')
                )

            print ('Applying existing channel contrast settings.')

            with napari.gui_qt():
                viewer = napari.view_image(
                    dna, rgb=False,
                    blending='additive', colormap='gray',
                    name=f'{dna1}:{self.view_sample}'
                    )

                viewer.layers[f'{dna1}:{self.view_sample}'].contrast_limits = (
                    contrast_limits[dna1][0], contrast_limits[dna1][1]
                    )

                for ch in abx_channels:
                    ch = ch.split(f'_{self.mask_object}')[0]
                    channel_number = markers['channel_number'][
                                markers['marker_name'] == ch]

                    # read antibody image
                    img = imread(
                        f'{self.in_dir}/tif/{self.view_sample}.*tif',
                        key=(channel_number-1)
                        )

                    viewer.add_image(
                        img, rgb=False, blending='additive',
                        colormap='green', visible=False,
                        name=ch
                        )

                    viewer.layers[ch].contrast_limits = (
                        contrast_limits[ch][0], contrast_limits[ch][1]
                        )

        else:

            print ('Channel contrast settings have not been defined.')

            with napari.gui_qt():
                viewer = napari.view_image(
                    dna, rgb=False,
                    blending='additive', colormap='gray',
                    name=f'{dna1}:{self.view_sample}'
                    )

                for ch in abx_channels:
                    ch = ch.split(f'_{self.mask_object}')[0]
                    channel_number = markers['channel_number'][
                                markers['marker_name'] == ch]

                    # read antibody image
                    img = imread(
                        f'{self.in_dir}/tif/{self.view_sample}.*tif',
                        key=(channel_number-1)
                        )

                    viewer.add_image(
                        img, rgb=False, blending='additive',
                        colormap='green', visible=False,
                        name=ch
                        )

        # create channel settings configuration file,
        # update with chosen constrast limits
        contrast_limits = {}
        for ch in [dna1] + abx_channels:
            if ch == dna1:
                contrast_limits[ch] = (
                    viewer.layers[f'{dna1}:{self.view_sample}'].contrast_limits
                    )
            else:
                ch = ch.split(f'_{self.mask_object}')[0]
                contrast_limits[ch] = viewer.layers[ch].contrast_limits

        with open(f'{self.out_dir}/contrast_limits.yml', 'w') as file:
            yaml.dump(contrast_limits, file)
        print()

        return data

    @module
    def curateThumbnails1(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering1
            ]

        thumbnails_dir = os.path.join(
            self.out_dir, 'clustering/clustering1/thumbnails'
            )
        if not os.path.exists(thumbnails_dir):
            os.mkdir(thumbnails_dir)

        # sort ome.tif files from largest to smallest in size
        os.chdir(f'{self.in_dir}/tif/')
        ome_tifs = os.listdir(os.getcwd())
        ome_tifs.sort(key=lambda f: os.stat(f).st_size, reverse=True)
        ome_tifs = [
            i for i in ome_tifs if i.split('.')[0] in data['Sample'].unique()
            ]

        # grab image contrast limit settings
        if os.path.exists(f'{self.out_dir}/contrast_limits.yml'):
            contrast_limits = yaml.safe_load(
                open(f'{self.out_dir}/contrast_limits.yml')
                )

        df = data[data['cluster'] != -1]

        # store the numbers of those clusters that have already been run to
        # pick up where left off
        if os.path.exists(os.path.join(thumbnails_dir, 'thmbnls_to_run.pkl')):
            f = open(os.path.join(thumbnails_dir, 'thmbnls_to_run.pkl'), 'rb')
            thmbnls_to_run = pickle.load(f)
            total_clusters = set(df['cluster'].unique())
            completed_clusters = set(thmbnls_to_run)
            clusters_to_run = natsorted(
                total_clusters.difference(completed_clusters)
                )
            print(f'Clusters to run: {len(clusters_to_run)}')
            print()
        else:
            clusters_to_run = natsorted(df['cluster'].unique())
            print(f'Clusters to run: {len(clusters_to_run)}')
            print()
            thmbnls_to_run = []

        for cluster in clusters_to_run:

            print(f'Cluster: {cluster}')

            markers_to_show = cluster_expression(
                df=data, markers=abx_channels,
                cluster=cluster, num_proteins=3,
                across_or_within='within',
                )

            markers_to_show = [
                '_'.join(i.split('_')[0:-1]) if '_' in i else
                i for i in [dna1] + markers_to_show
                ]

            color_dict = {}
            for i, j, k in zip(
              markers_to_show,

              [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0),
               (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
              ['gray', 'green', 'red', 'blue']
              ):
                color_dict[i] = j

            long_table = pd.DataFrame()

            for sample_name in ome_tifs:
                print(f'Sample: {sample_name}')

                # import cycle1 dna channel, convert to float and rgb

                dna = imread(
                    os.path.join(f'{self.in_dir}/tif/', sample_name),
                    key=0
                    )

                dna = img_as_float(dna)
                dna = gray2rgb(dna)

                # loop over the channels to create a dict of images
                for e, marker in enumerate(markers_to_show):
                    if e != 0:
                        channel_number = markers['channel_number'][
                                    markers['marker_name'] == marker]

                        # read antibody image
                        img = imread(
                            os.path.join(
                                f'{self.in_dir}/tif/', sample_name),
                            key=(channel_number-1)
                            )

                        img = img_as_float(img)

                        img -= (contrast_limits[marker][0]/65535)
                        img /= (
                            (contrast_limits[marker][1]/65535)
                            - (contrast_limits[marker][0]/65535)
                            )
                        img = np.clip(img, 0, 1)

                        img = gray2rgb(img)

                        img = (img * color_dict[marker])

                        # loop over ab channels, add to cycle1 dna
                        dna += img

                        print(f'overlayed {marker} image')
                        clearRAM(print_usage=False)
                        del img
                        clearRAM(print_usage=False)

                # crop out thumbnail images
                sample_cluster_subset = data[
                    (data['Sample'] == sample_name.split('.')[0])
                    & (data['cluster'] == cluster)
                    ]

                sample_cluster_subset.reset_index(
                    drop=True, inplace=True
                    )

                if self.numThumbnails > len(sample_cluster_subset):
                    dif = (
                        self.numThumbnails
                        - len(sample_cluster_subset)
                        )

                    extra_rows = pd.DataFrame(
                        data=0,
                        index=list(range(dif)),
                        columns=sample_cluster_subset.columns
                        )
                    sample_cluster_subset = (
                        sample_cluster_subset.append(extra_rows)
                        )
                    sample_cluster_subset.reset_index(
                        drop=True, inplace=True
                        )
                else:
                    sample_cluster_subset = (
                        sample_cluster_subset.sample(
                            n=self.numThumbnails, random_state=3)
                        )
                # add centroid mask to image overlay
                centroids = sample_cluster_subset[
                    ['X_centroid', 'Y_centroid']
                    ]

                clearRAM(print_usage=False)
                del sample_cluster_subset
                clearRAM(print_usage=False)

                centroid_img = np.zeros(
                    (dna.shape[0],
                     dna.shape[1]))

                centroid_dist = 1  # in pixels
                for example, centroid in enumerate(
                  centroids.iterrows()):

                    ystart_centroid = int(
                        centroid[1]['Y_centroid'] - centroid_dist
                        )
                    ystop_centroid = int(
                        centroid[1]['Y_centroid'] + centroid_dist
                        )

                    xstart_centroid = int(
                        centroid[1]['X_centroid'] - centroid_dist
                        )
                    xstop_centroid = int(
                        centroid[1]['X_centroid'] + centroid_dist
                        )

                    centroid_img[
                        ystart_centroid:ystop_centroid,
                        xstart_centroid:xstop_centroid
                        ] = 1

                # convert to rgb and colorize
                centroid_img = gray2rgb(centroid_img)
                centroid_img = (centroid_img * (1.0, 1.0, 1.0))

                # add to overlay
                dna += centroid_img

                print('overlayed centroids')
                clearRAM(print_usage=False)
                del centroid_img
                clearRAM(print_usage=False)

                # crop thumbnails
                for example, centroid in enumerate(
                  centroids.iterrows()):

                    if (
                        (centroid[1]['X_centroid'] == 0.0) &
                        (centroid[1]['Y_centroid'] == 0.0)
                    ):

                        blank_img = np.ones(
                            (self.squareWindowDimension,
                             self.squareWindowDimension))

                        long_table = long_table.append(
                            {'sample': sample_name.split('.')[0],
                             'example': example,
                             'image': blank_img},
                            ignore_index=True
                            )

                    else:

                        # specify window x, y ranges
                        ystart_window = int(
                            centroid[1]['Y_centroid']
                            - self.squareWindowDimension
                            )
                        ystop_window = int(
                            centroid[1]['Y_centroid']
                            + self.squareWindowDimension
                            )

                        xstart_window = int(
                            centroid[1]['X_centroid']
                            - self.squareWindowDimension
                            )
                        xstop_window = int(
                            centroid[1]['X_centroid']
                            + self.squareWindowDimension
                            )

                        # for centroids falling within
                        # self.squareWindowDimension pixels of the edge of the
                        # dna image, ensure that the thumbnail image is not
                        # cropped using negative slicing values, as it will
                        # return an empty array and lead to a runtime error
                        # during plotting.
                        window_list = [
                            ystart_window, ystop_window,
                            xstart_window, xstop_window
                            ]
                        (ystart_window,
                         ystop_window,
                         xstart_window,
                         xstop_window) = [
                            0 if i < 0 else i for i in window_list
                            ]

                        # crop overlay image to window size
                        thumbnail = dna[
                            ystart_window:ystop_window,
                            xstart_window:xstop_window
                            ]

                        long_table = long_table.append(
                            {'sample': sample_name.split('.')[0],
                             'example': example,
                             'image': thumbnail.copy()},
                            ignore_index=True
                            )
                        # thumbnail.copy() so overlay image can be gc'd

                # print('before dell dna')
                # clearRAM(print_usage=False)
                # del dna
                # print('after dell dna')
                # clearRAM(print_usage=False)

            long_table['example'] = [
                int(i) for i in long_table['example']
                ]

            # plot facet grid
            fig, ax = plt.subplots()

            g = sns.FacetGrid(
                long_table, row='sample', col='example',
                sharex=False, sharey=False,
                gridspec_kws={'hspace': 0.1, 'wspace': 0.05})

            g.map(
                lambda x, **kwargs: (
                    plt.imshow(x.values[0]), plt.grid(False)), 'image')

            for ax in g.axes.flatten():
                ax.get_xaxis().set_ticks([])
                ax.set_xlabel('')
                ax.get_yaxis().set_ticks([])
                ax.set_ylabel('')

            g.set_titles(
                col_template="Ex. {col_name}",
                row_template="Smpl. {row_name}",
                fontweight='bold', size=8)

            custom_lines = []
            for k, v in color_dict.items():
                custom_lines.append(
                    Line2D([0], [0], color=v, lw=6))

            ax.legend(
                custom_lines,
                list(color_dict.keys()), prop={'size': 12},
                bbox_to_anchor=(
                    1.05, len(long_table['sample'].unique()) + 0.3),
                loc='upper left'
                )

            plt.savefig(
                os.path.join(
                    thumbnails_dir,
                    'cluster' + str(cluster) + '_thumbnails.pdf')
                    )
            plt.close('all')

            thmbnls_to_run.append(cluster)

            os.chdir(thumbnails_dir)
            f = open(os.path.join(thumbnails_dir, 'thmbnls_to_run.pkl'), 'wb')
            pickle.dump(thmbnls_to_run, f)
            f.close()

            # clearRAM(print_usage=False)
            # del g
            # clearRAM(print_usage=False)
            print()

        return data

    @module
    def dropClusters(data, self, args):

        print(f'Dropping clusters {self.clustersToDrop} from the analysis.')
        data = data[~data['cluster'].isin(self.clustersToDrop)]

        return data

    @module
    def clustering2(data, self, args):

        clustering_dir = os.path.join(self.out_dir, 'clustering/clustering2')
        if not os.path.exists(clustering_dir):
            os.makedirs(clustering_dir)

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering2
            ]

        if os.path.exists(os.path.join(clustering_dir, 'embedding.npy')):

            # recapitulate df index at the point of embedding

            data = data[~data['Sample'].isin(self.samplesToRemoveClustering2)]

            if self.normalizeTissueCounts2:

                # calculate per tissue cell-count weighted random sample
                groups = data.groupby('Sample')
                sample_weights = pd.DataFrame({
                    'weights': 1 / (groups.size() * len(groups))
                })
                weights = pd.merge(
                    data[['Sample']], sample_weights,
                    left_on='Sample', right_index=True
                    )

                df = data.sample(
                    frac=self.fracForEmbedding2, replace=False,
                    weights=weights['weights'], random_state=5, axis=0
                    )
                print('Tissue counts normalized')

            else:

                df = data.sample(frac=self.fracForEmbedding2, random_state=5)

            df.reset_index(drop=True, inplace=True)
            print(df[abx_channels])

            embedding = np.load(os.path.join(clustering_dir, 'embedding.npy'))
            df['emb1'] = embedding[:, 0]
            df['emb2'] = embedding[:, 1]

        else:
            startTime = datetime.now()

            data = data[~data['Sample'].isin(self.samplesToRemoveClustering2)]

            if self.normalizeTissueCounts2:

                # calculate per tissue cell-count weighted random sample
                groups = data.groupby('Sample')
                sample_weights = pd.DataFrame({
                    'weights': 1 / (groups.size() * len(groups))
                })
                weights = pd.merge(
                    data[['Sample']], sample_weights,
                    left_on='Sample', right_index=True
                    )

                df = data.sample(
                    frac=self.fracForEmbedding2, replace=False,
                    weights=weights['weights'], random_state=5, axis=0
                    )
                print('Tissue counts normalized')

            else:

                df = data.sample(frac=self.fracForEmbedding2, random_state=5)

            df.reset_index(drop=True, inplace=True)
            print(df[abx_channels])

            if self.embeddingAlgorithm2 == 'TSNE':
                print('Computing TSNE embedding.')
                embedding = TSNE(
                    n_components=self.dimensionEmbedding2,
                    perplexity=self.perplexity2,
                    early_exaggeration=self.earlyExaggeration2,
                    learning_rate=self.learningRateTSNE2,
                    metric=self.metric2,
                    random_state=self.random_state2,
                    init='pca', n_jobs=-1).fit_transform(df[abx_channels])

            elif self.embeddingAlgorithm2 == 'UMAP':
                print('Computing UMAP embedding.')
                embedding = UMAP(
                    n_components=self.dimensionEmbedding2,
                    n_neighbors=self.nNeighbors2,
                    learning_rate=self.learningRateUMAP2,
                    output_metric=self.metric2,
                    min_dist=self.minDist2,
                    repulsion_strength=self.repulsionStrength2,
                    random_state=3,
                    n_epochs=1000,
                    init='spectral',
                    metric='euclidean',
                    metric_kwds=None,
                    output_metric_kwds=None,
                    n_jobs=-1,
                    low_memory=False,
                    spread=1.0,
                    local_connectivity=1.0,
                    set_op_mix_ratio=1.0,
                    negative_sample_rate=5,
                    transform_queue_size=4.0,
                    a=None,
                    b=None,
                    angular_rp_forest=False,
                    target_n_neighbors=-1,
                    target_metric='categorical',
                    target_metric_kwds=None,
                    target_weight=0.5,
                    transform_seed=42,
                    transform_mode='embedding',
                    force_approximation_algorithm=False,
                    verbose=False,
                    unique=False,
                    densmap=False,
                    dens_lambda=2.0,
                    dens_frac=0.6,
                    dens_var_shift=0.1,
                    disconnection_distance=None,
                    output_dens=False).fit_transform(df[abx_channels])

            print('Embedding completed in ' + str(datetime.now() - startTime))

            np.save(os.path.join(clustering_dir, 'embedding'), embedding)
            df['emb1'] = embedding[:, 0]
            df['emb2'] = embedding[:, 1]

        sns.set_style('white')

        def submit(text):

            markers, dna1, dna_moniker, abx_channels = read_markers(
                markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
                mask_object=self.mask_object,
                markers_to_exclude=self.markers_to_exclude,
                )
            abx_channels = [
                i for i in abx_channels
                if i not in self.channelExclusionsClustering2
                ]

            numerical_input = text.split('.')[0].strip()
            tup = tuple(map(int, numerical_input.split('-')))

            if len(tup) == 1:

                mylist = [tup[0]]

                for i in mylist:

                    min_cluster_size = i

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=min_cluster_size, min_samples=None,
                        metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,
                        memory=Memory(
                            location=None),
                        approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False).fit(
                            df[['emb1', 'emb2']]
                            )
                    df['cluster'] = clustering.labels_

                    print(
                        f'min_cluster_size={i}', np.unique(clustering.labels_)
                        )

                    plt.rcParams['figure.figsize'] = (15, 7)
                    fig, (ax1, ax2) = plt.subplots(1, 2)

                    plt.subplots_adjust(
                        wspace=0.7,
                        left=0.04
                        )

                    # PLOT TSNE
                    for color_by in ['cluster', 'Condition']:

                        highlight = 'none'

                        if color_by == 'cluster':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(df[color_by].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            # make black the first color to specify
                            # cluster outliers (i.e. cluster -1 cells)
                            # colors_list.insert(0, (0.0, 0.0, 0.0))
                            cmap = ListedColormap(
                                np.insert(
                                    arr=cmap.colors, obj=0,
                                    values=[0, 0, 0], axis=0)
                                    )

                            # trim qualitative cmap to number of unique samples
                            trim = (
                                len(cmap.colors) - len(df[color_by].unique())
                                )
                            cmap = ListedColormap(
                                cmap.colors[:-trim]
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df[color_by].unique()),
                                    list(range(len(df[color_by].unique()))))
                                    )

                            c = [sample_dict[i] for i in df[color_by]]

                            ax1.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=105000/len(df),
                                ec=[
                                    'k' if i == highlight else 'none' for
                                    i in df[color_by]
                                    ],
                                linewidth=0.1
                                )

                            ax1.axis('equal')
                            ax1.tick_params(labelsize=5)
                            ax1.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(df[color_by].unique())
                              ):

                                hi_markers = cluster_expression(
                                    df=df, markers=abx_channels,
                                    cluster=i, num_proteins=3,
                                    across_or_within='within',
                                    )

                                legend_elements.append(
                                    Line2D([0], [0], marker='o', color='none',
                                           label=f'Cluster: {i} {hi_markers}',
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor='none', lw=0.001,
                                           markersize=4)
                                           )

                            cluster_lgd = ax1.legend(
                                handles=legend_elements,
                                prop={'size': 7},
                                bbox_to_anchor=[1.02, 1.0]
                                )

                        elif color_by == 'Condition':

                            # label PCA plot points by condition and replicate
                            df['label'] = (
                                df[color_by] + '_' + df['Replicate'].map(str)
                                )

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(df['label'].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df['label'].unique()),
                                    list(range(len(df['label'].unique()))))
                                    )

                            c = [sample_dict[i] for i in df['label']]

                            ax2.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=105000/len(df),
                                ec=[
                                    'k' if i == highlight else 'none' for
                                    i in df['label']
                                    ],
                                linewidth=0.1
                                )

                            ax2.axis('equal')
                            ax2.tick_params(labelsize=5)
                            ax2.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(df['label'].unique())
                              ):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                uv = unmicst_version(self.sample_conditions)

                                sample_to_map = df['Sample'][
                                    df['label'] == i].unique()[0]
                                abbr = self.sample_abbreviations[
                                    f"{uv}-{sample_to_map }"
                                    ]

                                legend_elements.append(
                                    Line2D([0], [0], marker='o', color='none',
                                           label=f'Sample: {abbr} ({i})',
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001,
                                           markersize=4)
                                           )

                            sample_lgd = ax2.legend(
                                handles=legend_elements,
                                prop={'size': 7},
                                bbox_to_anchor=[1.27, 1.0]
                                )

                    plt.tight_layout()

                    if '.save' in text:

                        plt.savefig(
                            os.path.join(
                                clustering_dir,
                                f'{self.embeddingAlgorithm2}_'
                                f'{min_cluster_size}.png'),
                            bbox_extra_artists=(cluster_lgd, sample_lgd),
                            bbox_inches='tight', dpi=1000
                            )
                        plt.close('all')

                    plt.show(block=False)

            else:

                # df = df.sample(frac=0.01, random_state=22)

                mylist = list(range(tup[0], tup[1] + 1, 1))
                mylist.reverse()  # run higher sizes first for plot order

                for i in mylist:

                    min_cluster_size = i

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=min_cluster_size, min_samples=None,
                        metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,
                        memory=Memory(
                            location=None),
                        approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False).fit(
                            df[['emb1', 'emb2']]
                            )
                    df['cluster'] = clustering.labels_

                    print(
                        f'min_cluster_size={i}', np.unique(clustering.labels_)
                        )

        plt.rcParams['figure.figsize'] = (6, 2)
        axbox = plt.axes([0.6, 0.525, 0.35, 0.15])
        text_box = TextBox(
            axbox,
            'min_cluster_size (single # or range #-# .save)',
            initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.label.set_size(12)
        text_box.on_submit(submit)
        plt.show(block=True)
        print()

        df.drop(columns='label', inplace=True)

        return df

    @module
    def getClustermap2(data, self, args):

        clustering_dir = os.path.join(self.out_dir, 'clustering/clustering2')
        if not os.path.exists(clustering_dir):
            os.makedirs(clustering_dir)

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )

        clustermap_input = data[data['cluster'] != -1]

        cluster_heatmap_input = clustermap_input[
            abx_channels + ['cluster']].groupby('cluster').mean()

        sns.set(font_scale=0.8)
        g = sns.clustermap(
            cluster_heatmap_input, cmap='viridis', standard_scale=1,
            square=False, yticklabels=1, linewidth=0.1, cbar=True
            )

        plt.gcf().set_size_inches(8.0, 8.0)

        plt.savefig(
            os.path.join(
                clustering_dir, 'clustermap.pdf'), bbox_inches='tight')

        plt.show(block=True)
        print()

        return data

    @module
    def frequencyStats2(data, self, args):

        stats_input = data[['Sample', 'Replicate', 'cluster']][
            data['cluster'] >= 0]

        for i in range(
          len(list(self.sample_statuses.values())[0].split(', '))):

            comparison = set(
                [j.split(', ')[i] for j in self.sample_statuses.values()
                 if '-UNK' not in j.split(', ')[i]]
                )

            test = [i for i in comparison if i not in self.controlGroups][0]
            control = [i for i in comparison if i in self.controlGroups][0]

            frequency_dir = os.path.join(
                self.out_dir, 'clustering/clustering2/frequency_stats',
                f"{test}_v_{control}"
                )
            if not os.path.exists(frequency_dir):
                os.makedirs(frequency_dir)

            # create single-column dataFrame containing all sample names
            # to pad counts tables with zeros if a celltype is not in a tissue
            pad = pd.DataFrame(
                natsorted(stats_input['Sample'].unique())).rename(
                    columns={0: 'Sample'}
                    )

            uv = unmicst_version(self.sample_conditions)

            cluster_list = []
            ratio_list = []
            dif_list = []
            pval_list = []

            # intialize a dataframe to collect catplot data
            catplot_input = pd.DataFrame()

            # loop over clusters
            for w, group in stats_input.groupby('cluster'):

                print(
                    f'Calculating log2({test}/{control})'
                    f' of mean cell density for cluster {str(w)}.'
                    )

                group = (
                    group.groupby(['Sample', 'Replicate', 'cluster'])
                    .size()
                    .reset_index(drop=False)
                    .rename(columns={0: 'count'})
                    )

                group = (
                    group
                    .merge(pad, how='right', on='Sample')
                    .sort_values(by='count', ascending=False)
                    )

                # guard against NaNs induced by the absence
                # of a given cluster in one or
                # more of the tissue samples
                group['count'] = [
                    0 if np.isnan(i) else int(i) for i in group['count']]

                group['status'] = [
                    self.sample_statuses[f"{uv}-{j}"]
                    .split(', ')[i] for j in group['Sample']
                    ]

                group['Replicate'] = [
                    self.sample_replicates[f"{uv}-{i}"]
                    for i in group['Sample']
                    ]

                group['cluster'] = w

                group = group[~group['status'].str.contains('-UNK')]

                group.reset_index(drop=True, inplace=True)

                # get denominator cell count for each sample
                if self.denominatorCluster2 is None:
                    group['tissue_count'] = [
                        len(stats_input[stats_input['Sample'] == i]) for
                        i in group['Sample']]
                else:
                    group['tissue_count'] = [
                        len(stats_input[(stats_input['Sample'] == i) &
                            (stats_input['cluster'] ==
                             self.denominatorCluster2)])
                        for i in group['Sample']
                        ]

                # compute density of cells per sample
                group['density'] = group['count']/group['tissue_count']

                # append group data to catplot_input
                catplot_input = catplot_input.append(group)

                cnd1_values = group['density'][group['status'] == test]
                cnd2_values = group['density'][group['status'] == control]

                # Welch's t-test (equal_var=False)
                stat, pval = ttest_ind(
                    cnd1_values, cnd2_values,
                    axis=0, equal_var=False, nan_policy='propagate')

                stat = round(stat, 3)
                pval = round(pval, 3)

                cnd1_mean = np.mean(cnd1_values)
                cnd2_mean = np.mean(cnd2_values)

                ratio = np.log2((cnd1_mean + 0.000001)/(cnd2_mean + 0.000001))

                dif = cnd1_mean-cnd2_mean

                cluster_list.append(w)
                ratio_list.append(ratio)
                dif_list.append(dif)
                pval_list.append(pval)

            statistics = pd.DataFrame(
                list(zip(cluster_list, ratio_list, dif_list, pval_list)),
                columns=['cluster', 'ratio', 'dif', 'pval']).sort_values(
                    by='cluster')

            statistics.to_csv(
                os.path.join(
                    frequency_dir, 'stats_total.csv'), index=False)

            stats = importr('stats')

            p_adjust = stats.p_adjust(
                FloatVector(statistics['pval'].tolist()),
                method='BH')

            statistics['qval'] = p_adjust

            if self.FDRCorrection:
                stat = 'qval'

            else:
                stat = 'pval'

            significant = statistics[
                statistics[stat] <= 0.05].sort_values(by=stat)

            significant.to_csv(
                os.path.join(
                    frequency_dir, 'stats_sig.csv'), index=False)

            sns.set_style('whitegrid')
            fig, ax = plt.subplots()
            plt.scatter(abs(significant['dif']), significant['ratio'])

            for label, qval, x, y in zip(
              significant['cluster'], significant[stat],
              abs(significant['dif']), significant['ratio']):

                plt.annotate(
                    (label, f'{stat[0]}=' + str(qval)), size=3,
                    xy=(x, y), xytext=(10, 10),
                    textcoords='offset points', ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.1', fc='yellow',
                              alpha=0.0))

            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(5)

            plt.title(
                f'{test} vs. {control} ({stat[0]}<0.05)',
                fontsize=12
                )
            plt.xlabel(
                f'abs({test} - {control})', fontsize=10
                )
            plt.ylabel(
                f'log2({test} / {control})',
                fontsize=10
                )
            plt.savefig(os.path.join(frequency_dir, 'plot.pdf'))
            plt.close()

            catplot_input.reset_index(drop=True, inplace=True)

            catplot_input[stat] = [
                 'ns' if i not in significant['cluster'].unique() else
                 significant[stat][significant['cluster'] == i].values[0]
                 for i in catplot_input['cluster']]

            # build cmap
            cmap = categorical_cmap(
                numUniqueSamples=len(catplot_input['Sample'].unique()),
                numCatagories=10,
                cmap='tab10',
                continuous=False
                )

            sample_color_dict = dict(
                zip(
                    natsorted(catplot_input['Sample'].unique()),
                    cmap.colors)
                    )

            catplot_input.sort_values(
                by=['cluster', 'status', 'density'],
                ascending=[True, False, True], inplace=True
                )

            catplot_input['cluster'] = (
                catplot_input['cluster'].astype(str) + f'; {stat} = ' +
                catplot_input[stat].astype(str)
                )
            catplot_input['cluster'] = [
                i.split(f'; {stat} = ns')[0] for i in catplot_input['cluster']
                ]

            sns.set(font_scale=0.4)
            g = sns.catplot(
                x='status', y='density',
                hue=catplot_input['Sample'], col='cluster', col_wrap=6,
                data=catplot_input, kind='bar', palette=sample_color_dict,
                height=2, aspect=0.8, sharex=True, sharey=False,
                edgecolor='k', linewidth=0.1,
                legend=False
                )

            g.set(ylim=(0.0, None))

            sample_conds = [
                self.sample_conditions[f'{uv}-{i}']
                for i in natsorted(catplot_input['Sample'].unique())
                ]

            sample_abbrs = [
                self.sample_abbreviations[f'{uv}-{i}']
                for i in natsorted(catplot_input['Sample'].unique())
                ]

            cond_abbr = [
                f'{i}-{j}' for i, j in zip(sample_conds, sample_abbrs)
                ]

            handles_dict = dict(
                zip(natsorted(catplot_input['Sample'].unique()),
                    cond_abbr)
                    )

            legend_handles = []
            for k, v in handles_dict.items():
                legend_handles.append(
                    Line2D([0], [0], marker='o', color='none',
                           label=v, markerfacecolor=sample_color_dict[k],
                           markeredgecolor='k', markeredgewidth=0.2,
                           markersize=5.0)
                           )

            plt.legend(
                handles=legend_handles,
                prop={'size': 5.0},
                bbox_to_anchor=[1.03, 1.0]
                )

            plt.savefig(
                os.path.join(frequency_dir, 'catplot.pdf'),
                bbox_inches='tight'
                )
            plt.close('all')
            print()

        return data

    @module
    def curateThumbnails2(data, self, args):

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
            mask_object=self.mask_object,
            markers_to_exclude=self.markers_to_exclude,
            )
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering2
            ]

        thumbnails_dir = os.path.join(
            self.out_dir, 'clustering/clustering2/thumbnails'
            )
        if not os.path.exists(thumbnails_dir):
            os.mkdir(thumbnails_dir)

        # sort ome.tif files from largest to smallest in size
        os.chdir(f'{self.in_dir}/tif/')
        ome_tifs = os.listdir(os.getcwd())
        ome_tifs.sort(key=lambda f: os.stat(f).st_size, reverse=True)
        ome_tifs = [
            i for i in ome_tifs if i.split('.')[0] in data['Sample'].unique()
            ]

        # grab image contrast limit settings
        if os.path.exists(f'{self.out_dir}/contrast_limits.yml'):
            contrast_limits = yaml.safe_load(
                open(f'{self.out_dir}/contrast_limits.yml')
                )

        df = data[data['cluster'] != -1]

        # store the numbers of those clusters that have already been run to
        # pick up where left off
        if os.path.exists(os.path.join(thumbnails_dir, 'thmbnls_to_run.pkl')):
            f = open(os.path.join(thumbnails_dir, 'thmbnls_to_run.pkl'), 'rb')
            thmbnls_to_run = pickle.load(f)
            total_clusters = set(df['cluster'].unique())
            completed_clusters = set(thmbnls_to_run)
            clusters_to_run = natsorted(
                total_clusters.difference(completed_clusters)
                )
            print(f'Clusters to run: {len(clusters_to_run)}')
            print()
        else:
            clusters_to_run = natsorted(df['cluster'].unique())
            print(f'Clusters to run: {len(clusters_to_run)}')
            print()
            thmbnls_to_run = []

        for cluster in clusters_to_run:

            print(f'Cluster: {cluster}')

            markers_to_show = cluster_expression(
                df=data, markers=abx_channels,
                cluster=cluster, num_proteins=3,
                across_or_within='within',
                )

            markers_to_show = [
                '_'.join(i.split('_')[0:-1]) if '_' in i else
                i for i in [dna1] + markers_to_show
                ]

            color_dict = {}
            for i, j, k in zip(
              markers_to_show,

              [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0),
               (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
              ['gray', 'green', 'red', 'blue']
              ):
                color_dict[i] = j

            long_table = pd.DataFrame()

            for sample_name in ome_tifs:
                print(f'Sample: {sample_name}')

                # import cycle1 dna channel, convert to float and rgb

                dna = imread(
                    os.path.join(f'{self.in_dir}/tif/', sample_name),
                    key=0
                    )

                dna = img_as_float(dna)
                dna = gray2rgb(dna)

                # loop over the channels to create a dict of images
                for e, marker in enumerate(markers_to_show):
                    if e != 0:
                        channel_number = markers['channel_number'][
                                    markers['marker_name'] == marker]

                        # read antibody image
                        img = imread(
                            os.path.join(
                                f'{self.in_dir}/tif/', sample_name),
                            key=(channel_number-1)
                            )

                        img = img_as_float(img)

                        img -= (contrast_limits[marker][0]/65535)
                        img /= (
                            (contrast_limits[marker][1]/65535)
                            - (contrast_limits[marker][0]/65535)
                            )
                        img = np.clip(img, 0, 1)

                        img = gray2rgb(img)

                        img = (img * color_dict[marker])

                        # loop over ab channels, add to cycle1 dna
                        dna += img

                        print(f'overlayed {marker} image')
                        clearRAM(print_usage=False)
                        del img
                        clearRAM(print_usage=False)

                # crop out thumbnail images
                sample_cluster_subset = data[
                    (data['Sample'] == sample_name.split('.')[0])
                    & (data['cluster'] == cluster)
                    ]

                sample_cluster_subset.reset_index(
                    drop=True, inplace=True
                    )

                if self.numThumbnails > len(sample_cluster_subset):
                    dif = (
                        self.numThumbnails
                        - len(sample_cluster_subset)
                        )

                    extra_rows = pd.DataFrame(
                        data=0,
                        index=list(range(dif)),
                        columns=sample_cluster_subset.columns
                        )
                    sample_cluster_subset = (
                        sample_cluster_subset.append(extra_rows)
                        )
                    sample_cluster_subset.reset_index(
                        drop=True, inplace=True
                        )
                else:
                    sample_cluster_subset = (
                        sample_cluster_subset.sample(
                            n=self.numThumbnails, random_state=3)
                        )
                # add centroid mask to image overlay
                centroids = sample_cluster_subset[
                    ['X_centroid', 'Y_centroid']
                    ]

                clearRAM(print_usage=False)
                del sample_cluster_subset
                clearRAM(print_usage=False)

                centroid_img = np.zeros(
                    (dna.shape[0],
                     dna.shape[1]))

                centroid_dist = 1  # in pixels
                for example, centroid in enumerate(
                  centroids.iterrows()):

                    ystart_centroid = int(
                        centroid[1]['Y_centroid'] - centroid_dist
                        )
                    ystop_centroid = int(
                        centroid[1]['Y_centroid'] + centroid_dist
                        )

                    xstart_centroid = int(
                        centroid[1]['X_centroid'] - centroid_dist
                        )
                    xstop_centroid = int(
                        centroid[1]['X_centroid'] + centroid_dist
                        )

                    centroid_img[
                        ystart_centroid:ystop_centroid,
                        xstart_centroid:xstop_centroid
                        ] = 1

                # convert to rgb and colorize
                centroid_img = gray2rgb(centroid_img)
                centroid_img = (centroid_img * (1.0, 1.0, 1.0))

                # add to overlay
                dna += centroid_img

                print('overlayed centroids')
                clearRAM(print_usage=False)
                del centroid_img
                clearRAM(print_usage=False)

                # crop thumbnails
                for example, centroid in enumerate(
                  centroids.iterrows()):

                    if (
                        (centroid[1]['X_centroid'] == 0.0) &
                        (centroid[1]['Y_centroid'] == 0.0)
                    ):

                        blank_img = np.ones(
                            (self.squareWindowDimension,
                             self.squareWindowDimension))

                        long_table = long_table.append(
                            {'sample': sample_name.split('.')[0],
                             'example': example,
                             'image': blank_img},
                            ignore_index=True
                            )

                    else:

                        # specify window x, y ranges
                        ystart_window = int(
                            centroid[1]['Y_centroid']
                            - self.squareWindowDimension
                            )
                        ystop_window = int(
                            centroid[1]['Y_centroid']
                            + self.squareWindowDimension
                            )

                        xstart_window = int(
                            centroid[1]['X_centroid']
                            - self.squareWindowDimension
                            )
                        xstop_window = int(
                            centroid[1]['X_centroid']
                            + self.squareWindowDimension
                            )

                        # for centroids falling within
                        # self.squareWindowDimension pixels of the edge of the
                        # dna image, ensure that the thumbnail image is not
                        # cropped using negative slicing values, as it will
                        # return an empty array and lead to a runtime error
                        # during plotting.
                        window_list = [
                            ystart_window, ystop_window,
                            xstart_window, xstop_window
                            ]
                        (ystart_window,
                         ystop_window,
                         xstart_window,
                         xstop_window) = [
                            0 if i < 0 else i for i in window_list
                            ]

                        # crop overlay image to window size
                        thumbnail = dna[
                            ystart_window:ystop_window,
                            xstart_window:xstop_window
                            ]

                        long_table = long_table.append(
                            {'sample': sample_name.split('.')[0],
                             'example': example,
                             'image': thumbnail.copy()},
                            ignore_index=True
                            )
                        # thumbnail.copy() so overlay image can be gc'd

                # print('before dell dna')
                # clearRAM(print_usage=False)
                # del dna
                # print('after dell dna')
                # clearRAM(print_usage=False)

            long_table['example'] = [
                int(i) for i in long_table['example']
                ]

            # plot facet grid
            fig, ax = plt.subplots()

            g = sns.FacetGrid(
                long_table, row='sample', col='example',
                sharex=False, sharey=False,
                gridspec_kws={'hspace': 0.1, 'wspace': 0.05})

            g.map(
                lambda x, **kwargs: (
                    plt.imshow(x.values[0]), plt.grid(False)), 'image')

            for ax in g.axes.flatten():
                ax.get_xaxis().set_ticks([])
                ax.set_xlabel('')
                ax.get_yaxis().set_ticks([])
                ax.set_ylabel('')

            g.set_titles(
                col_template="Ex. {col_name}",
                row_template="Smpl. {row_name}",
                fontweight='bold', size=8)

            custom_lines = []
            for k, v in color_dict.items():
                custom_lines.append(
                    Line2D([0], [0], color=v, lw=6))

            ax.legend(
                custom_lines,
                list(color_dict.keys()), prop={'size': 12},
                bbox_to_anchor=(
                    1.05, len(long_table['sample'].unique()) + 0.3),
                loc='upper left'
                )

            plt.savefig(
                os.path.join(
                    thumbnails_dir,
                    'cluster' + str(cluster) + '_thumbnails.pdf')
                    )
            plt.close('all')

            thmbnls_to_run.append(cluster)

            os.chdir(thumbnails_dir)
            f = open(os.path.join(thumbnails_dir, 'thmbnls_to_run.pkl'), 'wb')
            pickle.dump(thmbnls_to_run, f)
            f.close()

            # clearRAM(print_usage=False)
            # del g
            # clearRAM(print_usage=False)
            print()

        return data

    # @module
    # def makeZarrs(self, args):
    #
    #     if self.start_module == 'makeZarrs':
    #         df = read_dataframe(
    #             modules_list=self.modules,
    #             start_module=self.start_module,
    #             outDir=self.out_dir,
    #             )
    #     else:
    #         df = self.data
    #
    #     # if self.start_module == 'makeZarrs':
    #     #     df = pd.read_csv(
    #     #         os.path.join(
    #     #             self.out_dir,
    #     #             f'dataframe_archive/getSingleCellData.csv'),
    #     #         index_col=0
    #     #             )
    #     # else:
    #     #     df = self.data
    #
    #     markers, dna1, dna_moniker, abx_channels = read_markers(
    #         markers_filepath=os.path.join(self.in_dir, 'markers.csv'),
    #         mask_object=self.mask_object,
    #         markers_to_exclude=self.markers_to_exclude,
    #         )
    #
    #     # make directory to store zarr arrays
    #     zarrs_dir = os.path.join(self.out_dir, 'zarrs')
    #     if not os.path.exists(zarrs_dir):
    #         os.mkdir(zarrs_dir)
    #
    #     # initialize a dictionary to index and access zarr arrays
    #     zs = {}
    #
    #     # select data compression algorithm for zarr arrays
    #     compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE)
    #
    #     # loop over tissue samples
    #     for sample_name in df['Sample'].unique():
    #
    #         # read DNA ashlar_output image
    #         img = imread(
    #             f'{self.in_dir}/tif/{sample_name}.*tif',
    #             key=0  # key 0 is first dna cycle in markers.csv file
    #             )
    #
    #         # initialize zarr array for DNA image, append to a dictionary
    #         zs[f'{sample_name}_{dna1}'] = zarr.open(
    #             f'{zarrs_dir}/{sample_name}_{dna1}.zarr', mode='w',
    #             shape=(
    #                 img.shape[0], img.shape[1]),
    #             chunks=(1000, 1000), dtype='uint16', compressor=compressor
    #             )
    #
    #         # update zarr array with DNA image data
    #         zs[f'{sample_name}_{dna1}'][:] = img
    #
    #         # update zarr array with all antibody images
    #         # for k, v in cycle_map.items():
    #         abx_channels = [
    #             i for i in markers['marker_name'] if
    #             dna_moniker not in i
    #             ]
    #
    #         for ab in abx_channels:
    #
    #             print(
    #                 f'Generating Zarr array for sample {sample_name} {ab}.'
    #                 )
    #
    #             channel_number = markers['channel_number'][
    #                 markers['marker_name'] == ab].iloc[0]
    #
    #             # read antibody image
    #             img = imread(
    #                 f'{self.in_dir}/tif/{sample_name}.*tif',
    #                 key=(channel_number-1)
    #                 )
    #
    #             # initialize zarr array for image, append to a dictionary
    #             zs[f'{sample_name}_{ab}'] = zarr.open(
    #                 f'{zarrs_dir}/{sample_name}_{ab}.zarr', mode='w',
    #                 shape=(
    #                     img.shape[0], img.shape[1]),
    #                 chunks=(1000, 1000), dtype='uint16',
    #                 compressor=compressor
    #                 )
    #
    #             # update zarr array with image data
    #             zs[f'{sample_name}_{ab}'][:] = img
    #     print()
    #
    #     # apply image rendering settings if available and desired
    #     if os.path.exists(f'{self.out_dir}/contrast_limits.yml'):
    #
    #         # load settings
    #         contrast_limits = yaml.safe_load(
    #             open(f'{self.out_dir}/contrast_limits.yml'))
    #
    #         for k, v in zs.items():
    #             if k in contrast_limits.keys():
    #                 print(f'Applying channel signal intensity cutoffs for {k}')
    #
    #                 bottom_omero_cutoff = contrast_limits[k][0]
    #                 top_omero_cutoff = contrast_limits[k][1]
    #
    #                 temp = img_as_float(zs[k])
    #                 temp -= (bottom_omero_cutoff/65535)
    #                 temp /= (
    #                     (top_omero_cutoff/65535)-(bottom_omero_cutoff/65535)
    #                     )
    #                 temp = np.clip(temp, 0, 1)
    #
    #                 zs[k][:] = img_as_uint(temp)
    #
    #     save_dataframe(
    #         df=df, outDir=self.out_dir, moduleName='makeZarrs'
    #         )
    #     print()
    #     # self.data = df
    #     return data
    #
    # @module
    # def clusterBoxplots(self, args):
    #
    #     if self.start_module == 'clusterBoxplots':
    #         df = read_dataframe(
    #             modules_list=self.modules,
    #             start_module=self.start_module,
    #             outDir=self.out_dir,
    #             )
    #     else:
    #         df = self.data
    #
    #     cmap = categorical_cmap(
    #         numCatagories=10, numSubcatagories=2,
    #         cmap='tab10', continuous=False
    #         )
    #
    #     boxplot_input = df[df['cluster'] >= 0]
    #
    #     # get tidy input data
    #     boxplot_input = (
    #         boxplot_input[[
    #             'cluster', 'sample', 'condition'] + self.markers]
    #         .melt(
    #             id_vars=['cluster', 'sample', 'condition'],
    #             var_name='protein', value_name='log10(intensity)'))
    #
    #     boxplot_input = boxplot_input.sample(frac=1.0)
    #
    #     boxplot_input['label'] = (
    #         boxplot_input['protein'] + '_'
    #         + boxplot_input['cluster'].map(str) + '_'
    #         + boxplot_input['condition']
    #         )
    #
    #     boxplot_input.sort_values(by='label', inplace=True)
    #
    #     boxplot_input.rename(
    #         columns={'log10(intensity)': '$log_{10}(intensity)$'},
    #         inplace=True)
    #
    #     sns.set(font_scale=0.27)
    #     sns.set_style('whitegrid')
    #
    #     g = sns.FacetGrid(
    #         boxplot_input, row='cluster', col='protein',
    #         sharex=False, sharey=False, height=1, aspect=1.5)
    #
    #     hue_dict = dict(
    #         zip(boxplot_input['cluster'].unique(), cmap.colors)
    #         )
    #
    #     g.map(
    #         lambda x, y, z, color:
    #             sns.boxplot(
    #                 data=boxplot_input, x=x, y=y,
    #                 hue=z, palette=hue_dict,
    #                 linewidth=0.95, width=0.75,
    #                 fliersize=0.95),
    #             'condition', '$log_{10}(intensity)$', 'cluster')
    #
    #     def statistics(label):
    #
    #         conditions = sorted(list(set(self.samples.values())))
    #         cond_dict = {
    #             conditions[0]: conditions[1], conditions[1]: conditions[0]
    #             }
    #         label = label.values[0]
    #         cond_name = label.split('_')[-1]
    #         opposite_label = label.replace(cond_name, cond_dict[cond_name])
    #
    #         means1 = boxplot_input[
    #                 boxplot_input['label'] == label].groupby(
    #                     'sample').mean()['$log_{10}(intensity)$']
    #
    #         means2 = boxplot_input[
    #                 boxplot_input['label'] == opposite_label].groupby(
    #                     'sample').mean()['$log_{10}(intensity)$']
    #
    #         # perform Welch's unequal variances t-test
    #         stat, pval = ttest_ind(
    #             means1, means2,
    #             axis=0, equal_var=False, nan_policy='propagate')
    #
    #         if self.bonferroniCorrection:
    #
    #             # perform Bonferroni correction
    #             p_adj = pval * len(boxplot_input['label'].unique())/2
    #
    #             if p_adj <= 0.05:
    #                 ax = plt.gca()
    #                 ax.text(
    #                     0.5, 0.85,
    #                     r'$p_{adj} = $' + '%.1E' % Decimal(str(p_adj)),
    #                     fontweight='normal', fontsize=11.0,
    #                     color='k', ha='center', va='center',
    #                     transform=ax.transAxes
    #                     )
    #         else:
    #             # DO NOT perform Bonferroni correction
    #             if pval <= 0.05:
    #                 ax = plt.gca()
    #                 ax.text(
    #                     0.5, 0.85, r'$p = $' + '%.1E' % Decimal(str(pval)),
    #                     fontweight='normal', fontsize=11.0,
    #                     color='k', ha='center', va='center',
    #                     transform=ax.transAxes
    #                     )
    #
    #     g.map(
    #         lambda label, color:
    #             statistics(label=label), 'label')
    #
    #     plt.savefig(
    #         os.path.join(self.out_dir, 'cluster_boxplots.pdf'), bbox='tight')
    #     plt.close('all')
    #     print()
    #     # self.data = df
    #     return data
    #
    # @module
    # def spatialAnalysis(self, args):
    #
    #     if self.start_module == 'spatialAnalysis':
    #         df = read_dataframe(
    #             modules_list=self.modules,
    #             start_module=self.start_module,
    #             outDir=self.out_dir,
    #             )
    #     else:
    #         df = self.data
    #
    #     zs = loadZarrs(
    #         df=df, outDir=self.out_dir,
    #         markers_filepath=os.path.join(self.in_dir, 'markers.csv')
    #         )
    #
    #     spatial_dir = os.path.join(self.out_dir, 'spatial_analysis')
    #     if not os.path.exists(spatial_dir):
    #         os.makedirs(spatial_dir)
    #
    #     stats = pd.DataFrame(
    #         columns=[
    #             'protein', 'sample', 'celltype', 'ratio_r',
    #             'centroids', 'ratio_f', 'points']
    #         )
    #
    #     stats_row_idx = 0
    #     for protein, binary_cutoff in self.spatialDict1.items():
    #         for sample in df['sample'].unique():
    #
    #             if sample in self.cropDict.keys():
    #                 section = self.cropDict[sample][0]
    #                 cut_point = self.cropDict[sample][1]
    #
    #                 if section == 'bottom':
    #                     dna = zs[f'{sample}_dna'][cut_point:]
    #                     img = zs[f'{sample}_{protein}'][cut_point:]
    #
    #                 elif section == 'top':
    #                     dna = zs[f'{sample}_dna'][:cut_point]
    #                     img = zs[f'{sample}_{protein}'][:cut_point]
    #
    #             else:
    #                 dna = zs[f'{sample}_dna'][:]
    #                 img = zs[f'{sample}_{protein}'][:]
    #
    #             for celltype, cluster in self.spatialDict2.items():
    #
    #                 print(celltype, protein, sample)
    #
    #                 total_centroids = df[['x', 'y']][
    #                     (df['sample'] == sample) &
    #                     (df['cluster'] == cluster)].astype(int)
    #
    #                 if sample in self.cropDict.keys():
    #                     section = self.cropDict[sample][0]
    #                     cut_point = self.cropDict[sample][1]
    #
    #                     if section == 'bottom':
    #                         total_centroids = total_centroids[
    #                             total_centroids['y'] >= cut_point]
    #
    #                         total_centroids['y'] = (
    #                             total_centroids['y']-cut_point
    #                             )
    #
    #                     elif section == 'top':
    #                         total_centroids = total_centroids[
    #                             total_centroids['y'] < cut_point]
    #
    #                 if len(total_centroids) > 1:
    #                     y_min = total_centroids['y'].min()
    #                     y_max = total_centroids['y'].max()
    #                     y_range = y_max - y_min
    #
    #                     print(y_min, y_max)
    #
    #                     x_min = total_centroids['x'].min()
    #                     x_max = total_centroids['x'].max()
    #                     x_range = x_max - x_min
    #
    #                     dna_blurred = gaussian(dna, sigma=12)
    #                     dna_mask = np.where(dna_blurred > 0.05, 1, 0)
    #
    #                     area_mask = dna_mask[y_min:y_max, :]
    #                     area_mask = area_mask[:, x_min:x_max]
    #
    #                     inside_tumor = np.argwhere(area_mask == 1)
    #                     outside_tumor = np.argwhere(area_mask == 0)
    #
    #                     frac_area_out = outside_tumor.shape[0]/(
    #                         outside_tumor.shape[0] + inside_tumor.shape[0])
    #
    #                     img_blurred = gaussian(img, sigma=12)
    #                     img_mask_real = np.where(
    #                         img_blurred > binary_cutoff, 1, 0
    #                         )
    #                     img_mask_fake = img_mask_real.copy()
    #
    #                     (img_mask_real[total_centroids['y'],
    #                      total_centroids['x']]
    #                      ) += 10
    #
    #                     outside_centroids = np.argwhere(img_mask_real == 10)
    #                     outside_centroids = pd.DataFrame(
    #                         outside_centroids, columns=['y', 'x'])
    #
    #                     inside_centroids = np.argwhere(img_mask_real == 11)
    #                     inside_centroids = pd.DataFrame(
    #                         inside_centroids, columns=['y', 'x'])
    #
    #                     radii = []
    #                     num_points = []
    #                     for radius in range(
    #                         self.radiusRange[0], self.radiusRange[1]
    #                       ):
    #                         total_poisson_points = pd.DataFrame(
    #                             poisson_disc_samples(
    #                                 width=x_range, height=y_range, r=radius),
    #                             columns=['x', 'y'])
    #                         print(radius, len(total_poisson_points))
    #                         radii.append(radius)
    #                         num_points.append(len(total_poisson_points))
    #
    #                     def closest(lst, K):
    #                         return lst[
    #                             min(
    #                                 range(
    #                                     len(lst)),
    #                                 key=lambda i: abs(lst[i]-K))]
    #
    #                     optimal_points = closest(
    #                         num_points, len(total_centroids) +
    #                         (int(len(total_centroids) * frac_area_out)))
    #                     optimal_radius = radii[
    #                         num_points.index(optimal_points)-2
    #                         ]
    #
    #                     total_poisson_points = pd.DataFrame(
    #                         poisson_disc_samples(
    #                             width=x_range, height=y_range,
    #                             r=optimal_radius),
    #                         columns=['x', 'y']).astype(int)
    #
    #                     # ensure simulation contains at least 2 points
    #                     while len(total_poisson_points) < 2:
    #                         optimal_radius -= 1
    #                         total_poisson_points = pd.DataFrame(
    #                             poisson_disc_samples(
    #                                 width=x_range, height=y_range,
    #                                 r=optimal_radius),
    #                             columns=['x', 'y']).astype(int)
    #
    #                     total_poisson_points['x'] = (
    #                         total_poisson_points['x'] + x_min
    #                         )
    #                     total_poisson_points['y'] = (
    #                         total_poisson_points['y'] + y_min
    #                         )
    #
    #                     (dna_mask[total_poisson_points['y'],
    #                               total_poisson_points['x']]
    #                      ) += 10
    #
    #                     total_poisson_points = np.argwhere(dna_mask == 11)
    #                     total_poisson_points = pd.DataFrame(
    #                         total_poisson_points, columns=['y', 'x'])
    #
    #                     (img_mask_fake[total_poisson_points['y'],
    #                      total_poisson_points['x']]
    #                      ) += 10
    #
    #                     outside_poisson_points = np.argwhere(
    #                         img_mask_fake == 10
    #                         )
    #                     outside_poisson_points = pd.DataFrame(
    #                         outside_poisson_points, columns=['y', 'x'])
    #
    #                     inside_poisson_points = np.argwhere(
    #                         img_mask_fake == 11
    #                         )
    #                     inside_poisson_points = pd.DataFrame(
    #                         inside_poisson_points, columns=['y', 'x'])
    #
    #                     rgb_img = img_as_float(img)
    #                     rgb_img = gray2rgb(rgb_img)
    #                     rgb_img = (rgb_img * (0.5, 0.5, 0.5))
    #
    #                     plt.imshow(rgb_img)
    #                     plt.scatter(
    #                         inside_centroids['x'],
    #                         inside_centroids['y'], s=0.5, ec='none', c='g'
    #                         )
    #
    #                     plt.scatter(
    #                         outside_centroids['x'],
    #                         outside_centroids['y'], s=0.5, ec='none', c='r')
    #
    #                     legend_elements = []
    #                     legend_elements.append(
    #                         Line2D([0], [0], marker='o', color='none',
    #                                label='inside',
    #                                markerfacecolor='g',
    #                                markeredgecolor='none', lw=0.001,
    #                                markersize=6)
    #                                )
    #                     legend_elements.append(
    #                         Line2D([0], [0], marker='o', color='none',
    #                                label='outside',
    #                                markerfacecolor='r',
    #                                markeredgecolor='none', lw=0.001,
    #                                markersize=6)
    #                                )
    #
    #                     plt.legend(
    #                         handles=legend_elements, prop={'size': 6},
    #                         bbox_to_anchor=[1.0, 1.0])
    #
    #                     ratio_real = str(
    #                         round(
    #                             len(inside_centroids)/len(total_centroids), 2)
    #                             )
    #                     total_cells = str(len(total_centroids))
    #                     title_statement = (
    #                         f'{sample}-{celltype}, inside/total: {ratio_real},'
    #                         + f' total cells={total_cells}'
    #                         )
    #                     plt.title(title_statement, fontsize=10)
    #
    #                     plt.grid(False)
    #                     plt.savefig(
    #                         os.path.join(
    #                             spatial_dir,
    #                             f'{protein}_{celltype}_{sample}.png'), dpi=800)
    #                     plt.close('all')
    #
    #                     plt.imshow(rgb_img)
    #
    #                     plt.scatter(
    #                         inside_poisson_points['x'],
    #                         inside_poisson_points['y'],
    #                         s=0.5, ec='none', c='g'
    #                         )
    #
    #                     plt.scatter(
    #                         outside_poisson_points['x'],
    #                         outside_poisson_points['y'],
    #                         s=0.5, ec='none', c='r'
    #                         )
    #
    #                     legend_elements = []
    #                     legend_elements.append(
    #                         Line2D([0], [0], marker='o', color='none',
    #                                label='inside',
    #                                markerfacecolor='g',
    #                                markeredgecolor='none', lw=0.001,
    #                                markersize=6)
    #                                )
    #                     legend_elements.append(
    #                         Line2D([0], [0], marker='o', color='none',
    #                                label='outside',
    #                                markerfacecolor='r',
    #                                markeredgecolor='none', lw=0.001,
    #                                markersize=6)
    #                                )
    #
    #                     plt.legend(
    #                         handles=legend_elements, prop={'size': 6},
    #                         bbox_to_anchor=[1.0, 1.0])
    #
    #                     if not len(total_poisson_points) == 0:
    #                         ratio_fake = str(
    #                             round(
    #                                 len(inside_poisson_points) /
    #                                 len(total_poisson_points), 2))
    #                     else:
    #                         ratio_fake = 'divide by zero'
    #
    #                     poisson_points = str(len(total_poisson_points))
    #                     title_statement = (
    #                         f'{sample}-{celltype}_Poisson-disc,' +
    #                         f' inside/total: {ratio_fake},' +
    #                         f' total points={poisson_points}'
    #                         )
    #                     plt.title(title_statement, fontsize=10)
    #
    #                     plt.grid(False)
    #                     plt.savefig(
    #                         os.path.join(
    #                             spatial_dir,
    #                             f'{protein}_{celltype}_Poisson-disc_' +
    #                             f'{sample}.png'),
    #                         dpi=800
    #                         )
    #                     plt.close('all')
    #
    #                     stats.loc[stats_row_idx] = (
    #                         protein, sample, celltype, float(ratio_real),
    #                         len(total_centroids), float(ratio_fake),
    #                         len(total_poisson_points)
    #                         )
    #                     stats_row_idx += 1
    #
    #     stats.to_csv(os.path.join(spatial_dir, 'stats.csv'))
    #     # stats = pd.read_csv(
    #     #     os.path.join(spatial_dir, 'stats.csv'), index_col=0)
    #
    #     stats = stats[stats['ratio_f'] != 'divide by zero']
    #     stats['ratio_f'] = [float(i) for i in stats['ratio_f']]
    #
    #     stats['condition'] = [
    #         re.findall(r"[^\W\d_]+|\d+", i)[0] for
    #         i in stats['sample']
    #         ]
    #
    #     simulation = stats.copy()
    #
    #     stats.rename(
    #         columns={
    #             'centroids': '# measured cells',
    #             'points': '# simulated cells'},
    #         inplace=True
    #         )
    #
    #     sns.set_style('whitegrid')
    #     sns.scatterplot(
    #         x=stats['# simulated cells'],
    #         y=stats['# measured cells'],
    #         color='k',
    #         data=stats
    #         )
    #     plt.savefig(
    #         os.path.join(
    #             spatial_dir, 'points_v_centroids.pdf'),
    #         dpi=600, bbox_inches='tight')
    #     plt.close('all')
    #
    #     stats[
    #         '% cells overlapping immunomarker (normalized to simulation)'] = (
    #         stats['ratio_r'] - stats['ratio_f']
    #         )
    #
    #     stats.sort_values(
    #         by=['protein', 'celltype', 'condition'],
    #         ascending=[True, False, True],
    #         inplace=True
    #         )
    #     stats.reset_index(drop=True, inplace=True)
    #
    #     stats['label'] = (
    #         stats['protein'] + '_' +
    #         stats['celltype']
    #         )
    #
    #     sns.set_style('whitegrid')
    #     sns.swarmplot(
    #         x='label',
    #         y='% cells overlapping immunomarker (normalized to simulation)',
    #         hue='condition',
    #         size=3,
    #         dodge=True,
    #         palette=['lightgray', 'firebrick'],
    #         data=stats
    #         )
    #     plt.xticks(rotation=90, size=5)
    #
    #     plt.savefig(
    #         os.path.join(
    #             spatial_dir,
    #             'percent cells overlapping immunomarker' +
    #             ' (normalized to simulation).pdf'),
    #         dpi=600, bbox_inches='tight')
    #     plt.close('all')
    #
    #     condition_sig = {}
    #     for name, group in stats.groupby(['protein', 'celltype']):
    #         a = group[
    #             '% cells overlapping immunomarker (normalized to simulation)'][
    #                 group['condition'] == 'cd']
    #         b = group[
    #             '% cells overlapping immunomarker (normalized to simulation)'][
    #                 group['condition'] == 'hfd']
    #         t, pval = ttest_ind(
    #             a, b, axis=0, equal_var=True, nan_policy='propagate')
    #         condition_sig[f'{name}'] = round(pval, 3)
    #
    #     condition_sig_df = pd.DataFrame(condition_sig, index=range(0, 27)).T[0]
    #     condition_sig_df.rename('pval', inplace=True)
    #     condition_sig_df.to_csv(
    #         os.path.join(spatial_dir, 'treatment_sig.csv'), header=True)
    #
    #     stats.to_csv(
    #         os.path.join(spatial_dir, 'stats_normalized.csv'), index=False
    #         )
    #
    #     simulation['dtype'] = 'simulation'
    #
    #     measurement = simulation.copy()
    #     measurement['dtype'] = 'measurement'
    #
    #     cols = [
    #         col for col in simulation.columns if
    #         col not in ['ratio_r', 'centroids']
    #         ]
    #     simulation = simulation[cols]
    #     simulation.rename(
    #         columns={'ratio_f': '% cells overlapping immunomarker'},
    #         inplace=True)
    #
    #     cols = [
    #         col for col in measurement.columns if
    #         col not in ['ratio_f', 'points']
    #         ]
    #     measurement = measurement[cols]
    #
    #     measurement.rename(
    #         columns={
    #             'ratio_r': '% cells overlapping immunomarker',
    #             'centroids': 'points'}, inplace=True
    #             )
    #
    #     stats = simulation.append(measurement, sort=True, ignore_index=True)
    #     del simulation
    #     del measurement
    #
    #     stats['% cells overlapping immunomarker'] = (
    #         stats['% cells overlapping immunomarker'] * 100
    #         )
    #
    #     stats.sort_values(
    #         by=['protein', 'celltype', 'condition', 'dtype'],
    #         ascending=[True, False, True, False],
    #         inplace=True
    #         )
    #     stats.reset_index(drop=True, inplace=True)
    #
    #     stats['label'] = (
    #         stats['protein'] + '_' +
    #         stats['celltype'] + '_' +
    #         stats['condition'] + '_' +
    #         stats['dtype']
    #         )
    #
    #     cols = [
    #         'protein', 'celltype', 'condition', 'dtype',
    #         '% cells overlapping immunomarker', 'points',
    #         'sample', 'label'
    #         ]
    #
    #     stats = stats[cols]
    #
    #     sns.set_style('whitegrid')
    #     sns.swarmplot(
    #         x='label',
    #         y='% cells overlapping immunomarker',
    #         hue='dtype',
    #         size=3,
    #         dodge=False,
    #         palette=['lightgray', 'firebrick'],
    #         data=stats
    #         )
    #     plt.xticks(rotation=90, size=5)
    #
    #     plt.savefig(
    #         os.path.join(
    #             spatial_dir, 'percent cells overlapping immunomarker.pdf'),
    #         dpi=600, bbox_inches='tight')
    #     plt.close('all')
    #
    #     model_sig = {}
    #     for name, group in stats.groupby(['protein', 'celltype', 'condition']):
    #         a = group[
    #             '% cells overlapping immunomarker'][
    #                 group['dtype'] == 'simulation']
    #         b = group[
    #             '% cells overlapping immunomarker'][
    #                 group['dtype'] == 'measurement']
    #         t, pval = ttest_ind(
    #             a, b, axis=0, equal_var=True, nan_policy='propagate')
    #         model_sig[f'{name}'] = round(pval, 3)
    #
    #     model_sig_df = pd.DataFrame(model_sig, index=range(0, 27)).T[0]
    #     model_sig_df.rename('pval', inplace=True)
    #     model_sig_df.to_csv(
    #         os.path.join(spatial_dir, 'model_sig.csv'), header=True
    #         )
    #
    #     stats.to_csv(os.path.join(
    #         spatial_dir, 'stats_unnormalized.csv'), index=False
    #         )
    #     print()
    #     # self.data = df
    #     return data
   #     return data
