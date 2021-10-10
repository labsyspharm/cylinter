import logging
import functools

import os
import re
import glob
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
from skimage import img_as_ubyte

from matplotlib.lines import Line2D
from matplotlib.widgets import Slider, Button
from matplotlib.widgets import TextBox
from matplotlib.colors import ListedColormap
from matplotlib.path import Path
import matplotlib.patheffects as path_effects
from matplotlib.patches import Ellipse

import napari
from tifffile import imread
from tifffile import TiffFile
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

from .utils import (
    dataset_files, log_banner, log_multiline,
    SelectFromCollection, read_dataframe, save_dataframe, read_markers,
    marker_channel_number, categorical_cmap, cluster_expression, clearRAM,
    single_channel_pyramid, matplotlib_warnings, napari_warnings,
    fdrcorrection
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
                 inDir=None,
                 outDir=None,
                 randomSampleSize=None,
                 startModule=None,
                 sampleNames=None,
                 sampleConditions=None,
                 sampleConditionAbbrs=None,
                 sampleStatuses=None,
                 sampleReplicates=None,
                 samplesToExclude=None,
                 markersToExclude=None,

                 # setContrast -
                 viewSample=None,

                 # selectROIs -
                 delintMode=None,
                 showAbChannels=None,

                 # crossCycleCorrelation -
                 yAxisGating=None,
                 logRatioRnge=None,

                 # pruneOutliers -
                 hexbins=None,
                 hexbinGridSize=None,

                 # metQC -
                 metaQC=None,
                 cleanReclassCutoff=None,
                 noisyReclassCutoff=None,

                 # PCA module —
                 channelExclusionsPCA=None,
                 samplesToRemovePCA=None,
                 dimensionPCA=None,
                 pointSize=None,
                 labelPoints=None,
                 distanceCutoff=None,
                 conditionsToSilhouette=None,

                 # clustering module —
                 embeddingAlgorithmQC=None,
                 embeddingAlgorithm=None,
                 channelExclusionsClusteringQC=None,
                 channelExclusionsClustering=None,
                 samplesToRemoveClusteringQC=None,
                 samplesToRemoveClustering=None,
                 normalizeTissueCounts=None,
                 fracForEmbeddingQC=None,
                 fracForEmbedding=None,
                 dimensionEmbeddingQC=None,
                 dimensionEmbedding=None,

                 perplexityQC=None,
                 perplexity=None,
                 earlyExaggerationQC=None,
                 earlyExaggeration=None,
                 learningRateTSNEQC=None,
                 learningRateTSNE=None,
                 metricQC=None,
                 metric=None,
                 randomStateQC=None,
                 randomState=None,

                 nNeighborsQC=None,
                 nNeighbors=None,
                 learningRateUMAPQC=None,
                 learningRateUMAP=None,
                 minDistQC=None,
                 minDist=None,
                 repulsionStrengthQC=None,
                 repulsionStrength=None,

                 # frequencyStats —
                 controlGroups=None,
                 denominatorCluster=None,
                 FDRCorrection=None,

                 # curateThumbnails —
                 numThumbnails=None,
                 squareWindowDimension=None,
                 segOutlines=None,

                 # textbox switch —
                 tbEntry=None,

                 # reclassified data chunks
                 reclassClean=None,
                 reclassNoisy=None,
                 ):

        # assert(SOMETHING)  # placeholder

        self.inDir = inDir
        self.outDir = outDir
        self.randomSampleSize = randomSampleSize
        self.startModule = startModule
        self.sampleNames = sampleNames
        self.sampleConditions = sampleConditions
        self.sampleConditionAbbrs = sampleConditionAbbrs
        self.sampleStatuses = sampleStatuses
        self.sampleReplicates = sampleReplicates
        self.samplesToExclude = samplesToExclude
        self.markersToExclude = markersToExclude

        self.viewSample = viewSample

        self.delintMode = delintMode
        self.showAbChannels = showAbChannels

        self.yAxisGating = yAxisGating
        self.logRatioRnge = logRatioRnge

        self.hexbins = hexbins
        self.hexbinGridSize = hexbinGridSize

        self.metaQC = metaQC
        self.cleanReclassCutoff = cleanReclassCutoff
        self.noisyReclassCutoff = noisyReclassCutoff

        self.channelExclusionsPCA = channelExclusionsPCA
        self.samplesToRemovePCA = samplesToRemovePCA
        self.dimensionPCA = dimensionPCA
        self.pointSize = pointSize
        self.labelPoints = labelPoints
        self.distanceCutoff = distanceCutoff
        self.conditionsToSilhouette = conditionsToSilhouette

        self.embeddingAlgorithmQC = embeddingAlgorithmQC
        self.embeddingAlgorithm = embeddingAlgorithm
        self.channelExclusionsClusteringQC = channelExclusionsClusteringQC
        self.channelExclusionsClustering = channelExclusionsClustering
        self.samplesToRemoveClusteringQC = samplesToRemoveClusteringQC
        self.samplesToRemoveClustering = samplesToRemoveClustering
        self.normalizeTissueCounts = normalizeTissueCounts
        self.fracForEmbeddingQC = fracForEmbeddingQC
        self.fracForEmbedding = fracForEmbedding
        self.dimensionEmbeddingQC = dimensionEmbeddingQC
        self.dimensionEmbedding = dimensionEmbedding

        self.perplexityQC = perplexityQC
        self.perplexity = perplexity
        self.earlyExaggerationQC = earlyExaggerationQC
        self.earlyExaggeration = earlyExaggeration
        self.learningRateTSNEQC = learningRateTSNEQC
        self.learningRateTSNE = learningRateTSNE
        self.metricQC = metricQC
        self.metric = metric
        self.randomStateQC = randomStateQC
        self.randomState = randomState

        self.nNeighborsQC = nNeighborsQC
        self.nNeighbors = nNeighbors
        self.learningRateUMAPQC = learningRateUMAPQC
        self.learningRateUMAP = learningRateUMAP
        self.minDistQC = minDistQC
        self.minDist = minDist
        self.repulsionStrengthQC = repulsionStrengthQC
        self.repulsionStrength = repulsionStrength

        self.controlGroups = controlGroups
        self.denominatorCluster = denominatorCluster
        self.FDRCorrection = FDRCorrection

        self.numThumbnails = numThumbnails
        self.squareWindowDimension = squareWindowDimension
        self.segOutlines = segOutlines

        self.tbEntry = tbEntry

        self.reclassClean = reclassClean
        self.reclassNoisy = reclassNoisy

    @module
    def aggregateData(data, self, args):

        files = natsorted(dataset_files(f'{self.inDir}/csv'))

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            )

        df_list = []
        raw_sample_names_dict = {}
        channel_setlist = []
        for file in files:
            if not file.startswith('.'):

                file_name = file.split('.csv')[0]

                # REMOVE: FOR UPDATING CSV FILE NAMES FOR COMPATIBILITY #####
                # csv = pd.read_csv(
                #     os.path.join(f'{self.inDir}/csv', file)
                #     )
                # keys = [i for i in csv.columns]
                # vals = [i.split('_cellMask')[0] for i in csv.columns]
                # mydict = dict(zip(keys, vals))
                # csv.rename(columns=mydict, inplace=True)
                # csv.to_csv(os.path.join(f'{self.inDir}/csv2', file))

                #################################################
                # disregard samples specified in
                # "samplesToExclude" config parameter
                sample_name = self.sampleNames[file_name]
                if sample_name not in self.samplesToExclude:

                    print(f'IMPORTING sample {sample_name}')

                    csv = pd.read_csv(os.path.join(f'{self.inDir}/csv', file))

                    # drop markers specified in
                    # "markersToExclude" config param
                    csv.drop(
                        columns=[i for i in self.markersToExclude
                                 if i in csv.columns], axis=1, inplace=True)

                    # select boilerplate columns and use specific
                    # mask quantifications for different antibodies
                    # mask_dict = {
                    #     'Hoechst0': 'nucleiRingMask',
                    #     'Hoechst1': 'nucleiRingMask',
                    #     'Hoechst2': 'nucleiRingMask',
                    #     'anti_CD3': 'cytoRingMask',
                    #     'anti_CD45RO': 'cytoRingMask',
                    #     'Hoechst3': 'nucleiRingMask',
                    #     'Keratin_570': 'cellRingMask',
                    #     'aSMA_660': 'cellRingMask',
                    #     'Hoechst4': 'nucleiRingMask',
                    #     'CD4_488': 'cytoRingMask',
                    #     'CD45_PE': 'cytoRingMask',
                    #     'PD1_647': 'cytoRingMask',
                    #     'Hoechst5': 'nucleiRingMask',
                    #     'CD20_488': 'cytoRingMask',
                    #     'CD68_555': 'cellRingMask',
                    #     'CD8a_660': 'cytoRingMask',
                    #     'Hoechst6': 'nucleiRingMask',
                    #     'CD163_488': 'cellRingMask',
                    #     'FOXP3_570': 'nucleiRingMask',
                    #     'PDL1_647': 'cytoRingMask',
                    #     'Hoechst7': 'nucleiRingMask',
                    #     'Ecad_488': 'cellRingMask',
                    #     'Vimentin_555': 'cellRingMask',
                    #     'CDX2_647': 'cellRingMask',
                    #     'Hoechst8': 'nucleiRingMask',
                    #     'LaminABC_488': 'nucleiRingMask',
                    #     'Desmin_555': 'cellRingMask',
                    #     'CD31_647': 'nucleiRingMask',
                    #     'Hoechst9': 'nucleiRingMask',
                    #     'PCNA_488': 'nucleiRingMask',
                    #     'CollagenIV_647': 'cellRingMask'
                    #     }
                    # mask_object_cols = (
                    #     ['CellID', 'X_centroid', 'Y_centroid', 'Area',
                    #      'MajorAxisLength', 'MinorAxisLength',
                    #      'Eccentricity', 'Solidity', 'Extent',
                    #      'Orientation'] +
                    #     [f'{i}_{mask_dict[i]}' for i
                    #      in markers['marker_name']]
                    #      )

                    # select boilerplate columns
                    # and drop mask object columns not specified by
                    # the "maskObject" parameter in config.yml
                    cols = (
                        ['CellID', 'X_centroid', 'Y_centroid', 'Area',
                         'MajorAxisLength', 'MinorAxisLength',
                         'Eccentricity', 'Solidity', 'Extent',
                         'Orientation'] +
                        [i for i in markers['marker_name']]
                         )

                    csv = csv[cols]

                    # remap mask objects to a common name for convenience
                    # use with mask object dict above
                    # csv.columns = (
                    #     ['CellID', 'X_centroid', 'Y_centroid', 'Area',
                    #      'MajorAxisLength', 'MinorAxisLength',
                    #      'Eccentricity', 'Solidity', 'Extent',
                    #      'Orientation'] +
                    #     [i for i in markers['marker_name']]
                    #      )

                    # add sample column
                    csv['Sample'] = sample_name

                    # add condition column
                    csv['Condition'] = self.sampleConditions[file_name]

                    # add replicate column
                    csv['Replicate'] = self.sampleReplicates[file_name]

                    # append dataframe to list
                    df_list.append(csv)

                    # append the set of csv columns for sample to a list
                    # this will be used to select columns shared among samples
                    channel_setlist.append(set(csv.columns))
                else:
                    print(f'censoring sample {sample_name}')

        print()

        # stack dataframes row-wise
        data = pd.concat(df_list, axis=0)
        del df_list

        # organize columns
        cols = (
            ['CellID', 'Sample', 'Condition', 'Replicate', 'X_centroid',
             'Y_centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength',
             'Eccentricity', 'Solidity', 'Extent', 'Orientation'] +
            [i for i in markers['marker_name']]
             )
        data = data[cols]

        # ONLY SELECT CHANNELS COMMON TO ALL SAMPLES
        channels_set = list(set.intersection(*channel_setlist))

        print(f'{len(data.columns)} total features')
        print(f'{len(channels_set)} features in common between all samples')

        before = set(data.columns)
        after = set(channels_set)
        if len(before.difference(after)) == 0:
            pass
        else:
            markers_to_drop = list(before.difference(after))
            print(
                f'Features {markers_to_drop} were not probed in all' +
                ' samples and will be dropped from the analysis.'
                )
        data = data[channels_set].copy()

        # perform data subsetting
        data = data.sample(frac=self.randomSampleSize, random_state=1)
        data.sort_values(by=['Sample', 'CellID'], inplace=True)

        # assign global index
        data.reset_index(drop=True, inplace=True)

        print()
        print()
        return data

    @module
    def selectROIs(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            )

        # create ROIs directory if it doesn't already exist
        roi_dir = os.path.join(self.outDir, 'ROIs')
        if not os.path.exists(roi_dir):
            os.makedirs(roi_dir)

        # load polygon dictionary if it exists,
        # identify samples in need of ROI selection
        if os.path.exists(os.path.join(roi_dir, 'polygon_dict.pkl')):
            f = open(os.path.join(roi_dir, 'polygon_dict.pkl'), 'rb')
            polygon_dict = pickle.load(f)
            completed_samples = set(polygon_dict.keys())
            total_samples = set(data['Sample'].unique())
            samples_to_draw = total_samples.difference(completed_samples)
            print(f'Samples requiring ROI selection: {len(samples_to_draw)}')

        # if polygon dictionary does not exist, create it
        else:
            samples_to_draw = data['Sample'].unique()
            print(f'Samples requiring ROI selection: {len(samples_to_draw)}')
            polygon_dict = {}

        # if the polygon dict doesn't have vertices for all samples
        # or their are samples which don't have a corresponding text file
        # with cell IDs of selected cells,
        # loop over remaining samples
        if (
          (len(samples_to_draw) > 0) or
          (len([name for name in os.listdir(roi_dir) if
           name.endswith('.txt')]) < len(data['Sample'].unique()))):

            napari_warnings()

            for sample_name in natsorted(samples_to_draw):

                polygons = []

                if self.showAbChannels:
                    for e, ch in enumerate(reversed(abx_channels)):
                        channel_number = marker_channel_number(markers, ch)

                        # read antibody image
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{sample_name}*.tif'):
                            img = single_channel_pyramid(
                                file_path, channel=channel_number.item() - 1)

                        if e == 0:
                            viewer = napari.view_image(
                                img, rgb=False, blending='additive',
                                colormap='green', visible=False,
                                name=ch
                                )

                        else:
                            viewer.add_image(
                                img, rgb=False, blending='additive',
                                colormap='green', visible=False,
                                name=ch
                                )

                # read segmentation outlines
                file_path = f'{self.inDir}/seg/{sample_name}.ome.tif'
                seg = single_channel_pyramid(
                    glob.glob(file_path)[0], channel=0)

                viewer.add_image(
                    seg, rgb=False, blending='additive',
                    opacity=0.5, colormap='red', visible=False,
                    name='segmentation')

                # read DNA1 channel
                for file_path in glob.glob(
                  f'{self.inDir}/tif/{sample_name}*.tif'):
                    dna = single_channel_pyramid(file_path, channel=0)

                viewer.add_image(
                    dna, rgb=False, blending='additive',
                    colormap='gray', visible=True,
                    name=f'{dna1}: {sample_name}'
                    )

                selection_layer = viewer.add_shapes(
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

                napari.run()

                # store lists vertices per sample as a dictionary
                for roi in selection_layer.data:
                    polygons.append((selection_layer.shape_type[0], roi))

                polygon_dict[sample_name] = polygons

                os.chdir(roi_dir)
                f = open(os.path.join(roi_dir, 'polygon_dict.pkl'), 'wb')
                pickle.dump(polygon_dict, f)
                f.close()
            print()

            # make a set of samples requiring cell ID selection
            os.chdir(roi_dir)
            samples_for_cell_selection = set(
                data['Sample'].unique()).difference(
                    set([i.split('.txt')[0] for i in os.listdir()
                         if i.endswith('.txt')]))

            # create txt files per sample with cell IDs to drop
            for sample_name in natsorted(samples_for_cell_selection):

                print(f'Applying ROIs to sample {sample_name}')

                sample_data = data[['X_centroid', 'Y_centroid', 'CellID']][
                    data['Sample'] == sample_name].astype(int)

                sample_data['tuple'] = list(
                    zip(sample_data['Y_centroid'],
                        sample_data['X_centroid'])
                    )

                for file_path in glob.glob(
                  f'{self.inDir}/tif/{sample_name}*.tif'):
                    dna = imread(file_path, key=0)

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

                cell_ids_to_drop = set()
                mask_coords = set()
                if polygon_dict[sample_name]:
                    for shape_type, verts in polygon_dict[sample_name]:

                        selection_verts = np.round(verts).astype(int)

                        if shape_type == 'ellipse':
                            col_min = selection_verts[0][1]
                            col_max = selection_verts[1][1]
                            row_min = selection_verts[0][0]
                            row_max = selection_verts[2][0]
                            col_center = col_min + ((col_max-col_min)/2)
                            row_center = row_min + ((row_max-row_min)/2)
                            ellipse = Ellipse(
                                (row_center, col_center),
                                width=(row_max-row_min),
                                height=(col_max-col_min)
                                )
                            grid = ellipse.contains_points(pixel_coords)
                            mask = grid.reshape(
                                dna.shape[0], dna.shape[1])
                            mask_coords.update(
                                [tuple(i) for i in np.argwhere(mask)]
                                )
                        else:
                            polygon = Path(selection_verts)
                            grid = polygon.contains_points(pixel_coords)
                            mask = grid.reshape(
                                dna.shape[0], dna.shape[1])
                            mask_coords.update(
                                [tuple(i) for i in np.argwhere(mask)]
                                )

                    clearRAM(print_usage=False)
                    del grid, mask, dna, pixel_coords
                    clearRAM(print_usage=False)

                inter = mask_coords.intersection(cell_coords)

                if self.delintMode is True:
                    cell_ids_to_drop.update(
                        [i[1]['CellID'] for i in sample_data.iterrows() if
                         i[1]['tuple'] in inter]
                         )
                else:
                    if polygon_dict[sample_name]:
                        cell_ids_to_drop.update(
                            [i[1]['CellID'] for i in sample_data.iterrows()
                             if i[1]['tuple'] not in inter]
                             )
                    # if not polygons selected for a sample,
                    # no cell IDs will be dropped

                clearRAM(print_usage=False)
                del sample_data, inter, cell_coords
                clearRAM(print_usage=False)

                # pickle and save selected cell IDs for the sample
                os.chdir(roi_dir)
                f = open(f'{sample_name}.txt', 'wb')
                pickle.dump(cell_ids_to_drop, f)
                f.close()
            print()

            # drop cells from samples
            for file in natsorted(os.listdir(roi_dir)):
                if file.endswith('.txt'):
                    file_name = file.split('.txt')[0]
                    print(f'Dropping cells from sample {file_name}')
                    f = open(os.path.join(roi_dir, file), 'rb')
                    cell_ids = pickle.load(f)
                    cell_ids = set(cell_ids)
                    global_idxs_to_drop = data[
                        (data['Sample'] == file_name)
                        & (data['CellID'].isin(cell_ids))
                        ].index
                    data.drop(global_idxs_to_drop, inplace=True)

        else:
            print()

            # drop cells from samples
            for file in natsorted(os.listdir(roi_dir)):
                if file.endswith('.txt'):
                    file_name = file.split('.txt')[0]
                    print(f'Dropping cells from {file_name}')
                    f = open(os.path.join(roi_dir, file), 'rb')
                    cell_ids = pickle.load(f)
                    cell_ids = set(cell_ids)
                    global_idxs_to_drop = data[
                        (data['Sample'] == file_name)
                        & (data['CellID'].isin(cell_ids))
                        ].index
                    data.drop(global_idxs_to_drop, inplace=True)
        print()

        # save images of tissue with selected data points
        image_dir = os.path.join(roi_dir, 'images')
        if not os.path.exists(image_dir):
            os.mkdir(image_dir)

        for sample_name in natsorted(data['Sample'].unique()):

            print(f'Plotting ROI selections for sample {sample_name}')

            for file_path in glob.glob(
              f'{self.inDir}/tif/{sample_name}*.tif'):
                dna = imread(file_path, key=0)

            fig, ax = plt.subplots()
            ax.imshow(dna, cmap='gray')
            ax.grid(False)
            coords = data[['X_centroid', 'Y_centroid', 'Area']][
                data['Sample'] == sample_name]
            sp = ax.scatter(
                coords['X_centroid'], coords['Y_centroid'],
                s=0.35, lw=0.0,
                c=coords['Area'], cmap='viridis'
                )
            plt.title(
                f'Sample {sample_name}. ' +
                'Selected cells colored by segmentation area')
            plt.colorbar(sp)
            plt.savefig(
                os.path.join(
                    image_dir, f'{sample_name}.png'), dpi=1000)
            plt.close('all')

        print()
        print()
        return data

    @module
    def intensityFilter(data, self, args):

        napari_warnings()

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            )

        # create intensity directory if it doesn't already exist
        intensity_dir = os.path.join(self.outDir, 'intensity')
        if not os.path.exists(intensity_dir):
            os.makedirs(intensity_dir)

        # set histogram bin size
        bins = 100

        # set histogram type
        histtype = 'stepfilled'

        # pick up where samples-loop left off
        if os.path.exists(
          os.path.join(intensity_dir, 'idxs_to_drop.pkl')):

            # read dictionary of sample indices to drop
            f = open(os.path.join(intensity_dir, 'idxs_to_drop.pkl'), 'rb')
            idxs_to_drop = pickle.load(f)

            # get names of samples already run
            previously_run_samples = idxs_to_drop.keys()

            # get names of samples remaining
            samples_remaining = (
                len(data['Sample'].unique())
                - len(previously_run_samples)
                )

            print(f'Samples to threshold: {samples_remaining}')

            # drop samples previously run
            df = data[~data['Sample'].isin(previously_run_samples)].copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique())
                )
            df.sort_values('Sample', inplace=True)
            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        else:
            # get names of samples remaining (total samples in this case)
            samples_remaining = len(data['Sample'].unique())

            print(f'Samples to threshold: {samples_remaining}')

            # initialize dictionary of sample indices to drop
            idxs_to_drop = {}

            # drop samples previously run (total samples in this case)
            df = data.copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique())
                )
            df.sort_values('Sample', inplace=True)
            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        # loop over samples
        for name, group in df.groupby('Sample', sort=False):

            sns.set_style('whitegrid')
            fig, ax = plt.subplots()
            matplotlib_warnings(fig)

            plt.subplots_adjust(left=0.25, bottom=0.25)

            n, bins, patches = plt.hist(
                group[dna1], bins=bins,
                density=False, color='grey', ec='none',
                alpha=0.75, histtype=histtype,
                range=None, label='before'
                )

            plt.title(
                f'Sample={name}  mean DNA intensity', size=10)

            plt.ylabel('count')

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

                # remove current lines
                [i.remove() for i in ax.get_lines()]

                # new cutoffs
                lowerCutoff = sLower.val
                upperCutoff = sUpper.val

                # update plot with cutoffs
                blueLine = ax.axvline(
                    x=lowerCutoff, c='b', linewidth=2.5)
                redLine = ax.axvline(
                    x=upperCutoff, c='r', linewidth=2.5)

                return lowerCutoff, upperCutoff

            sLower.on_changed(update)
            sUpper.on_changed(update)

            resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
            button = Button(
                resetax, 'Reset', color=axcolor, hovercolor='0.975'
                    )

            def reset(event):
                sLower.reset()
                sUpper.reset()

            button.on_clicked(reset)

            def submit(text):

                if text == self.tbEntry:
                    # reset tbEntry
                    self.tbEntry = None
                else:
                    # assign current textbox entry to tbEntry
                    # if text and tbEntry are different
                    self.tbEntry = text

                    # get current cutoffs
                    lowerCutoff, upperCutoff = update(val=None)

                    # apply lower and upper cutoffs
                    group_update = group[
                        (group[dna1] > lowerCutoff) &
                        (group[dna1] < upperCutoff)
                        ]

                    if text == group['Sample'].unique():

                        # read DNA1 channel
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{text}*.tif'):
                            dna = single_channel_pyramid(file_path, channel=0)

                        # read segmentation outlines
                        file_path = f'{self.inDir}/seg/{text}*.ome.tif'

                        seg = single_channel_pyramid(
                            glob.glob(file_path)[0], channel=0
                            )

                        centroids = group_update[['Y_centroid', 'X_centroid']]

                        dna_intensity = group_update[dna1].values
                        point_properties = {'dna_intensity': dna_intensity}

                        viewer = napari.view_image(
                            dna, rgb=False, name=f'{dna1}: {text}')

                        viewer.add_image(
                            seg, rgb=False, blending='additive',
                            opacity=0.5, colormap='red', visible=False,
                            name='segmentation'
                            )

                        viewer.add_points(
                            centroids, name='Intensity',
                            properties=point_properties,
                            face_color='dna_intensity',
                            face_colormap='viridis',
                            edge_width=0.0, size=4.0
                            )

                        napari.run()

                    else:
                        print('Must enter name of current sample.')

            axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
            text_box = TextBox(
                axbox, 'evaluation sample name', initial='',
                color='0.95',
                hovercolor='1.0',
                label_pad=0.05
                )

            text_box.on_submit(submit)

            plt.show(block=True)

            # ensure tbEntry is set to default (None) for future use
            self.tbEntry = None

            lowerCutoff, upperCutoff = update(val=None)

            # plot DNA intensity histogram BEFORE filtering
            fig, ax = plt.subplots()
            plt.hist(
                group[dna1], bins=bins,
                density=False, color='b', ec='none',
                alpha=0.5, histtype=histtype,
                range=None, label='before'
                )

            if lowerCutoff == upperCutoff:
                lowerCutoff = group[dna1].min()
                upperCutoff = group[dna1].max()

            # apply lower and upper cutoffs
            group_update = group[
                (group[dna1] > lowerCutoff) &
                (group[dna1] < upperCutoff)
                ]

            # plot DNA intensity histogram AFTER filtering
            plt.hist(
                group_update[dna1], bins=bins, color='r', ec='none',
                alpha=0.5, histtype=histtype, range=None,
                label='after')
            plt.xlabel('mean DNA intensity')
            plt.ylabel('count')
            plt.title(
                f'Sample={name}  mean DNA intensity)',
                size=10
                )

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
            plt.savefig(
                os.path.join(intensity_dir, f'{name}.pdf')
                )
            plt.close()

            # isolate sample data to drop
            data_to_drop = group.copy()[
                (group[dna1] < lowerCutoff) |
                (group[dna1] > upperCutoff)
                ]

            # create unique IDs for cells to drop in current sample
            data_to_drop['handle'] = (
                data_to_drop['CellID'].map(str) + '_' +
                data_to_drop['Sample']
                )

            # add sample indices to drop to idxs_to_drop dictionary
            idxs_to_drop[name] = [i for i in data_to_drop['handle']]

            # save updated idxs_to_drop dictionary as a pickle
            os.chdir(intensity_dir)
            f = open('idxs_to_drop.pkl', 'wb')
            pickle.dump(idxs_to_drop, f)
            f.close()

        # create unique IDs for cells across all samples
        data['handle'] = data['CellID'].map(str) + '_' + data['Sample']

        # create a single list of indices to drop
        total_indices_to_drop = []
        for k, v in idxs_to_drop.items():
            total_indices_to_drop.extend(v)

        # isolate cells not in total_indices_to_drop
        data = data[~data['handle'].isin(total_indices_to_drop)].copy()

        # drop unique ID column
        data.drop(columns='handle', inplace=True)

        print()
        print()
        return data

    @module
    def areaFilter(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            )

        # create area directory if it doesn't already exist
        area_dir = os.path.join(self.outDir, 'area')
        if not os.path.exists(area_dir):
            os.makedirs(area_dir)

        # set histogram bin size
        bins = 50

        # set histogram type
        histtype = 'stepfilled'

        # pick up where samples-loop left off
        if os.path.exists(
          os.path.join(area_dir, 'idxs_to_drop.pkl')):

            # read dictionary of sample indices to drop
            f = open(os.path.join(area_dir, 'idxs_to_drop.pkl'), 'rb')
            idxs_to_drop = pickle.load(f)

            # get names of samples already run
            previously_run_samples = idxs_to_drop.keys()

            # get names of samples remaining
            samples_remaining = (
                len(data['Sample'].unique())
                - len(previously_run_samples)
                )

            print(f'Samples to threshold: {samples_remaining}')

            # drop samples previously run
            df = data[~data['Sample'].isin(previously_run_samples)].copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique())
                )
            df.sort_values('Sample', inplace=True)
            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        else:
            # get names of samples remaining (total samples in this case)
            samples_remaining = len(data['Sample'].unique())

            print(f'Samples to threshold: {samples_remaining}')

            # initialize dictionary of sample indices to drop
            idxs_to_drop = {}

            # drop samples previously run (total samples in this case)
            df = data.copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique())
                )
            df.sort_values('Sample', inplace=True)
            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        # loop over remaining samples
        for name, group in df.groupby('Sample', sort=False):

            sns.set_style('whitegrid')
            fig, ax = plt.subplots()
            matplotlib_warnings(fig)

            plt.subplots_adjust(left=0.25, bottom=0.25)

            n, bins, patches = plt.hist(
                group['Area'], bins=bins,
                density=False, color='grey', ec='none',
                alpha=0.75, histtype=histtype,
                range=None, label='before'
                )

            plt.title(f'Sample={name}  cell segmentation area', size=10)

            plt.ylabel('count')

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

                # remove current lines
                [i.remove() for i in ax.get_lines()]

                # new cutoffs
                lowerCutoff = sLower.val
                upperCutoff = sUpper.val

                # update plot with cutoffs
                blueLine = ax.axvline(
                    x=lowerCutoff, c='b', linewidth=2.5)
                redLine = ax.axvline(
                    x=upperCutoff, c='r', linewidth=2.5)

                return lowerCutoff, upperCutoff

            sLower.on_changed(update)
            sUpper.on_changed(update)

            resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
            button = Button(
                resetax, 'Reset', color=axcolor, hovercolor='0.975'
                )

            def reset(event):
                sLower.reset()
                sUpper.reset()

            button.on_clicked(reset)

            def submit(text):

                if text == self.tbEntry:
                    # reset tbEntry
                    self.tbEntry = None
                else:
                    # assign current textbox entry to tbEntry
                    # if text and tbEntry are different
                    self.tbEntry = text

                    lowerCutoff, upperCutoff = update(val=None)

                    # apply lower and upper cutoffs
                    group_update = group[
                        (group['Area'] > lowerCutoff) &
                        (group['Area'] < upperCutoff)
                        ]

                    if text in group['Sample'].unique():

                        # read DNA1 channel
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{text}*.tif'):
                            dna = single_channel_pyramid(file_path, channel=0)

                        # read segmentation outlines
                        file_path = f'{self.inDir}/seg/{text}.ome.tif'
                        seg = single_channel_pyramid(
                            glob.glob(file_path)[0], channel=0
                            )

                        centroids = group_update[['Y_centroid', 'X_centroid']]

                        cell_area = group_update['Area'].values
                        point_properties = {'cell_area': cell_area}

                        viewer = napari.view_image(
                            dna, rgb=False, name=dna1
                            )

                        viewer.add_image(
                            seg, rgb=False, blending='additive',
                            opacity=0.5, colormap='red', visible=False,
                            name='segmentation'
                            )

                        viewer.add_points(
                            centroids, name='Area',
                            properties=point_properties,
                            face_color='cell_area',
                            face_colormap='viridis',
                            edge_width=0.0, size=4.0
                            )

                        napari.run()
                    else:
                        print('Must enter name of current sample.')

            axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
            text_box = TextBox(
                axbox, 'evaluation sample name', initial='',
                color='0.95',
                hovercolor='1.0',
                label_pad=0.05
                )

            text_box.on_submit(submit)
            plt.show(block=True)

            # ensure tbEntry is set to default (None) for future use
            self.tbEntry = None

            lowerCutoff, upperCutoff = update(val=None)

            # plot DNA area histogram BEFORE filtering
            fig, ax = plt.subplots()
            plt.hist(
                group['Area'], bins=bins,
                density=False, color='b', ec='none',
                alpha=0.5, histtype=histtype,
                range=None, label='before'
                )

            if lowerCutoff == upperCutoff:
                lowerCutoff = group['Area'].min()
                upperCutoff = group['Area'].max()

            # apply lower and upper cutoffs
            group_update = group[
                (group['Area'] > lowerCutoff) &
                (group['Area'] < upperCutoff)
                ]

            # plot DNA area histogram AFTER filtering
            plt.hist(
                group_update['Area'], bins=bins, color='r', ec='none',
                alpha=0.5, histtype=histtype, range=None, label='after'
                )
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
            plt.savefig(os.path.join(area_dir, f'{name}.pdf'))
            plt.close()

            # isolate sample data to drop
            data_to_drop = group.copy()[
                (group['Area'] < lowerCutoff) |
                (group['Area'] > upperCutoff)
                ]

            # create unique IDs for cells to drop in current sample
            data_to_drop['handle'] = (
                data_to_drop['CellID'].map(str) + '_' +
                data_to_drop['Sample']
                )

            # add sample indices to drop to idxs_to_drop dictionary
            idxs_to_drop[name] = [i for i in data_to_drop['handle']]

            # save updated idxs_to_drop dictionary as a pickle
            os.chdir(area_dir)
            f = open('idxs_to_drop.pkl', 'wb')
            pickle.dump(idxs_to_drop, f)
            f.close()

        # create unique IDs for cells across all samples
        data['handle'] = data['CellID'].map(str) + '_' + data['Sample']

        # create a single list of indices to drop
        total_indices_to_drop = []
        for k, v in idxs_to_drop.items():
            total_indices_to_drop.extend(v)

        # isolate cells not in total_indices_to_drop
        data = data[~data['handle'].isin(total_indices_to_drop)].copy()

        # drop unique ID column
        data.drop(columns='handle', inplace=True)

        print()
        print()
        return data

    @module
    def cycleCorrelation(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            )

        # create cycles directory if it doesn't already exist
        cycles_dir = os.path.join(self.outDir, 'cycles')
        if not os.path.exists(cycles_dir):
            os.makedirs(cycles_dir)

        # get ordered list of DNA cycles
        dna_cycles = natsorted(
            data.columns[data.columns.str.contains(dna_moniker)]
            )

        # compute log(cycle 1/n) ratios
        ratios = pd.DataFrame(
            [np.log10((data[dna1] + 0.00001) /
             (data[i] + 0.00001)) for i in dna_cycles]).T

        # computing ratios changes columns headers,
        # create new ratio column headers
        unnamed_headers = [
            i for i in ratios.columns if i.startswith('Unnamed')
            ]
        ratio_headers = [f'{dna1}/{i}' for i in [j for j in dna_cycles[1:]]]
        ratio_columns = dict(zip(unnamed_headers, ratio_headers))
        ratio_columns[dna1] = f'{dna1}/{dna1}'
        ratios.rename(columns=ratio_columns, inplace=True)
        ratios['sample'] = data['Sample']

        # # compute log(cycle n/n+1) ratios
        # ratios = pd.DataFrame(
        #     [np.log10((data[i] + 0.00001) /
        #      (data[dna_moniker +
        #       str(int([re.findall(r'(\w+?)(\d+)', i)[0]][0][1]) + 1)] +
        #       0.00001)) for i in natsorted(
        #         data.columns[
        #             data.columns.str.contains(dna_moniker)])[0:-1]]).T
        #
        # # computing ratios changes columns headers,
        # # create new ratio column headers
        # ratio_column_headers1 = [i for i in ratios.columns]
        # ratio_column_headers2 = [
        #     f'{dna_moniker}{i}/{dna_moniker}{i+1}'
        #     for i in range(1, len(ratio_column_headers1)+1)]
        # ratio_columns = dict(
        #     zip(ratio_column_headers1, ratio_column_headers2))
        # ratios.rename(columns=ratio_columns, inplace=True)
        # ratios['sample'] = data['Sample']

        # melt ratios dataframe
        ratios_melt = (
            ratios
            .reset_index()
            .melt(id_vars=['sample', 'index'], var_name='cycle',
                  value_name='log10(ratio)')
            )

        # convert sample and cycle columns to ordered categoricals
        # and sort naturally on sample, cycle, and index
        ratios_melt['sample'] = pd.Categorical(
            ratios_melt['sample'], ordered=True,
            categories=natsorted(
                ratios_melt['sample'].unique()))
        ratios_melt['cycle'] = pd.Categorical(
            ratios_melt['cycle'], ordered=True,
            categories=natsorted(
                ratios_melt['cycle'].unique()))
        ratios_melt = ratios_melt.sort_values(['sample', 'cycle', 'index'])
        # convert columns back to strings
        ratios_melt['sample'] = ratios_melt['sample'].astype('str')
        ratios_melt['cycle'] = ratios_melt['cycle'].astype('str')

        # plot log(cycle 1/n) ratio histograms for all samples for evaluation
        # if cycle_correlation(logRatio).png doesn't already exist
        if not os.path.exists(
          os.path.join(cycles_dir, 'cycle_correlation(logRatio).png')):

            sns.set(font_scale=0.5)
            sns.set_style('whitegrid')

            g = sns.FacetGrid(
                ratios_melt, row='sample',
                col='cycle', sharey=False
                )

            g = g.map(
                plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
                ec='none', range=self.logRatioRnge, bins=200, density=True
                )

            plt.savefig(
                os.path.join(cycles_dir, 'cycle_correlation(logRatio).pdf')
                )
            plt.close()

        if self.yAxisGating is True:

            # grab selected y-axis count cutoff if one was entered
            if os.path.exists(
              os.path.join(cycles_dir, 'count_cutoff.pkl')):

                pickle_in = open(os.path.join(
                    cycles_dir, 'count_cutoff.pkl'), 'rb')
                count_cutoff = pickle.load(pickle_in)

            else:

                subprocess.call(
                    ['open', '-a', 'Preview',
                     os.path.join(
                        cycles_dir, 'cycle_correlation(logRatio).pdf')]
                     )

                def submit(text):

                    # save selected y-axis count cutoff as pickle
                    count_cutoff = float(text)
                    os.chdir(cycles_dir)
                    f = open('count_cutoff.pkl', 'wb')
                    pickle.dump(count_cutoff, f)
                    f.close()

                    ###########################################################
                    # visualize log(cycle n/n+1) ratios in Napari

                    # sample_to_inspect = str(text.split(', ')[1])
                    #
                    # if sample_to_inspect in data['Sample'].unique():
                    #     sample_df = data[data['Sample'] == sample_to_inspect]
                    #     sample_centroids = sample_df[
                    #         ['Y_centroid', 'X_centroid']]
                    #
                    #     with napari.gui_qt():
                    #
                    #         # add dna images
                    #         for e, i in enumerate(sorted(markers.index[
                    #           markers['marker_name'].str.contains(
                    #               dna_moniker)], reverse=True)):
                    #
                    #             dna = imread(
                    #                 f'{self.inDir}/tif/' +
                    #                 f'{sample_to_inspect}.*tif',
                    #                 key=i
                    #                 )
                    #
                    #             name = markers.loc[i]['marker_name']
                    #             cycle_num = markers.loc[i]['cycle_number']
                    #
                    #             if name == dna1:
                    #                 visible = True
                    #             else:
                    #                 visible = False
                    #
                    #             if e == 0:
                    #                 viewer = napari.view_image(
                    #                     dna, rgb=False, blending='opaque',
                    #                     visible=visible,
                    #                     name=name
                    #                     )
                    #             else:
                    #                 viewer.add_image(
                    #                     dna, rgb=False, blending='opaque',
                    #                     visible=visible,
                    #                     name=name
                    #                     )
                    #
                    #             # log(cycle n/n+1) ratios
                    #             if cycle_num < markers['cycle_number'].max():
                    #                 sample_ratios = np.log10(
                    #                     (sample_df[
                    #                         f'{name}_{self.maskObject}']
                    #                         + 0.00001)
                    #                     /
                    #                     (sample_df[dna_moniker
                    #                      + str(int(
                    #                         [re.findall(
                    #                             r'(\w+?)(\d+)', name)[0]
                    #                          ][0][1]
                    #                         ) + 1) + '_' + self.maskObject]
                    #                         + 0.00001)
                    #                     )
                    #
                    #             # log(cycle 1/n) ratios
                    #             # if cycle_num != 1:
                    #             #     sample_ratios = np.log10(
                    #             #         (sample_df[
                    #             #             f'{dna1}_{self.maskObject}']
                    #             #             + 0.00001)
                    #             #         /
                    #             #         (sample_df[
                    #             #             f'{name}_{self.maskObject}']
                    #             #             + 0.00001)
                    #             #         )
                    #
                    #                 sample_ratios = np.clip(
                    #                     sample_ratios,
                    #                     a_min=np.percentile(
                    #                         sample_ratios, 1.0),
                    #                     a_max=np.percentile(
                    #                         sample_ratios, 99.0)
                    #                     )
                    #
                    #                 point_properties = {
                    #                     'face_color': sample_ratios
                    #                     }
                    #
                    #                 viewer.add_points(
                    #                     sample_centroids,
                    #                     name=f'log({cycle_num}/' +
                    #                     f'{cycle_num + 1})',
                    #                     visible=False,
                    #                     properties=point_properties,
                    #                     face_color='face_color',
                    #                     face_colormap='PiYG',
                    #                     edge_color='k',
                    #                     edge_width=0.0, size=7.0
                    #                     )
                    #
                    #         # rearrange viwer layers (DNA images first)
                    #         layer_names = [str(i) for i in viewer.layers]
                    #
                    #         current_order = tuple(
                    #             [layer_names.index(i) for i in layer_names]
                    #             )
                    #
                    #         dna_idxs = [
                    #             layer_names.index(i) for i in layer_names
                    #             if dna_moniker in i
                    #             ]
                    #         log_idxs = [
                    #             layer_names.index(i) for i in layer_names
                    #             if 'log(' in i
                    #             ]
                    #
                    #         target_order = tuple(log_idxs + dna_idxs)
                    #
                    #         viewer.layers[current_order] = viewer.layers[
                    #             target_order
                    #             ]
                    ###########################################################

                    # render log(cycle 1/n) histograms with y-axis count cutoff
                    sns.set(font_scale=0.5)
                    sns.set_style('whitegrid')
                    g = sns.FacetGrid(
                        ratios_melt, row='sample',
                        col='cycle', sharey=False
                        )
                    g = g.map(
                        plt.hist, 'log10(ratio)', color='r',
                        histtype='stepfilled', ec='none',
                        range=self.logRatioRnge,
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
                    axbox, 'countCutoff', initial='',
                    color='0.95',
                    hovercolor='1.0',
                    label_pad=0.05
                    )
                text_box.label.set_size(12)

                text_box.on_submit(submit)

                plt.show(block=True)

                # grab selected y-axis count cutoff if one was entered
                if os.path.exists(
                  os.path.join(cycles_dir, 'count_cutoff.pkl')):

                    pickle_in = open(os.path.join(
                        cycles_dir, 'count_cutoff.pkl'), 'rb')
                    count_cutoff = pickle.load(pickle_in)

                else:
                    # save count cutoff as pickle as 0.0
                    count_cutoff = 0.0
                    f = open(os.path.join(
                        cycles_dir, 'count_cutoff.pkl'), 'wb')
                    pickle.dump(count_cutoff, f)
                    f.close()

            # initialize set of cell indices to drop
            indices_to_drop = set()

            # loop over samples and cycles except (cycle1/cycle1)
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
                        count_indices = np.where(counts >= count_cutoff)
                        bin_values = [
                            bins[i] for i in np.append(
                                count_indices[0], count_indices[0].max() + 1)
                                ]

                        if len(bin_values) > 1:
                            min_bin_val = min(bin_values)
                            max_bin_val = max(bin_values)

                            # get indices in log(ratio) series outside
                            # min_bin_val and max_bin_val
                            idxs = list(
                                cycle_data['index'][
                                    (cycle_data['log10(ratio)'] < min_bin_val)
                                    |
                                    (cycle_data['log10(ratio)'] > max_bin_val)]
                                    )

                            # append indices of uncorrelated
                            # log(ratios) to idx_list
                            indices_to_drop.update(set(idxs))
            print()

            # filter dataframe by selecting indices NOT in the
            # indices_to_drop list
            print('Y-axis gating: Dropping unstable cells from all samples...')
            df = data.loc[~data.index.isin(indices_to_drop)]
            plt.close('all')

        elif self.yAxisGating is False:

            # pick up where samples loop left off
            if os.path.exists(
              os.path.join(cycles_dir, 'sample_drop_idxs.pkl')):
                f = open(
                    os.path.join(cycles_dir, 'sample_drop_idxs.pkl'), 'rb')
                sample_drop_idxs = pickle.load(f)
                samples_to_threshold = (
                    len(data['Sample'].unique())
                    - len(sample_drop_idxs.keys())
                    )
                print(f'Samples to threshold: {samples_to_threshold}')

            else:
                # initialize a dictionary to append indices to drop
                samples_to_threshold = len(data['Sample'].unique())
                print(f'Samples to threshold: {samples_to_threshold }')

                sample_drop_idxs = {}

            # drop samples previously run
            ratios_melt = ratios_melt[
                ~ratios_melt['sample'].isin(
                    [i for i in sample_drop_idxs.keys()])]

            if samples_to_threshold > 0:

                subprocess.call(
                    ['open', '-a', 'Preview', os.path.join(
                        cycles_dir, 'cycle_correlation(logRatio).pdf')])

                for name, group in ratios_melt.groupby(['sample']):

                    # loop over all cycles
                    # for cycle_ratio in group['cycle'].unique():
                    #     if not (
                    #       cycle_ratio.split('/')[0]
                    #       == cycle_ratio.split('/')[1]):

                    # loop over last cycle only
                    cycle_ratio = group['cycle'].unique()[-1]

                    # isolate ratio data
                    cycle_data = group[group['cycle'] == cycle_ratio]

                    cycle_num = cycle_ratio.split('/')[1]

                    # plot log(cycle 1/n) histogram for current sample
                    num_bins = 300
                    histtype = 'stepfilled'
                    sns.set_style('whitegrid')

                    # save dataframe of current cycle
                    cycle_data.to_parquet(
                        os.path.join(cycles_dir, 'cycle_data.parquet')
                        )

                    fig, ax = plt.subplots()
                    matplotlib_warnings(fig)

                    plt.subplots_adjust(left=0.25, bottom=0.25)
                    counts, bins, patches = plt.hist(
                        cycle_data['log10(ratio)'], bins=num_bins,
                        density=False, color='grey', ec='none',
                        alpha=0.75, histtype=histtype,
                        range=None, label='before'
                        )

                    log10 = '$log_{10}$'
                    plt.title(
                        f'Sample={name}   {log10}({dna1}/{cycle_num})',
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

                        # remove current lines
                        [i.remove() for i in ax.get_lines()]

                        # new cutoffs
                        lowerCutoff = sLower.val
                        upperCutoff = sUpper.val

                        # update plot with cutoffs
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

                        if text == self.tbEntry:
                            # reset tbEntry
                            self.tbEntry = None
                        else:
                            # assign current textbox entry to tbEntry
                            # if text and tbEntry are different
                            self.tbEntry = text

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
                                channel_number = marker_channel_number(
                                    markers, cycle_num)

                                # read first DNA channel
                                for file_path in glob.glob(
                                  f'{self.inDir}/tif/{text}*.tif'):
                                    dna_first = single_channel_pyramid(
                                        file_path, channel=0)

                                # read last DNA channel
                                for file_path in glob.glob(
                                  f'{self.inDir}/tif/{text}*.tif'):
                                    dna_last = single_channel_pyramid(
                                        file_path,
                                        channel=channel_number.item() - 1)

                                # read segmentation outlines
                                file_path = f'{self.inDir}/seg/{text}*.ome.tif'
                                seg = single_channel_pyramid(
                                    glob.glob(file_path)[0], channel=0
                                    )

                                # filter group data by selecting
                                # indices NOT in idxs
                                sample_data = data[data['Sample'] == name]
                                drop_df = sample_data.index.isin(idxs)
                                sample_centroids = sample_data[
                                    ['Y_centroid', 'X_centroid']][~drop_df]

                                napari_warnings()
                                viewer = napari.view_image(
                                    dna_last, rgb=False,
                                    blending='additive',
                                    colormap='magenta',
                                    name=f'{cycle_num}'
                                    )

                                viewer.add_image(
                                    dna_first, rgb=False, blending='additive',
                                    colormap='green',
                                    name=f'{dna1}'
                                    )

                                viewer.add_image(
                                    seg, rgb=False, blending='additive',
                                    opacity=0.5, colormap='red', visible=False,
                                    name='segmentation'
                                    )

                                viewer.add_points(
                                    sample_centroids,
                                    name='Selected Cells',
                                    properties=None,
                                    face_color='yellow',
                                    edge_color='k',
                                    edge_width=0.0, size=4.0
                                    )

                                napari.run()

                            else:
                                print('Must enter name of current sample.')

                    axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
                    text_box = TextBox(
                        axbox, 'evaluation sample name', initial='',
                        color='0.95',
                        hovercolor='1.0',
                        label_pad=0.05
                        )
                    text_box.on_submit(submit)
                    plt.show(block=True)

                    # ensure tbEntry is set to default (None) for future use
                    self.tbEntry = None

                    # get final lower and upper cutoffs
                    lowerCutoff, upperCutoff = update(val=None)

                    # read dataframe of current cycle and apply final cutoffs
                    cycle_data = pd.read_parquet(
                        os.path.join(cycles_dir, 'cycle_data.parquet')
                        )

                    # take all data if sliders not moved
                    if lowerCutoff == upperCutoff:
                        lowerCutoff = cycle_data['log10(ratio)'].min()
                        upperCutoff = cycle_data['log10(ratio)'].max()

                    # otherwise, get indices outside lower and upper cutoffs
                    # and append to sample_drop_idxs dictionary
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
            print()

            # filter data with indices to frop from all samples
            print('X-axis gating: Dropping unstable cells from all samples...')
            indices_to_drop = set()
            for k, v in sample_drop_idxs.items():
                for idx in v:
                    indices_to_drop.update(set([idx]))

            df = data.loc[~data.index.isin(indices_to_drop)]
            plt.close('all')

        # grab dna and sample columns of filtered dataframe
        facet_input = df.loc[
            :, df.columns.str.contains(f'{dna_moniker}|Sample')].copy()

        # melt filtered dataframe
        facet_per_cycle_melt = (
            facet_input
            .sample(frac=1.0)
            .reset_index()
            .melt(id_vars=['Sample', 'index'], var_name='cycle')
            )

        # convert sample and cycle columns to ordered categorical
        # with sorted catgories by natsorted
        facet_per_cycle_melt['Sample'] = pd.Categorical(
            facet_per_cycle_melt['Sample'], ordered=True,
            categories=natsorted(
                facet_per_cycle_melt['Sample'].unique()))
        facet_per_cycle_melt['cycle'] = pd.Categorical(
            facet_per_cycle_melt['cycle'], ordered=True,
            categories=natsorted(
                facet_per_cycle_melt['cycle'].unique()))

        # sort melt on sample, cycle, and index
        facet_per_cycle_melt = facet_per_cycle_melt.sort_values(
            ['Sample', 'cycle', 'index']
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
                    == dna1], y,
                s=0.05, alpha=0.1, linewidth=None,
                marker='o', c='r'), 'value')

        plt.savefig(
            os.path.join(
                cycles_dir, 'cycle_correlation(perCycle).png'), dpi=600
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

        sample_color_dict = dict(zip(
            natsorted(facet_per_cycle_melt['Sample'].unique()),
            cmap.colors))

        g = sns.FacetGrid(
            facet_per_cycle_melt, col='cycle', hue='Sample',
            col_wrap=5, sharex=True, sharey=True
            )

        g.map(
            lambda sam, y, color, **kwargs: plt.scatter(
                facet_per_cycle_melt.loc[
                    (facet_per_cycle_melt['Sample'] ==
                     sam.unique()[0])
                    & (facet_per_cycle_melt['cycle'] == dna1),
                    'value'], y,
                c=np.reshape(sample_color_dict[sam.unique()[0]], (-1, 3)),
                s=0.05, linewidth=None, marker='o', **kwargs),
            'Sample', 'value'
            )

        plt.legend(markerscale=10, bbox_to_anchor=(1.1, 1.05))

        plt.savefig(
            os.path.join(
                cycles_dir, 'cycle_correlation(perSample).png'), dpi=600,
            bbox_inches='tight'
            )
        plt.close('all')
        print()

        # remove last sample groupby dataframe
        if os.path.exists(os.path.join(cycles_dir, 'cycle_data.parquet')):
            os.remove(os.path.join(cycles_dir, 'cycle_data.parquet'))

        print()
        print()
        return df

    @module
    def logTransform(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude)

        # log10 transform immunomarker signals
        data.loc[:, abx_channels] += 0.00000000001
        data.loc[:, abx_channels] = np.log10(data[abx_channels].copy())

        print()
        print()
        return data

    @module
    def pruneOutliers(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            )

        # create pruning directory if it doesn't already exist
        pruning_dir = os.path.join(self.outDir, 'pruning')
        if not os.path.exists(pruning_dir):
            os.mkdir(pruning_dir)

        # if percentile cutoffs have already been assigned for some antibodies
        # open pruning dict, get list of remaining markers to prune
        # sorted by markers.csv, and read partially-pruned data
        if os.path.exists(os.path.join(pruning_dir, 'pruning_dict.pkl')):

            f = open(os.path.join(pruning_dir, 'pruning_dict.pkl'), 'rb')
            pruning_dict = pickle.load(f)

            total_markers = set(abx_channels)
            completed_markers = set(pruning_dict.keys())
            remaining_markers = total_markers.difference(completed_markers)

            marker_idxs = {}
            for marker in remaining_markers:
                channel_number = marker_channel_number(markers, marker)
                marker_idxs[marker] = channel_number.item()
            markers_to_prune = [
                k for k, v in
                sorted(marker_idxs.items(), key=lambda item: item[1])]

            if os.path.exists(os.path.join(pruning_dir, 'data_copy1.parquet')):
                data_copy1 = pd.read_parquet(
                    os.path.join(pruning_dir, 'data_copy1.parquet'))

            print(f'Immunomarker channels to prune: {len(markers_to_prune)}')
            print()

        # if percentile cutoffs have not yet been assigned for any antibodies,
        # assign all markers for pruning, create pruning_dict, and make a
        # copy of input data to be serially redcated
        else:
            pruning_dict = {}
            total_markers = set(abx_channels)
            marker_idxs = {}
            for marker in total_markers:
                channel_number = marker_channel_number(markers, marker)
                marker_idxs[marker] = channel_number.item()
            markers_to_prune = [
                k for k, v in
                sorted(marker_idxs.items(), key=lambda item: item[1])]

            data_copy1 = data.copy()

            print(f'Immunomarker channels to prune: {len(markers_to_prune)}')
            print()

        # view raw signal intensities for all samples
        for ab in markers_to_prune:

            print(ab)

            hist_facet = (
                data_copy1[['Sample', 'Condition', 'Area'] + [ab]]
                .sample(frac=1.0)
                .melt(id_vars=['Sample', 'Condition', 'Area'],
                      var_name='channel', value_name='signal'))

            # naturally sort hist_facet by sample
            hist_facet['Sample'] = pd.Categorical(
                hist_facet['Sample'], ordered=True,
                categories=natsorted(hist_facet['Sample'].unique()))
            hist_facet.sort_values('Sample', inplace=True)
            hist_facet['Sample'] = hist_facet['Sample'].astype('str')

            # create column for facet labels
            hist_facet['for_plot'] = (
                hist_facet['Sample'] + ', ' +
                hist_facet['Condition']
                )

            # plot facets
            sns.set_style('white')
            g = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=15,
                height=3.0, aspect=1.0, sharex=True, sharey=False,
                )

            # use hexbins for plotting
            if self.hexbins:
                g.map(
                    lambda x, y, color: plt.hexbin(
                        x, y, gridsize=self.hexbinGridSize,
                        linewidths=0.02, color='dimgrey'),
                    'signal', 'Area')

            # use scatter plots for plotting
            else:
                g.map(
                    lambda x, y, color: plt.scatter(
                        x, y, s=1.0, linewidths=0.0, color='k'),
                    'signal', 'Area')

            g.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=9.0, pad=2.0)

            g.fig.suptitle(ab, y=1.1, size=20.0)

            for ax in g.axes.flatten():
                ax.tick_params(
                    axis='both', which='major',
                    labelsize=5.0, pad=-2)

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
                top=0.90, hspace=0.4, wspace=0.4)

            plt.savefig(
                os.path.join(
                    pruning_dir,
                    f'{ab}_raw.png'), dpi=300,
                bbox_inches='tight')
            plt.close()

            # show raw signal intensity distributions
            subprocess.call(
                ['open', '-a', 'Preview', os.path.join(
                    pruning_dir,
                    f'{ab}_raw.png')]
                    )

            def submit(text):

                lowerPercentileCutoff = float(text.split(', ')[0])
                upperPercentileCutoff = float(text.split(', ')[1])

                # check whether a sample name was entered for cutoff fitting
                num_inputs = len(text.split(', '))
                if num_inputs == 3:
                    sample_name = str(text.split(', ')[2])
                else:
                    sample_name = None

                # add entered cutoffs to pruning_dict
                pruning_dict[ab] = (
                    lowerPercentileCutoff,
                    upperPercentileCutoff)

                # create a copy of partially-pruned data for testing cutoffs
                data_copy2 = data_copy1[[
                    'Sample', 'Condition', 'Area'] + [ab]].copy()

                # create running lists of indices removed according
                # to lower and upper cutoffs (used for data viz in Napari)
                total_low_idxs = []
                total_high_idxs = []

                # apply current percentile cutoffs to individual samples
                for sample in natsorted(data_copy2['Sample'].unique()):

                    # drop cells < lower cutoff and > than upper cutoff
                    indices_to_drop = []

                    sample_channel_data = data_copy2[
                        data_copy2['Sample'] == sample][ab]

                    low_drop_idxs = sample_channel_data.index[
                        sample_channel_data < np.percentile(
                            sample_channel_data, lowerPercentileCutoff)]
                    indices_to_drop.extend(low_drop_idxs)

                    high_drop_idxs = sample_channel_data.index[
                        sample_channel_data > np.percentile(
                            sample_channel_data, upperPercentileCutoff)]
                    indices_to_drop.extend(high_drop_idxs)

                    data_copy2.drop(
                        labels=set(indices_to_drop), axis=0,
                        inplace=True, errors='raise')

                    # rescale residual signal intensities
                    pruned_data = data_copy2[
                        data_copy2['Sample'] == sample][ab]

                    scaler = (
                        MinMaxScaler(feature_range=(0, 1), copy=True)
                        .fit(pruned_data.values.reshape(-1, 1)))
                    rescaled_data = scaler.transform(
                        pruned_data.values.reshape(-1, 1))
                    rescaled_data = pd.DataFrame(
                        data=rescaled_data,
                        index=pruned_data.index,
                        ).rename(columns={0: ab})

                    # update residual signal intensities
                    data_copy2.update(rescaled_data)

                    # update running lists of total indices
                    total_low_idxs.extend(low_drop_idxs)
                    total_high_idxs.extend(high_drop_idxs)

                # melt pruned and rescaled data_copy2
                hist_facet = (
                    data_copy2
                    .sample(frac=1.0)
                    .melt(id_vars=['Sample', 'Condition', 'Area'],
                          var_name='channel', value_name='signal'))

                # naturally sort hist_facet by sample
                hist_facet['Sample'] = pd.Categorical(
                    hist_facet['Sample'], ordered=True,
                    categories=natsorted(hist_facet['Sample'].unique()))
                hist_facet.sort_values('Sample', inplace=True)
                hist_facet['Sample'] = hist_facet['Sample'].astype('str')

                # create column for facet labels
                hist_facet['for_plot'] = (
                    hist_facet['Sample'] + ', ' +
                    hist_facet['Condition'])

                # plot pruned facets
                g = sns.FacetGrid(
                    hist_facet, col='for_plot', col_wrap=15,
                    height=3.0, aspect=1.0, sharex=True, sharey=False
                    )

                # use hexbins for plotting
                if self.hexbins:
                    g.map(
                        lambda x, y, color: plt.hexbin(
                            x, y, gridsize=self.hexbinGridSize,
                            linewidths=0.02, color='dimgrey'),
                        'signal', 'Area')

                # use scatter plots for plotting
                else:
                    g.map(
                        lambda x, y, color: plt.scatter(
                            x, y, s=1.0, linewidths=0.0, color='k'),
                        'signal', 'Area')

                g.set_titles(
                    col_template="{col_name}", fontweight='bold',
                    size=9.0, pad=2.0
                    )

                g.fig.suptitle(ab, y=1.1, size=20.0)

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
                        f'{ab}_pruned_rescaled.png'), dpi=300,
                    bbox_inches='tight'
                        )
                plt.close()

                # show pruned and rescaled signal intensity distributions
                subprocess.call(
                    ['open', '-a', 'Preview', os.path.join(
                        pruning_dir,
                        f'{ab}_pruned_rescaled.png')]
                        )

                # if sample_name was entered,
                # plot low and high outliers in Napri
                if sample_name:

                    if text == self.tbEntry:
                        # reset tbEntry
                        self.tbEntry = None
                    else:
                        # assign current textbox entry to tbEntry
                        # if text and tbEntry are different
                        self.tbEntry = text

                        napari_warnings()

                        channel_number = marker_channel_number(markers, ab)

                        # read DNA1 channel
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{sample_name}*.tif'):
                            dna = single_channel_pyramid(
                                file_path, channel=0)

                        # read antibody channel
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{sample_name}*.tif'):
                            channel = single_channel_pyramid(
                                file_path, channel=channel_number.item() - 1)

                        # read segmentation outlines
                        file_path = f'{self.inDir}/seg/{sample_name}*.ome.tif'
                        seg = single_channel_pyramid(
                            glob.glob(file_path)[0], channel=0
                            )

                        # grab centroids of low signal intensity outliers
                        low_centroids = data_copy1[
                            ['Y_centroid', 'X_centroid']][
                            (data_copy1.index.isin(total_low_idxs)) &
                            (data_copy1['Sample'] == sample_name)]

                        # grab centroids of high signal intensity outliers
                        high_centroids = data_copy1[
                            ['Y_centroid', 'X_centroid']][
                            (data_copy1.index.isin(total_high_idxs)) &
                            (data_copy1['Sample'] == sample_name)]

                        # decorate Napari viewer
                        viewer = napari.view_image(
                            dna, rgb=False, opacity=0.25,
                            name=dna1
                            )

                        viewer.add_image(
                            channel, rgb=False, blending='additive',
                            colormap='green', visible=False,
                            name=ab)

                        viewer.add_image(
                            seg, rgb=False, blending='additive',
                            opacity=0.5, colormap='red', visible=False,
                            name='segmentation')

                        viewer.add_points(
                            low_centroids,
                            name='low centroids',
                            properties=None,
                            face_color='magenta',
                            edge_color='k',
                            edge_width=0.0, size=8.0)

                        viewer.add_points(
                            high_centroids,
                            name='high centroids',
                            properties=None,
                            face_color='cyan',
                            edge_color='k',
                            edge_width=0.0, size=8.0)

                        napari.run()

            fig = plt.figure(figsize=(6, 2))
            matplotlib_warnings(fig)
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

            # ensure tbEntry is set to default (None) for future use
            self.tbEntry = None

            # take all data if cutoff window is closed without entering values
            if ab not in pruning_dict.keys():
                pruning_dict[ab] = (0.0, 100.0)

            # update data_copy1 (pruned dataframe) with selected cutoffs
            for sample in natsorted(data_copy1['Sample'].unique()):

                sample_channel_data = data_copy1[
                    data_copy1['Sample'] == sample][ab]

                # drop cells < lower cutoff and > than upper cutoff
                indices_to_drop = []

                indices_to_drop.extend(
                    sample_channel_data.index[
                        sample_channel_data < np.percentile(
                            sample_channel_data,
                            pruning_dict[ab][0])])

                indices_to_drop.extend(
                    sample_channel_data.index[
                        sample_channel_data > np.percentile(
                            sample_channel_data,
                            pruning_dict[ab][1])])

                data_copy1.drop(
                    labels=set(indices_to_drop), axis=0,
                    inplace=True, errors='raise')

                # rescale pruned antibody signal intensities
                pruned_data = data_copy1[data_copy1['Sample'] == sample][ab]

                scaler = (
                    MinMaxScaler(feature_range=(0, 1), copy=True)
                    .fit(pruned_data.values.reshape(-1, 1)))
                rescaled_data = scaler.transform(
                    pruned_data.values.reshape(-1, 1))
                rescaled_data = pd.DataFrame(
                    data=rescaled_data,
                    index=pruned_data.index,
                    ).rename(columns={0: ab})

                data_copy1.update(rescaled_data)

            # melt pruned and rescaled data_copy1
            hist_facet = (
                data_copy1[['Sample', 'Condition', 'Area'] + [ab]]
                .sample(frac=1.0)
                .melt(id_vars=['Sample', 'Condition', 'Area'],
                      var_name='channel', value_name='signal'))

            # naturally sort hist_facet by sample
            hist_facet['Sample'] = pd.Categorical(
                hist_facet['Sample'], ordered=True,
                categories=natsorted(hist_facet['Sample'].unique()))
            hist_facet.sort_values('Sample', inplace=True)
            hist_facet['Sample'] = hist_facet['Sample'].astype('str')

            # create column for facet labels
            hist_facet['for_plot'] = (
                hist_facet['Sample'] + ', ' +
                hist_facet['Condition']
                )

            # re-plot pruned facets
            g = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=15,
                height=3.0, aspect=1.0, sharex=True, sharey=False
                )

            # use hexbins for plotting
            if self.hexbins:
                g.map(
                    lambda x, y, color: plt.hexbin(
                        x, y, gridsize=self.hexbinGridSize,
                        linewidths=0.02, color='dimgrey'),
                    'signal', 'Area')

            # use scatter plots for plotting
            else:
                g.map(
                    lambda x, y, color: plt.scatter(
                        x, y, s=1.0, linewidths=0.0, color='k'),
                    'signal', 'Area')

            g.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=9.0, pad=2.0
                )

            g.fig.suptitle(ab, y=1.1, size=20.0)

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
                    f'{ab}_pruned_rescaled.png'), dpi=300,
                bbox_inches='tight'
                    )
            plt.close()

            # save updated pruning_dict
            f = open(os.path.join(pruning_dir, 'pruning_dict.pkl'), 'wb')
            pickle.dump(pruning_dict, f)
            f.close()

            # save updated (pruned) data_copy1
            data_copy1.to_parquet(
                os.path.join(pruning_dir, 'data_copy1.parquet'))
        print()

        # apply antibody percentile cutoffs to actual dataframe (i.e. data)
        # in the same order they were originally curated
        # (according to markers.csv order). This assumes pruning_dict
        # preserves the order in which antibodies were appended.
        # As of Python 3.6, for the CPython implementation of Python,
        # dictionaries remember the order of items inserted.
        for k, v in pruning_dict.items():
            print(f'Applying percentile cutoffs to the {k} channel.')

            for sample in natsorted(data_copy1['Sample'].unique()):

                sample_channel_data = data[data['Sample'] == sample][k]

                # drop cells < lower cutoff and > than upper cutoff
                indices_to_drop = []

                indices_to_drop.extend(
                    sample_channel_data.index[
                        sample_channel_data < np.percentile(
                            sample_channel_data, v[0])])

                indices_to_drop.extend(
                    sample_channel_data.index[
                        sample_channel_data > np.percentile(
                            sample_channel_data, v[1])])

                data.drop(
                    labels=set(indices_to_drop), axis=0,
                    inplace=True, errors='raise')

                # rescale pruned antibody signal intensities
                pruned_data = data[data['Sample'] == sample][k]

                scaler = (
                    MinMaxScaler(feature_range=(0, 1), copy=True)
                    .fit(pruned_data.values.reshape(-1, 1)))
                rescaled_data = scaler.transform(
                    pruned_data.values.reshape(-1, 1))
                rescaled_data = pd.DataFrame(
                    data=rescaled_data,
                    index=pruned_data.index,
                    ).rename(columns={0: k})

                data.update(rescaled_data)

        print()
        print()
        return data

    @module
    def metaQC(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude)

        # drop antibody channel exclusions for metaQC clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClusteringQC]

        # create metaQC directory if metaQC is to be performed and
        # the directory doesn't already exist
        reclass_dir = os.path.join(
            self.outDir, f'clustering/metaQC')
        if not os.path.exists(reclass_dir):
            os.makedirs(reclass_dir)

        # specify the names of modules in that pipeline that perform
        # data redaction leading up to the metQC modules
        modules = ['aggregateData', 'selectROIs', 'intensityFilter',
                   'areaFilter', 'cycleCorrelation', 'pruneOutliers']

        # build a dictionary of redacted data returned by each module (clean)
        module_dict = {}
        for module_idx, module in enumerate(modules):
            data = pd.read_parquet(
                os.path.join(
                    self.outDir, f'checkpoints/{module}.parquet'))
            module_dict[module_idx] = [module, data]

        #######################################################################
        # build QCData: QCData is a combination of clean
        # and noisy data

        if self.metaQC:
            # loop over modules to get data redacted by each module (noisy)
            for module_idx in [i for i in module_dict.keys()][:-1]:

                # isolate data
                noise = module_dict[module_idx][1][
                    ~module_dict[module_idx][1].index.isin(
                        module_dict[module_idx+1][1].index)].copy()

                # if noisy data exists, add a QC stage column
                if not noise.empty:
                    noise.loc[:, 'Filter'] = (
                        module_dict[module_idx+1][0])

                # append to module_dict
                module_dict[module_idx+1].append(noise)

            # combine noisy data from all modules into a single dataframe
            if self.delintMode:
                # if negative ROI selection, include gated cells as noisy data
                noisyData = pd.DataFrame()
                for module_idx, module in enumerate(modules):
                    # the number of dictionary value enteries should be 3:
                    # module name, clean data, noisy data if noisy data exists
                    if len(module_dict[module_idx]) == 3:
                        noisyData = noisyData.append(
                            module_dict[module_idx][2])
            else:
                # if postive ROI selection, exlcude cells not gated on
                noisyData = pd.DataFrame()
                for module_idx, module in enumerate(modules):
                    if len(module_dict[module_idx]) == 3:
                        if not module == 'selectROIs':
                            noisyData = noisyData.append(
                                module_dict[module_idx][2])

            if noisyData.empty:
                print('No data were filtered during prior QC steps. ' +
                      'Returning unfiltered data without ' +
                      'reclassification.')
                return data

            # create QC_status column for combined noisy data
            noisyData.loc[:, 'QC_status'] = 'noisy'

            # get raw version (untransformed, not rescaled) of data
            # from fully-redacted dataframe
            first_module_idx = modules.index(modules[0])
            last_module_idx = modules.index(modules[-1])

            cleanDataRescaled = module_dict[last_module_idx][1]

            cleanDataRaw = module_dict[first_module_idx][1][
                module_dict[first_module_idx][1].index.isin(
                    cleanDataRescaled.index)].copy()

            # create QC_status column for selected clean data
            cleanDataRaw.loc[:, 'QC_status'] = 'clean'

            # specify channels on which to perform metaQC clustering
            abx_channels_dna = [dna1] + abx_channels

            # if QCData.pkl exists, read the pickle
            if os.path.exists(os.path.join(reclass_dir, 'QCData.pkl')):
                f = open(os.path.join(
                    reclass_dir, 'QCData.pkl'), 'rb')
                QCData = pickle.load(f)

                # read current chunk index
                f = open(os.path.join(reclass_dir, 'chunk_index.pkl'), 'rb')
                chunk_index = pickle.load(f)

                # read current chunk if one exists
                if os.path.exists(os.path.join(reclass_dir, 'chunk.pkl')):
                    f = open(os.path.join(
                        reclass_dir, 'chunk.pkl'), 'rb')
                    chunk = pickle.load(f)
            else:
                # if QCData.pkl doesn't exist, append noisyData
                # to cleanDataRaw, row-wise
                QCData = cleanDataRaw.append(noisyData)

                # shuffle QCData row order to randomize cells
                # for batch clustering
                QCData = QCData.sample(frac=1.0, random_state=5)
                QCData.reset_index(drop=True, inplace=True)

                # log-transform QCData
                QCData.loc[:, abx_channels_dna] += 0.00000000001
                QCData.loc[:, abx_channels_dna] = np.log10(
                    QCData[abx_channels_dna])

                # rescale channel signal intensities (0-1)
                for marker in abx_channels_dna:
                    channel_data = QCData[marker]

                    scaler = (
                        MinMaxScaler(feature_range=(0, 1), copy=True)
                        .fit(channel_data.values.reshape(-1, 1))
                            )
                    rescaled_data = scaler.transform(
                        channel_data.values.reshape(-1, 1)
                        )

                    rescaled_data = pd.DataFrame(
                        data=rescaled_data,
                        index=channel_data.index,
                        ).rename(columns={0: marker})

                    QCData.update(rescaled_data)

                # initialize chunk index counter at 0
                chunk_index = 0
                f = open(os.path.join(reclass_dir, 'chunk_index.pkl'), 'wb')
                pickle.dump(chunk_index, f)
                f.close()

                # save QCData.pkl
                f = open(os.path.join(reclass_dir, 'QCData.pkl'), 'wb')
                pickle.dump(QCData, f)
                f.close()

            ###################################################################
            # chunk QCData

            # specify the number of cells in each clustering batch
            # (this limits clustering times and memory pressure)
            batch_size = 500000

            # ensure minimun batch size equals the batch_size variable
            # otherwise cluster all data at once
            if len(QCData) < (batch_size)*2:
                num_chunks = 1
                chunks = np.array_split(QCData, num_chunks)
            else:
                num_chunks = math.ceil(len(QCData)/batch_size)
                chunks = np.array_split(QCData, num_chunks)

            # replace the specific chunk in the chunks list
            # associated with the current chunk_index with
            # the chunk.pkl which has clustering labels.
            # This allows the program to avoid re-clustering the first chunk
            try:
                chunks[chunk_index] = chunk

            except UnboundLocalError:
                pass

            ###################################################################
            # initialize reclassification storage

            # if reclassified data storage dict exists, open it
            if os.path.exists(
              os.path.join(reclass_dir, 'reclass_storage_dict.pkl')):
                f = open(os.path.join(
                    reclass_dir, 'reclass_storage_dict.pkl'), 'rb')
                reclass_storage_dict = pickle.load(f)
            # else, initialize reclassified data storage dict for
            # reclassified clean and noisy data
            else:
                reclass_storage_dict = {}
                for reclass in ['clean', 'noisy']:
                    reclass_storage_dict[reclass] = pd.DataFrame()
                f = open(os.path.join(
                    reclass_dir, 'reclass_storage_dict.pkl'), 'wb')
                pickle.dump(reclass_storage_dict, f)
                f.close()
            print()

            # assign default HDBSCAN minimum cluster size
            default_range_tuple = 250

            # assign default reclassification cutoffs
            default_reclass_tuple = (0.75, 0.75)

            # assign default sample for viz
            default_sample = self.viewSample

            # define the degree of transparency for cells in the embedding
            alpha = 1.0

            ###################################################################
            # define nested plotting functions

            def do_plots():

                # initialize counter to allow for different
                # behavior on the first round of plotting
                iteration = 1

                sns.set_style('whitegrid')

                fig = plt.figure(figsize=(17, 6))
                matplotlib_warnings(fig)

                gs = plt.GridSpec(2, 4, figure=fig, height_ratios=[1, 10])

                # define text box axes
                ax_mcs_tb = fig.add_subplot(gs[0, 0])
                ax_reclass_tb = fig.add_subplot(gs[0, 1])
                ax_sample_tb = fig.add_subplot(gs[0, 3])

                # define embedding axes
                ax_sc_cluster = fig.add_subplot(gs[1, 0])
                ax_sc_status = fig.add_subplot(gs[1, 1])
                ax_sc_reclass = fig.add_subplot(gs[1, 2])
                ax_sc_sample = fig.add_subplot(gs[1, 3])

                # assign default inputs
                range_tuple = default_range_tuple
                reclass_tuple = default_reclass_tuple
                sample = default_sample
                selector = None

                def cluster_and_plot():
                    nonlocal iteration
                    nonlocal selector
                    nonlocal range_tuple
                    nonlocal reclass_tuple

                    if not isinstance(range_tuple, tuple):
                        range_tuple = (range_tuple, )

                    # loop for single-value MCS inputs
                    if iteration == 1:
                        if len(range_tuple) == 1:

                            mylist = [range_tuple[0]]

                            for i in mylist:

                                min_cluster_size = i

                                clustering = hdbscan.HDBSCAN(
                                    min_cluster_size=min_cluster_size,
                                    min_samples=None,
                                    metric='euclidean', alpha=1.0, p=None,
                                    algorithm='best', leaf_size=40,
                                    memory=Memory(
                                        location=None),
                                    approx_min_span_tree=True,
                                    gen_min_span_tree=False,
                                    core_dist_n_jobs=-1,
                                    cluster_selection_method='eom',
                                    allow_single_cluster=False,
                                    prediction_data=False,
                                    match_reference_implementation=False).fit(
                                        chunk[['emb1', 'emb2']]
                                        )
                                chunk['cluster'] = clustering.labels_

                                print(
                                    f'min_cluster_size = {i}',
                                    np.unique(clustering.labels_)
                                    )
                                print()

                                # PLOT embedding
                                for color_by in ['cluster', 'QC_status',
                                                 'Reclass', 'Sample']:

                                    highlight = 'none'

                                    if color_by == 'cluster':

                                        # build cmap
                                        cmap = categorical_cmap(
                                            numUniqueSamples=len(
                                                chunk[color_by].unique()),
                                            numCatagories=10,
                                            cmap='tab10',
                                            continuous=False
                                            )

                                        # make black the first color to specify
                                        # unclustered cells (cluster -1)
                                        cmap = ListedColormap(
                                            np.insert(
                                                arr=cmap.colors, obj=0,
                                                values=[0, 0, 0], axis=0)
                                                )

                                        # trim cmap to # unique samples
                                        trim = (
                                            len(cmap.colors) - len(
                                                chunk[color_by].unique())
                                            )
                                        cmap = ListedColormap(
                                            cmap.colors[:-trim]
                                            )

                                        sample_dict = dict(
                                            zip(
                                                natsorted(
                                                    chunk[color_by].unique()),
                                                list(range(len(chunk[color_by]
                                                     .unique()))))
                                                )

                                        c = [sample_dict[i] for i
                                             in chunk[color_by]]

                                        ax_sc_cluster.cla()

                                        cluster_paths = ax_sc_cluster.scatter(
                                            chunk['emb1'],
                                            chunk['emb2'],
                                            c=c,
                                            alpha=alpha,
                                            s=point_size,
                                            cmap=cmap,
                                            ec=[
                                                'k' if i == highlight
                                                else 'none' for i in
                                                chunk[color_by]
                                                ],
                                            linewidth=0.1
                                            )

                                        ax_sc_cluster.set_title('HDBSCAN')
                                        ax_sc_cluster.axis('equal')
                                        ax_sc_cluster.tick_params(labelsize=5)
                                        ax_sc_cluster.grid(False)

                                        legend_elements = []
                                        for e, i in enumerate(
                                            natsorted(chunk[color_by].unique())
                                          ):

                                            hi_markers = cluster_expression(
                                                df=chunk, markers=abx_channels,
                                                cluster=i, num_proteins=3,
                                                standardize='within',
                                                )

                                            legend_elements.append(
                                                Line2D([0], [0], marker='o',
                                                       color='none',
                                                       label=(
                                                        f'Cluster {i}: '
                                                        f'{hi_markers}'),
                                                       markerfacecolor=(
                                                        cmap.colors[e]),
                                                       markeredgecolor='none',
                                                       lw=0.001, markersize=6))

                                        # cluster_lgd = ax_sc_cluster.legend(
                                        #     handles=legend_elements,
                                        #     prop={'size': 8},
                                        #     bbox_to_anchor=[1.0, 1.0]
                                        #     )

                                    elif color_by == 'QC_status':

                                        # build cmap
                                        cmap = ListedColormap(
                                            np.array([[0.91, 0.29, 0.235],
                                                      [0.18, 0.16, 0.15]])
                                             )

                                        sample_dict = dict(
                                            zip(
                                                natsorted(
                                                    chunk[
                                                        'QC_status'].unique()),
                                                list(range(len(
                                                    chunk['QC_status']
                                                     .unique()))))
                                                )

                                        c = [sample_dict[i] for
                                             i in chunk['QC_status']]

                                        ax_sc_status.cla()

                                        ax_sc_status.scatter(
                                            chunk['emb1'],
                                            chunk['emb2'],
                                            c=c,
                                            cmap=cmap,
                                            alpha=alpha,
                                            s=point_size,
                                            ec=[
                                                'k' if i == highlight
                                                else 'none' for i in
                                                chunk['QC_status']
                                                ],
                                            linewidth=0.1
                                            )

                                        ax_sc_status.set_title('QC Status')
                                        ax_sc_status.axis('equal')
                                        ax_sc_status.tick_params(labelsize=5)
                                        ax_sc_status.grid(False)

                                        legend_elements = []
                                        for e, i in enumerate(
                                            natsorted(
                                                chunk['QC_status'].unique())
                                          ):

                                            if i == highlight:
                                                markeredgecolor = 'k'
                                            else:
                                                markeredgecolor = 'none'

                                            sample_to_map = (
                                                chunk['Sample'][
                                                    chunk['QC_status'] == i]
                                                .unique()[0]
                                                )

                                            legend_elements.append(
                                                Line2D([0], [0], marker='o',
                                                       color='none',
                                                       label=i,
                                                       markerfacecolor=(
                                                        cmap.colors[e]
                                                        ),
                                                       markeredgecolor=(
                                                        markeredgecolor),
                                                       lw=0.001,
                                                       markersize=6)
                                                       )

                                        qc_lgd = ax_sc_status.legend(
                                            handles=legend_elements,
                                            prop={'size': 8}, loc='best'
                                            # bbox_to_anchor=[1.0, 1.0]
                                            )

                                    elif color_by == 'Reclass':
                                        print(
                                            'Applying reclassification ' +
                                            f'cutoffs: {reclass_tuple}')

                                        clean = pd.DataFrame()
                                        noisy = pd.DataFrame()
                                        for name, cluster in chunk.groupby(
                                          'cluster'):
                                            if name != -1:
                                                # if a cluster contains
                                                # >= n% noisy data,
                                                # reclassify all clustering
                                                # cells as noisy
                                                if (
                                                  (len(cluster[cluster[
                                                    'QC_status'] == 'clean'])
                                                   / len(cluster)) >=
                                                  reclass_tuple[0]):
                                                    clean = clean.append(
                                                        cluster)
                                                # elif a cluster contains
                                                # >= n% clean data,
                                                # reclassify all clustering
                                                # cells as clean
                                                elif (
                                                  (len(cluster[cluster[
                                                    'QC_status'] == 'noisy'])
                                                   / len(cluster)) >=
                                                  reclass_tuple[1]):
                                                    noisy = noisy.append(
                                                        cluster)
                                                # else keep respective
                                                # cell statuses
                                                else:
                                                    noisy = noisy.append(
                                                        cluster[cluster[
                                                            'QC_status'] ==
                                                            'noisy'])
                                                    clean = clean.append(
                                                        cluster[cluster[
                                                            'QC_status'] ==
                                                            'clean'])

                                        # consider -1 cells from clean data
                                        # to be "clean"
                                        clean_outliers = chunk[
                                            (chunk['cluster'] == -1) &
                                            (chunk['QC_status'] == 'clean')
                                            ].copy()
                                        clean = clean.append(clean_outliers)

                                        # consider -1 cells from noisy data
                                        # to be "noisy"
                                        noisy_outliers = chunk[
                                            (chunk['cluster'] == -1) &
                                            (chunk['QC_status'] == 'noisy')
                                            ].copy()
                                        noisy = noisy.append(noisy_outliers)

                                        chunk['Reclass'] = 'init'
                                        chunk.loc[
                                            chunk.index.isin(clean.index),
                                            'Reclass'] = 'clean'
                                        chunk.loc[
                                            chunk.index.isin(noisy.index),
                                            'Reclass'] = 'noisy'

                                        # build cmap
                                        cmap = ListedColormap(
                                            np.array([[0.91, 0.29, 0.235],
                                                      [0.18, 0.16, 0.15]]))

                                        sample_dict = dict(
                                            zip(
                                                natsorted(
                                                    chunk['Reclass'].unique()),
                                                list(range(len(chunk['Reclass']
                                                     .unique()))))
                                                )

                                        c = [sample_dict[i] for
                                             i in chunk['Reclass']]

                                        ax_sc_reclass.cla()

                                        ax_sc_reclass.scatter(
                                            chunk['emb1'],
                                            chunk['emb2'],
                                            c=c,
                                            cmap=cmap,
                                            alpha=alpha,
                                            s=point_size,
                                            ec=[
                                                'k' if i == highlight
                                                else 'none' for i in
                                                chunk['Reclass']
                                                ],
                                            linewidth=0.1
                                            )

                                        ax_sc_reclass.set_title(
                                            'Reclassification')
                                        ax_sc_reclass.axis('equal')
                                        ax_sc_reclass.tick_params(labelsize=5)
                                        ax_sc_reclass.grid(False)

                                        legend_elements = []
                                        for e, i in enumerate(
                                            natsorted(
                                                chunk['Reclass'].unique())
                                          ):

                                            if i == highlight:
                                                markeredgecolor = 'k'
                                            else:
                                                markeredgecolor = 'none'

                                            sample_to_map = (
                                                chunk['Sample'][
                                                    chunk['Reclass'] == i]
                                                .unique()[0]
                                                )

                                            legend_elements.append(
                                                Line2D([0], [0], marker='o',
                                                       color='none',
                                                       label=i,
                                                       markerfacecolor=(
                                                        cmap.colors[e]
                                                        ),
                                                       markeredgecolor=(
                                                        markeredgecolor),
                                                       lw=0.001,
                                                       markersize=6)
                                                       )

                                        reclass_lgd = ax_sc_reclass.legend(
                                            handles=legend_elements,
                                            prop={'size': 8}, loc='best'
                                            # bbox_to_anchor=[1.0, 1.0]
                                            )

                                    elif color_by == 'Sample':

                                        # build cmap
                                        cmap = categorical_cmap(
                                            numUniqueSamples=len(
                                                chunk['Sample'].unique()),
                                            numCatagories=10,
                                            cmap='tab10',
                                            continuous=False
                                            )

                                        sample_dict = dict(
                                            zip(
                                                natsorted(
                                                    chunk['Sample'].unique()),
                                                list(range(len(chunk['Sample']
                                                     .unique()))))
                                                )

                                        c = [sample_dict[i] for
                                             i in chunk['Sample']]

                                        ax_sc_sample.cla()

                                        ax_sc_sample.scatter(
                                            chunk['emb1'],
                                            chunk['emb2'],
                                            c=c,
                                            cmap=cmap,
                                            alpha=alpha,
                                            s=point_size,
                                            ec=[
                                                'k' if i == highlight
                                                else 'none' for i in
                                                chunk['Sample']
                                                ],
                                            linewidth=0.1
                                            )
                                        ax_sc_sample.set_title('Sample')

                                        ax_sc_sample.axis('equal')
                                        ax_sc_sample.tick_params(labelsize=5)
                                        ax_sc_sample.grid(False)

                                        legend_elements = []
                                        for e, i in enumerate(
                                            natsorted(chunk['Sample'].unique())
                                          ):

                                            if i == highlight:
                                                markeredgecolor = 'k'
                                            else:
                                                markeredgecolor = 'none'

                                            sample_to_map = (
                                                chunk['Sample'][
                                                    chunk['Sample'] == i]
                                                .unique()[0]
                                                )

                                            legend_elements.append(
                                                Line2D([0], [0], marker='o',
                                                       color='none',
                                                       label=i,
                                                       markerfacecolor=(
                                                        cmap.colors[e]
                                                        ),
                                                       markeredgecolor=(
                                                        markeredgecolor),
                                                       lw=0.001,
                                                       markersize=6)
                                                       )

                                        sample_lgd = ax_sc_sample.legend(
                                            handles=legend_elements,
                                            prop={'size': 8}, loc='best'
                                            # bbox_to_anchor=[1.0, 1.0]
                                            )

                                fig.tight_layout()

                                # Must call draw() before creating selector,
                                # or alpha setting doesn't work.
                                fig.canvas.draw()
                                print('Done!')

                                if selector:
                                    selector.disconnect()
                                selector = SelectFromCollection(
                                    ax_sc_cluster, cluster_paths
                                    )

                                # allow for different behavior
                                # on first iteration of plotting
                                iteration += 1

                                return min_cluster_size, clean, noisy

                        else:
                            # loop for a range of input values
                            mylist = list(
                                range(range_tuple[0], range_tuple[1] + 1, 1)
                                )
                            mylist.reverse()  # run higher cluster sizes first

                            for i in mylist:

                                min_cluster_size = i

                                clustering = hdbscan.HDBSCAN(
                                    min_cluster_size=min_cluster_size,
                                    min_samples=None,
                                    metric='euclidean', alpha=1.0, p=None,
                                    algorithm='best', leaf_size=40,
                                    memory=Memory(
                                        location=None),
                                    approx_min_span_tree=True,
                                    gen_min_span_tree=False,
                                    core_dist_n_jobs=-1,
                                    cluster_selection_method='eom',
                                    allow_single_cluster=False,
                                    prediction_data=False,
                                    match_reference_implementation=False).fit(
                                        chunk[['emb1', 'emb2']]
                                        )
                                chunk['cluster'] = clustering.labels_

                                print(
                                    f'min_cluster_size = {i}',
                                    np.unique(clustering.labels_)
                                    )
                            print()

                            # acts to terminate loop
                            return min_cluster_size

                    if iteration != 1:
                        if len(range_tuple) == 1:

                            mylist = [range_tuple[0]]

                            for i in mylist:

                                min_cluster_size = i

                                clustering = hdbscan.HDBSCAN(
                                    min_cluster_size=min_cluster_size,
                                    min_samples=None,
                                    metric='euclidean', alpha=1.0, p=None,
                                    algorithm='best', leaf_size=40,
                                    memory=Memory(
                                        location=None),
                                    approx_min_span_tree=True,
                                    gen_min_span_tree=False,
                                    core_dist_n_jobs=-1,
                                    cluster_selection_method='eom',
                                    allow_single_cluster=False,
                                    prediction_data=False,
                                    match_reference_implementation=False).fit(
                                        chunk[['emb1', 'emb2']]
                                        )
                                chunk['cluster'] = clustering.labels_

                                print(
                                    f'min_cluster_size = {i}',
                                    np.unique(clustering.labels_)
                                    )
                                print()

                                # PLOT embedding
                                for color_by in ['cluster']:

                                    highlight = 'none'

                                    if color_by == 'cluster':

                                        # build cmap
                                        cmap = categorical_cmap(
                                            numUniqueSamples=len(
                                                chunk[color_by].unique()),
                                            numCatagories=10,
                                            cmap='tab10',
                                            continuous=False
                                            )

                                        # make black the first color to specify
                                        # unclustered cells (cluster -1)
                                        cmap = ListedColormap(
                                            np.insert(
                                                arr=cmap.colors, obj=0,
                                                values=[0, 0, 0], axis=0)
                                                )

                                        # trim cmap to # unique samples
                                        trim = (
                                            len(cmap.colors) - len(
                                                chunk[color_by].unique())
                                            )
                                        cmap = ListedColormap(
                                            cmap.colors[:-trim]
                                            )

                                        sample_dict = dict(
                                            zip(
                                                natsorted(
                                                    chunk[color_by].unique()),
                                                list(range(len(chunk[color_by]
                                                     .unique()))))
                                                )

                                        c = [sample_dict[i] for i
                                             in chunk[color_by]]

                                        ax_sc_cluster.cla()

                                        cluster_paths = ax_sc_cluster.scatter(
                                            chunk['emb1'],
                                            chunk['emb2'],
                                            c=c,
                                            alpha=alpha,
                                            s=point_size,
                                            cmap=cmap,
                                            ec=[
                                                'k' if i == highlight
                                                else 'none' for i in
                                                chunk[color_by]
                                                ],
                                            linewidth=0.1
                                            )

                                        ax_sc_cluster.set_title('HDBSCAN')
                                        ax_sc_cluster.axis('equal')
                                        ax_sc_cluster.tick_params(labelsize=5)
                                        ax_sc_cluster.grid(False)

                                        legend_elements = []
                                        for e, i in enumerate(
                                            natsorted(chunk[color_by].unique())
                                          ):

                                            hi_markers = cluster_expression(
                                                df=chunk, markers=abx_channels,
                                                cluster=i, num_proteins=3,
                                                standardize='within',
                                                )

                                            legend_elements.append(
                                                Line2D([0], [0], marker='o',
                                                       color='none',
                                                       label=(
                                                        f'Cluster {i}: '
                                                        f'{hi_markers}'),
                                                       markerfacecolor=(
                                                        cmap.colors[e]),
                                                       markeredgecolor='none',
                                                       lw=0.001, markersize=6))

                                        # cluster_lgd = ax_sc_cluster.legend(
                                        #     handles=legend_elements,
                                        #     prop={'size': 8},
                                        #     bbox_to_anchor=[1.0, 1.0]
                                        #     )

                                fig.tight_layout()

                                # Must call draw() before creating selector,
                                # or alpha setting doesn't work.
                                fig.canvas.draw()

                                if selector:
                                    selector.disconnect()
                                selector = SelectFromCollection(
                                    ax_sc_cluster, cluster_paths
                                    )

                                # acts to terminate loop
                                return min_cluster_size

                        else:
                            # loop for a range of input values
                            mylist = list(
                                range(range_tuple[0], range_tuple[1] + 1, 1)
                                )
                            mylist.reverse()  # run higher cluster sizes first

                            for i in mylist:

                                min_cluster_size = i

                                clustering = hdbscan.HDBSCAN(
                                    min_cluster_size=min_cluster_size,
                                    min_samples=None,
                                    metric='euclidean', alpha=1.0, p=None,
                                    algorithm='best', leaf_size=40,
                                    memory=Memory(
                                        location=None),
                                    approx_min_span_tree=True,
                                    gen_min_span_tree=False,
                                    core_dist_n_jobs=-1,
                                    cluster_selection_method='eom',
                                    allow_single_cluster=False,
                                    prediction_data=False,
                                    match_reference_implementation=False).fit(
                                        chunk[['emb1', 'emb2']]
                                        )
                                chunk['cluster'] = clustering.labels_

                                print(
                                    f'min_cluster_size = {i}',
                                    np.unique(clustering.labels_)
                                    )
                            print()

                            # acts to terminate loop
                            return min_cluster_size

                def mcs_tb_submit(text):
                    nonlocal range_tuple

                    numerical_input = text.split('.save')[0].strip()
                    new_range_tuple = tuple(map(
                        int, numerical_input.split('-')))

                    # ensure code doesn't refire on
                    # navigating away from text box
                    if new_range_tuple != range_tuple:
                        range_tuple = new_range_tuple
                        min_cluster_size = cluster_and_plot()

                    if '.save' in text:
                        if len(new_range_tuple) == 1:

                            print(
                                f'Saving current minimum cluster size: ' +
                                f'{new_range_tuple[0]}')
                            f = open(os.path.join(
                                reclass_dir, 'MCS.pkl'), 'wb')
                            pickle.dump(new_range_tuple[0], f)
                            f.close()

                            print('Saving current chunk...')
                            f = open(os.path.join(
                                reclass_dir, 'chunk.pkl'), 'wb')
                            pickle.dump(chunk, f)
                            f.close()

                            print('Saving current figure...')
                            fig.savefig(
                                os.path.join(
                                    chunk_dir,
                                    f'{self.embeddingAlgorithmQC}_'
                                    f'{new_range_tuple[0]}.png'),
                                bbox_inches='tight', dpi=1000
                                )

                            print('Proceeding with pipeline.')
                            plt.close('all')

                        else:
                            print('A single value must be entered to save an' +
                                  'optimal minimum cluster size.')

                mcs_text_box = TextBox(
                    ax_mcs_tb, "Minimum Cluster Size",
                    range_tuple, label_pad=0.1)
                mcs_text_box.on_submit(mcs_tb_submit)

                label = mcs_text_box.ax.get_children()[1]
                label.set_position([0.5, 1.3])
                label.set_horizontalalignment('center')

                def reclass_and_plot():
                    nonlocal iteration
                    nonlocal selector

                    print('Applying reclassification ' +
                          f'cutoffs: {reclass_tuple}')

                    clean = pd.DataFrame()
                    noisy = pd.DataFrame()
                    for name, cluster in chunk.groupby('cluster'):
                        if name != -1:

                            # if a cluster contains >= n% noisy data,
                            # reclassify all clustering cells as noisy
                            if (
                              (len(cluster[cluster['QC_status'] == 'clean'])
                               / len(cluster)) >= reclass_tuple[0]):
                                clean = clean.append(cluster)
                            # elif a cluster contains >= n% clean data,
                            # reclassify all clustering cells as clean
                            elif (
                              (len(cluster[cluster['QC_status'] == 'noisy'])
                               / len(cluster)) >= reclass_tuple[1]):
                                noisy = noisy.append(cluster)
                            # else keep respective cell statuses
                            else:
                                noisy = noisy.append(
                                    cluster[cluster['QC_status'] == 'noisy'])
                                clean = clean.append(
                                    cluster[cluster['QC_status'] == 'clean'])

                    # consider -1 cells from clean data to be "clean"
                    clean_outliers = chunk[
                        (chunk['cluster'] == -1)
                        & (chunk['QC_status'] == 'clean')].copy()
                    clean = clean.append(clean_outliers)

                    # consider -1 cells from noisy data to be "noisy"
                    noisy_outliers = chunk[
                        (chunk['cluster'] == -1)
                        & (chunk['QC_status'] == 'noisy')].copy()
                    noisy = noisy.append(noisy_outliers)

                    chunk['Reclass'] = 'init'
                    chunk.loc[
                        chunk.index.isin(clean.index), 'Reclass'] = 'clean'
                    chunk.loc[
                        chunk.index.isin(noisy.index), 'Reclass'] = 'noisy'

                    # PLOT embedding
                    for color_by in ['Reclass']:

                        highlight = 'none'

                        if color_by == 'Reclass':

                            # build cmap
                            cmap = ListedColormap(
                                np.array([[0.91, 0.29, 0.235],
                                          [0.18, 0.16, 0.15]]))

                            sample_dict = dict(
                                zip(
                                    natsorted(
                                        chunk['Reclass'].unique()),
                                    list(range(len(chunk['Reclass']
                                         .unique()))))
                                    )

                            c = [sample_dict[i] for
                                 i in chunk['Reclass']]

                            ax_sc_reclass.cla()

                            ax_sc_reclass.scatter(
                                chunk['emb1'],
                                chunk['emb2'],
                                c=c,
                                cmap=cmap,
                                alpha=alpha,
                                s=point_size,
                                ec=[
                                    'k' if i == highlight
                                    else 'none' for i in
                                    chunk['Reclass']
                                    ],
                                linewidth=0.1
                                )

                            ax_sc_reclass.set_title('Reclassification')
                            ax_sc_reclass.axis('equal')
                            ax_sc_reclass.tick_params(labelsize=5)
                            ax_sc_reclass.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(chunk['Reclass'].unique())
                              ):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                sample_to_map = (
                                    chunk['Sample'][
                                        chunk['Reclass'] == i]
                                    .unique()[0]
                                    )

                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none',
                                           label=i,
                                           markerfacecolor=(
                                            cmap.colors[e]
                                            ),
                                           markeredgecolor=(
                                            markeredgecolor),
                                           lw=0.001,
                                           markersize=6)
                                           )

                            reclass_lgd = ax_sc_reclass.legend(
                                handles=legend_elements,
                                prop={'size': 8}, loc='best'
                                # bbox_to_anchor=[1.0, 1.0]
                                )

                    fig.tight_layout()

                    # must call draw() to get update figure to show without
                    # first navigating mouse away from reclass text box.
                    fig.canvas.draw()
                    print('Done!')
                    print()

                    # acts to terminate loop
                    return clean, noisy

                def reclass_tb_submit(text):
                    nonlocal reclass_tuple

                    new_reclass_tuple = tuple(map(
                        float, text.split(', ')))

                    # ensure code doesn't refire on
                    # navigating away from text box
                    if new_reclass_tuple != reclass_tuple:
                        reclass_tuple = new_reclass_tuple
                        clean, noisy = reclass_and_plot()

                        # update self with current reclassified QCData slice
                        self.reclassClean = clean
                        self.reclassNoisy = noisy

                reclass_text_box = TextBox(
                    ax_reclass_tb, "ReclassCutoffs (%clean, %noisy)",
                    f'{reclass_tuple[0]}, {reclass_tuple[1]}', label_pad=0.1)
                reclass_text_box.on_submit(reclass_tb_submit)

                label = reclass_text_box.ax.get_children()[1]
                label.set_position([0.5, 1.3])
                label.set_horizontalalignment('center')

                def sample_tb_submit(text):
                    nonlocal sample

                    new_sample = str(text)

                    # ensure code doesn't refire on
                    # navigating away from text box
                    if new_sample != sample:
                        sample = new_sample
                        print(f'Setting sample-to-view to: {sample}')
                        print()

                sample_text_box = TextBox(
                    ax_sample_tb, "Sample Name",
                    f'{sample}', label_pad=0.1)
                sample_text_box.on_submit(sample_tb_submit)

                label = sample_text_box.ax.get_children()[1]
                label.set_position([0.5, 1.3])
                label.set_horizontalalignment('center')

                # kick things off with cluster_and_plot function
                min_cluster_size, clean, noisy = cluster_and_plot()

                # update self with current reclassified QCData slice
                self.reclassClean = clean
                self.reclassNoisy = noisy

                plt.show(block=True)

                # update new default MCS and reclassCutoffs
                # to the last ones entered so window reopens with
                # latest values (not original defaults)
                if len(range_tuple) == 1:
                    new_default_range_tuple = range_tuple[0]
                    new_default_reclass_tuple = reclass_tuple
                    new_default_sample = sample

                return (selector.ind, new_default_range_tuple,
                        new_default_reclass_tuple, new_default_sample)
            ###################################################################
            # loop over QCData chunks

            for chunk in chunks[chunk_index:]:

                print(
                    f'Clustering: {chunk_index + 1} of ' +
                    f'{len(chunks)} data chunks')

                # make directory for current chunk if it hasn't already been
                chunk_dir = os.path.join(reclass_dir, str(chunk_index + 1))
                if not os.path.exists(chunk_dir):
                    os.makedirs(chunk_dir)

                # if embedding for current chunk has already been computed
                # apply emb1 and emb2 to chunk dataframe
                if os.path.exists(os.path.join(chunk_dir, 'embedding.npy')):

                    # recapitulate chunk index at the point of embedding
                    chunk = chunk[~chunk['Sample'].isin(
                        self.samplesToRemoveClusteringQC)]
                    chunk = chunk.sample(
                        frac=self.fracForEmbeddingQC, random_state=5)

                    print(chunk[abx_channels_dna])
                    print()

                    embedding = np.load(
                        os.path.join(chunk_dir, 'embedding.npy'))
                    chunk['emb1'] = embedding[:, 0]
                    chunk['emb2'] = embedding[:, 1]

                else:
                    # compute embedding for chunk
                    startTime = datetime.now()

                    chunk = chunk[~chunk['Sample'].isin(
                        self.samplesToRemoveClusteringQC)]
                    chunk = chunk.sample(
                        frac=self.fracForEmbeddingQC, random_state=5)

                    print(chunk[abx_channels_dna])
                    print()

                    if self.embeddingAlgorithmQC == 'TSNE':
                        print('Computing TSNE embedding...')
                        embedding = TSNE(
                            n_components=self.dimensionEmbeddingQC,
                            perplexity=self.perplexityQC,
                            early_exaggeration=self.earlyExaggerationQC,
                            learning_rate=self.learningRateTSNEQC,
                            metric=self.metricQC,
                            random_state=self.randomStateQC,
                            init='pca', n_jobs=-1).fit_transform(
                                chunk[abx_channels_dna])

                    elif self.embeddingAlgorithmQC == 'UMAP':
                        print('Computing UMAP embedding...')
                        embedding = UMAP(
                            n_components=self.dimensionEmbeddingQC,
                            n_neighbors=self.nNeighborsQC,
                            learning_rate=self.learningRateUMAPQC,
                            output_metric=self.metricQC,
                            min_dist=self.minDistQC,
                            repulsion_strength=self.repulsionStrengthQC,
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
                            output_dens=False).fit_transform(
                                chunk[abx_channels_dna])

                    print(
                        'Embedding completed in ' +
                        str(datetime.now() - startTime)
                        )
                    print()

                    np.save(
                        os.path.join(chunk_dir, 'embedding'),
                        embedding
                        )
                    chunk['emb1'] = embedding[:, 0]
                    chunk['emb2'] = embedding[:, 1]

                # define the point size for cells in the embedding
                point_size = 50000/len(chunk)

                # interact with plots to identify optimal min cluster size
                while not os.path.isfile(os.path.join(reclass_dir, 'MCS.pkl')):

                    (selected_idx,
                     default_range_tuple,
                     default_reclass_tuple,
                     default_sample) = do_plots()

                    # if lassoed cells...
                    if len(selected_idx) >= 1:

                        # show highest expression channels
                        chunk_copy = chunk.copy()

                        # assign lassoed data a dummy cluster variable
                        # and get highest expressed markers
                        chunk_copy.loc[
                            chunk_copy.index.isin(selected_idx),
                            'cluster'] = 1000
                        hi_markers = cluster_expression(
                            df=chunk_copy, markers=abx_channels,
                            cluster=1000, num_proteins=3,
                            standardize='within')

                        print(f'Top three markers {hi_markers}')
                        print()

                        # superimpose centroids of lassoed noisy cells
                        # colored by stage removed over channel images
                        print('Opening Napari...')
                        print()

                        napari_warnings()
                        if self.showAbChannels:
                            for e, ch in enumerate(reversed(abx_channels)):
                                channel_number = marker_channel_number(
                                    markers, ch)

                                # read antibody image
                                for file_path in glob.glob(
                                  f'{self.inDir}/tif/{default_sample}*.tif'):
                                    img = single_channel_pyramid(
                                        file_path,
                                        channel=channel_number.item() - 1
                                        )

                                # initialize Napari viewer with first channel
                                if e == 0:
                                    viewer = napari.view_image(
                                        img, rgb=False,
                                        blending='additive',
                                        colormap='green',
                                        visible=False,
                                        name=ch)
                                else:
                                    viewer.add_image(
                                        img, rgb=False,
                                        blending='additive',
                                        colormap='green',
                                        visible=False,
                                        name=ch)

                        # color noisy data points by module used to redact them
                        cmap = categorical_cmap(
                            numUniqueSamples=len(modules[1:]),
                            numCatagories=10,
                            cmap='tab10',
                            continuous=False
                            )

                        # reverse module order so they appear in
                        # correct order in Napari
                        QC_color_dict = dict(zip(modules[1:], cmap.colors))

                        for module, color in reversed(QC_color_dict.items()):

                            centroids = chunk[['Y_centroid', 'X_centroid']][
                                    (chunk.index.isin(selected_idx))
                                    & (chunk['Sample'] == default_sample)
                                    & (chunk['QC_status'] == 'noisy')
                                    & (chunk['Filter'] == module)]

                            viewer.add_points(
                                centroids, name=module, visible=False,
                                face_color=color, edge_width=0.0, size=4.0)

                        # read segmentation outlines, add to Napari
                        for file_path in glob.glob(
                          f'{self.inDir}/seg/{default_sample}*.ome.tif'):
                            seg = single_channel_pyramid(
                                file_path, channel=0)
                        viewer.add_image(
                            seg, rgb=False,
                            blending='additive',
                            colormap='gray',
                            visible=False,
                            name='segmentation')

                        # read last DNA, add to Napari
                        last_dna_cycle = natsorted(
                            [i for i in chunk.columns if dna_moniker in i])[-1]
                        channel_number = marker_channel_number(
                            markers, last_dna_cycle)
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{default_sample}*.tif'):
                            dna_last = single_channel_pyramid(
                                file_path,
                                channel=channel_number.item() - 1)
                            viewer.add_image(
                                dna_last, rgb=False, blending='additive',
                                opacity=0.5, colormap='gray', visible=False,
                                name=f'{last_dna_cycle}: ' +
                                f'{default_sample}')

                        # read first DNA, add to Napari
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{default_sample}*.tif'):
                            dna_first = single_channel_pyramid(
                                file_path, channel=0)
                        viewer.add_image(
                            dna_first, rgb=False, blending='additive',
                            opacity=0.5, colormap='gray', visible=True,
                            name=f'{dna1}: ' +
                            f'{default_sample}')

                        napari.run()

                ###############################################################
                # once optimal MCS has been saved

                if os.path.exists(os.path.join(reclass_dir, 'MCS.pkl')):
                    f = open(os.path.join(
                        reclass_dir, 'MCS.pkl'), 'rb')
                    min_cluster_size = pickle.load(f)

                    ###########################################################
                    # apply last entered reclassification cutoffs

                    print(
                        f'Applying current reclassification ' +
                        f'cutoffs: {default_reclass_tuple}')
                    print()

                    # update clean and noisy storage dataframes and save
                    reclass_storage_dict['clean'] = (
                        reclass_storage_dict['clean'].append(
                            self.reclassClean))
                    reclass_storage_dict['noisy'] = (
                        reclass_storage_dict['noisy'].append(
                            self.reclassNoisy))

                    f = open(os.path.join(
                        reclass_dir, 'reclass_storage_dict.pkl'), 'wb')
                    pickle.dump(reclass_storage_dict, f)
                    f.close()

                    print('Reclassified clean tally: ' +
                          f"{len(reclass_storage_dict['clean'])}")
                    print('Reclassified noisy tally: ' +
                          f"{len(reclass_storage_dict['noisy'])}")
                    print()

                    ###########################################################
                    # cluster chunk using selected MCS
                    # (not applicable to first chunk, which gets
                    # clustered during plotting above)

                    if 'cluster' not in chunk.columns:

                        print(
                            'Applying current minimum cluster size: ' +
                            f'{min_cluster_size}')
                        print()

                        clustering = hdbscan.HDBSCAN(
                            min_cluster_size=min_cluster_size,
                            min_samples=None,
                            metric='euclidean', alpha=1.0, p=None,
                            algorithm='best', leaf_size=40,
                            memory=Memory(
                                location=None),
                            approx_min_span_tree=True,
                            gen_min_span_tree=False,
                            core_dist_n_jobs=-1,
                            cluster_selection_method='eom',
                            allow_single_cluster=False,
                            prediction_data=False,
                            match_reference_implementation=False).fit(
                                chunk[['emb1', 'emb2']]
                                )
                        chunk['cluster'] = clustering.labels_

                    ###########################################################
                    # get clustermap

                    # exit program if all cells are considered ambiguous by the
                    # clustering algorithm (likely too few cells per chunk)
                    if chunk['cluster'].eq(-1).all():
                        print(
                            f'WARNING: All cells in chunk {chunk_index + 1} ' +
                            'were deemed ambiguous by clustering algorithm ' +
                            '(i.e. assigned to cluster -1), exiting program. '
                            + 'Try using a larger batch size.')
                        exit()

                    clustermap_input = chunk[chunk['cluster'] != -1]

                    cluster_heatmap_input = clustermap_input[
                        abx_channels + ['cluster']].groupby('cluster').mean()

                    sns.set(font_scale=0.8)
                    for name, axis in zip(['within', 'across'], [0, 1]):
                        g = sns.clustermap(
                            cluster_heatmap_input, cmap='viridis',
                            standard_scale=axis, square=False, yticklabels=1,
                            linewidth=0.1, cbar=True
                            )

                        plt.gcf().set_size_inches(8.0, 8.0)

                        plt.savefig(
                            os.path.join(
                                chunk_dir, f'clustermap_{name}.pdf'),
                            bbox_inches='tight')
                        plt.close('all')

                    ###########################################################
                    # increment chunk_index
                    chunk_index = chunk_index + 1
                    f = open(os.path.join(
                        reclass_dir, 'chunk_index.pkl'), 'wb')
                    pickle.dump(chunk_index, f)
                    f.close()

                    # remove saved chunk.pkl
                    if os.path.exists(os.path.join(reclass_dir, 'chunk.pkl')):
                        os.remove(os.path.join(reclass_dir, 'chunk.pkl'))

            ###################################################################
            # perform data reclassification

            # create explicit global labels for raw data
            for module_idx, module in enumerate(modules):
                if module_dict[module_idx][0] == 'aggregateData':
                    pre_qc = module_dict[module_idx][1].copy()
                    pre_qc['handle'] = (
                        pre_qc['CellID'].map(str) + '_' +
                        pre_qc['Sample'])

            # create explicit global labels for
            # cleaned data (before reclassifiction)
            post_qc = module_dict[
                [i for i in module_dict.keys()][-1]][1].copy()
            post_qc['handle'] = (
                post_qc['CellID'].map(str) + '_' +
                post_qc['Sample'])

            # get raw values of cells in post_qc data
            cleaned_raw = pre_qc[pre_qc['handle'].isin(post_qc['handle'])]

            # convert clean data in predominantly noisy clusters to noisy
            # to yield final clean data
            drop = reclass_storage_dict['noisy'][
                reclass_storage_dict['noisy']['QC_status'] == 'clean'].copy()
            drop['handle'] = (
                drop['CellID'].map(str) + '_' +
                drop['Sample'])
            dropped = cleaned_raw[~cleaned_raw['handle'].isin(drop['handle'])]

            # convert noisy data in predominantly clean clusters to clean
            # to yield final replace data
            replace = reclass_storage_dict['clean'][
                reclass_storage_dict['clean']['QC_status'] == 'noisy'].copy()
            replace['handle'] = (
                replace['CellID'].map(str) + '_' +
                replace['Sample'])
            replaced = pre_qc[pre_qc['handle'].isin(replace['handle'])]

            data = dropped.append(replaced)

            ###################################################################
            # transform data

            data.loc[:, abx_channels] += 0.00000000001
            data.loc[:, abx_channels] = np.log10(data[abx_channels])

            ###################################################################
            # rescale antibody signal intensities (0-1)

            for ab in abx_channels:
                channel_data = data[ab]

                # rescale channel signal intensities
                scaler = (
                    MinMaxScaler(feature_range=(0, 1), copy=True)
                    .fit(channel_data.values.reshape(-1, 1))
                        )
                rescaled_data = scaler.transform(
                    channel_data.values.reshape(-1, 1)
                    )

                rescaled_data = pd.DataFrame(
                    data=rescaled_data,
                    index=channel_data.index,
                    ).rename(columns={0: ab})

                # data to return from metQC module
                data.update(rescaled_data)

        #######################################################################
        # compute number of cells remaining after each QC stage.

        # Consider reclassified data if self.metaQC is True.
        # Otherwise, only consider originally filtered data.

        qc_dict = {}
        for module_idx, module in enumerate(modules):

            mod_name = module_dict[module_idx][0]
            mod_data = len(module_dict[module_idx][1])

            if mod_name not in ['aggregateData', 'selectROIs']:
                if self.metaQC:
                    qc_dict[mod_name] = mod_data + len(
                        replace[replace['Filter'] == mod_name])
                else:
                    qc_dict[mod_name] = mod_data

            elif mod_name == 'selectROIs':
                if self.metaQC:
                    if self.delintMode:
                        qc_dict[mod_name] = mod_data + len(
                            replace[replace['Filter'] == mod_name])
                    else:
                        qc_dict[mod_name] = mod_data
                else:
                    qc_dict[mod_name] = mod_data

            elif mod_name == 'aggregateData':
                qc_dict[mod_name] = mod_data

        if self.metaQC:
            meta_qc = len(drop[drop['QC_status'] == 'clean'])
            qc_dict['metaQC'] = meta_qc
        else:
            qc_dict['metaQC'] = 0

        # using qc_dict, compute percentage of total data
        # redacted at each stage.
        percents = {}
        qc_keys = list(qc_dict)
        for i in range(len(modules) - 1):

            percent = (
                (qc_dict[qc_keys[i]]
                 - qc_dict[qc_keys[i + 1]]) / qc_dict[qc_keys[0]]) * 100

            percents[qc_keys[i + 1]] = percent

        # add percent of data redacted during metaQC
        percents['metaQC'] = (
            (qc_dict[qc_keys[-1]]) / qc_dict[qc_keys[0]]) * 100

        # compute fraction of data passing all QC stages (residual data)
        percents['residual'] = (100 - sum(percents.values()))

        # plot pie chart of percent data redacted at each QC stage
        labels = list(percents.keys())
        sizes = list(percents.values())
        explode = [0.0] * len(sizes)

        theme = plt.get_cmap('Set2')

        fig1, ax1 = plt.subplots()
        ax1.set_prop_cycle('color', [theme(1.0 * i / len(sizes))
                                     for i in range(len(sizes))])
        patches, texts, autotexts = ax1.pie(
            sizes, explode=explode, labels=labels, autopct='%1.2f%%',
            shadow=False, startangle=90
            )

        for i in range(len(texts)):
            texts[i].set_fontsize(4)
            autotexts[i].set_fontsize(4)

        ax1.axis('equal')

        plt.savefig(
            os.path.join(
                reclass_dir, 'censored_by_stage.pdf'), bbox_inches='tight'
            )
        plt.close('all')

        # drop indexing handle from data before returning
        if self.metaQC:
            data.drop('handle', axis=1, inplace=True)

        print()
        print()
        return data

    @module
    def PCA(data, self, args):

        if len(data['Sample'].unique()) > 1:

            # read marker metadata
            markers, dna1, dna_moniker, abx_channels = read_markers(
                markers_filepath=os.path.join(self.inDir, 'markers.csv'),
                markers_to_exclude=self.markersToExclude,
                )

            # drop antibody channel exclusions for PCA
            abx_channels = [
                i for i in abx_channels if i not in self.channelExclusionsPCA]

            # create PCA directory if it doesn't already exist
            pca_dir = os.path.join(self.outDir, f'PCA')
            if not os.path.exists(pca_dir):
                os.makedirs(pca_dir)

            # compute median antibody expression per sample
            # samples (rows) x features (columns)
            medians = data.groupby(['Sample']).median()[abx_channels]

            # drop sample exclusions for PCA
            medians = medians[~medians.index.isin(self.samplesToRemovePCA)]

            # mean-center data
            medians = medians-medians.mean(axis=0)

            # sort medians index (sample names) naturally
            medians = medians.reindex(natsorted(medians.index))

            # initialize PCA
            pca = PCA(self.dimensionPCA, random_state=1)

            # fit PCA to data
            projected = pca.fit_transform(medians)

            # generate dataframe of the projection coordinates
            # to be used as plot input
            scatter_input = pd.DataFrame(data=projected, index=medians.index)
            scatter_input.rename(columns={0: 'PC1', 1: 'PC2'}, inplace=True)

            # assign row index (sample names) as a column
            scatter_input.reset_index(drop=False, inplace=True)

            # get sample file names (i.e. sampleMetadata keys) from config.yml
            # based on "Sample" column (first elements of sampleMetadata vals)
            def get_key(val):
                for key, value in self.sampleNames.items():
                    if val == value:
                        return key

                return "key doesn't exist"
            file_names = [get_key(i) for i in scatter_input['Sample']]

            # create columns of full and abbreviated condition names
            scatter_input['condition'] = [
                self.sampleConditions[i] for i in file_names]
            scatter_input['condition_abbr'] = [
                self.sampleConditionAbbrs[i] for i in file_names]

            # set condition abbreviation names as row index
            scatter_input.set_index('condition_abbr', inplace=True)

            # create a naturally-sorted list of samples that will NOT be grayed
            ordered_conditions = natsorted(set(scatter_input.index.unique())
                                           .difference(
                                            set(self.conditionsToSilhouette)))

            # build cmap for samples to NOT be grayed
            cmap = categorical_cmap(
                numUniqueSamples=len(ordered_conditions),
                numCatagories=10,
                cmap='tab10',
                continuous=False)
            condition_color_dict = dict(zip(ordered_conditions, cmap.colors))
            # assign conditions for silhouetting the same color (e.g. gray)
            for s in self.conditionsToSilhouette:
                silhouette_color = [0.863, 0.863, 0.863]
                condition_color_dict[s] = silhouette_color

            # plot sample PCA scores
            sns.set_style('white')
            for e, (condition_name, sample_scores) in enumerate(
              scatter_input.iterrows()):

                point = pd.DataFrame(sample_scores).T

                if condition_name in self.conditionsToSilhouette:
                    edgecolor = silhouette_color
                    zorder = 2
                else:
                    edgecolor = 'k'
                    zorder = 3

                g = sns.scatterplot(
                    data=point, x='PC1', y='PC2', hue=point.index,
                    palette=condition_color_dict, edgecolor=edgecolor,
                    linewidth=0.2, zorder=zorder, s=self.pointSize,
                    alpha=1.0, legend=False)

            g.grid(color='gray', linewidth=0.05, linestyle='-', alpha=1.0)
            plt.setp(g.spines.values(), color='k', lw=0.5)

            # assign row index (condition abbreviations) as column
            scatter_input = scatter_input.reset_index().rename(
                columns={'index': 'abbreviation'}
                )

            # annotate data points
            if self.labelPoints is True:

                # generate squareform distance matrix
                sq = squareform(
                    pdist(scatter_input[['PC1', 'PC2']], metric='euclidean')
                    )

                # add numerical row and column indices
                df = pd.DataFrame(
                    sq, index=scatter_input.index,
                    columns=scatter_input.index
                    )

                # isolate values from upper triangle
                df1 = df.where(np.triu(np.ones(df.shape)).astype(np.bool))

                # filter upper triangle values according to distance cutoff
                df2 = df1[df1 < self.distanceCutoff]

                # flatten, set multi-index as columns (sample1, sample2)
                df3 = (
                    df2
                    .stack()
                    .reset_index()
                    .rename(columns={'level_0': 'sample_id1',
                                     'level_1': 'sample_id2',
                                     0: 'dist'})
                    )

                # get condition abbreviations for sample pairs
                df3['sample_id1_cond'] = [
                    scatter_input.loc[i, 'condition_abbr']
                    for i in df3['sample_id1']]
                df3['sample_id2_cond'] = [
                    scatter_input.loc[i, 'condition_abbr']
                    for i in df3['sample_id2']]

                # drop diagonal values
                df4 = df3[df3['dist'] != 0.0].dropna()

                # get proximal sample pairs of the same condition
                df5 = df4[df4['sample_id1_cond'] == df4['sample_id2_cond']]

                # create a set of proximal sample indices
                proximal_indices = set(
                    list(df5['sample_id1']) +
                    list(df5['sample_id2']))

                # set of conditions
                unique_conds = set()
                # set of neighbors
                neighbors_set = set()

                # loop over scatter_input to annoate plot points
                for e, (cond, x, y) in enumerate(zip(
                  scatter_input['condition_abbr'],
                  scatter_input['PC1'], scatter_input['PC2'])):

                    # if data point e has an associated condition
                    # which is not to be grayed out
                    if cond not in self.conditionsToSilhouette:
                        # if data point e has at least one neighbor
                        if e in proximal_indices:
                            # and hasn't already been accounted for as
                            # another sample's neighbor
                            if e not in neighbors_set:
                                # get data for all samples
                                # neighboring data point e
                                df6 = (
                                    df5[(df5['sample_id1'] == e) |
                                        (df5['sample_id2'] == e)]
                                    [['sample_id1', 'sample_id2']])

                                # get set of all indices
                                # neighboring data point e
                                neighbors = set(
                                    list(df6['sample_id1']) +
                                    list(df6['sample_id2']))

                                # add neighboring indices to overall
                                # neighbors_set
                                neighbors_set = neighbors_set.union(neighbors)

                                # slice scatter_input to get samples
                                # proximal to data point e
                                neighbors_df = scatter_input.loc[neighbors]

                                # generate centroid between samples
                                # neighboring data point e
                                pc1 = neighbors_df['PC1']
                                pc2 = neighbors_df['PC2']
                                centroid = (
                                    sum(pc1) / len(neighbors_df),
                                    sum(pc2) / len(neighbors_df))
                                x = centroid[0]
                                y = centroid[1]
                            else:
                                x = None
                                y = None
                        else:
                            pass

                        # if data point e has already been accounted
                        # for by a neighborhood centroid, do nothing
                        if x == y is None:
                            pass
                        else:
                            # annotate plot at the centroid between
                            # neighboring samples of the same condition
                            # or annotate data point itself
                            text = plt.annotate(
                                cond, xy=(x, y),
                                xytext=(0, 0), size=4.75,
                                fontweight='bold',
                                color=condition_color_dict[cond],
                                textcoords='offset points', ha='center',
                                va='center',
                                bbox=dict(boxstyle='round,pad=0.1',
                                          fc='yellow', alpha=0.0))

                            text.set_path_effects(
                                [path_effects.Stroke(
                                    linewidth=0.75, foreground='k'),
                                 path_effects.Normal()]
                                )

                            # add condition to set of conditions
                            unique_conds.add(cond)

            # get n per condition in naturally sorted order
            lgd_data = scatter_input[
                ['condition_abbr', 'condition']].value_counts()
            lgd_data = lgd_data.reindex(index=natsorted(lgd_data.index))

            # loop over legend data to create legend handles
            legend_handles = []
            for i in lgd_data.items():
                abbr = i[0][0]
                cond = i[0][1]
                n = i[1]

                if abbr in self.conditionsToSilhouette:
                    markerfacecolor = 'gainsboro'
                    markeredgecolor = 'gainsboro'
                else:
                    markerfacecolor = condition_color_dict[abbr]
                    markeredgecolor = 'k'

                legend_handles.append(
                    Line2D([0], [0], marker='o', color='none',
                           label=f'{abbr} ({cond}, n={n})',
                           markerfacecolor=markerfacecolor,
                           markeredgecolor=markeredgecolor,
                           markeredgewidth=0.2,
                           markersize=5.0))

            # add legend to plot
            g.legend(handles=legend_handles, prop={'size': 5.0}, loc='best')

            # update x and y axis labels
            plt.xlabel(
                f'PC1 ({round((pca.explained_variance_ratio_[0] * 100), 2)}'
                '% of variance)', fontsize=10, labelpad=7.0)
            plt.ylabel(
                f'PC2 ({round((pca.explained_variance_ratio_[1] * 100), 2)}'
                '% of variance)', fontsize=10, labelpad=4.0)

            # modify x and y axis ticks
            plt.tick_params(axis='both', which='major', labelsize=7.0)

            # save figure
            plt.savefig(
                os.path.join(pca_dir, 'pcaScoresPlot.pdf'),
                bbox_inches='tight')
            plt.close('all')

        print()
        print()
        return data

    @module
    def clustering(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude)

        # drop antibody channel exclusions for clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering]

        # create clustering directory if it hasn't already
        clustering_dir = os.path.join(self.outDir, 'clustering/final')
        if not os.path.exists(clustering_dir):
            os.makedirs(clustering_dir)

        # recapitulate df index at the point of embedding
        data = data[~data['Sample'].isin(self.samplesToRemoveClustering)]

        # assign default HDBSCAN minimum cluster size
        default_range_tuple = 250

        # assign default sample for viz
        default_sample = self.viewSample

        # define the degree of transparency for cells in the embedding
        alpha = 1.0

        # pick a random seed for reproducibility of
        # data sampling via self.fracForEmbedding
        random_state = 5

        if self.normalizeTissueCounts:
            print('Performing weighted sampling of cells from tissues...')
            print(
                'Check that resulting cell counts are similar across samples.')
            print('If not, try embedding a smaller fraction of data.')
            print()

            # calculate per tissue cell-count weighted random sample
            groups = data.groupby('Sample')

            sample_weights = pd.DataFrame({
                'weights': 1 / (groups.size() * len(groups))})

            weights = pd.merge(
                data[['Sample']], sample_weights,
                left_on='Sample', right_index=True)

            data = data.sample(
                frac=self.fracForEmbedding, replace=False,
                weights=weights['weights'], random_state=random_state,
                axis=0)

            print('Cell counts:')
            print(data.groupby(['Sample']).size())
            print()

        else:

            data = data.sample(
                frac=self.fracForEmbedding, random_state=random_state)

        data.reset_index(drop=True, inplace=True)
        print(data[abx_channels])

        # if embedding has already been computed
        if os.path.exists(os.path.join(clustering_dir, 'embedding.npy')):

            embedding = np.load(os.path.join(clustering_dir, 'embedding.npy'))
            data['emb1'] = embedding[:, 0]
            data['emb2'] = embedding[:, 1]

        else:
            startTime = datetime.now()

            if self.embeddingAlgorithm == 'TSNE':
                print('Computing TSNE embedding.')
                embedding = TSNE(
                    n_components=self.dimensionEmbedding,
                    perplexity=self.perplexity,
                    early_exaggeration=self.earlyExaggeration,
                    learning_rate=self.learningRateTSNE,
                    metric=self.metric,
                    random_state=self.randomState,
                    init='pca', n_jobs=-1).fit_transform(data[abx_channels])

            elif self.embeddingAlgorithm == 'UMAP':
                print('Computing UMAP embedding.')
                embedding = UMAP(
                    n_components=self.dimensionEmbedding,
                    n_neighbors=self.nNeighbors,
                    learning_rate=self.learningRateUMAP,
                    output_metric=self.metric,
                    min_dist=self.minDist,
                    repulsion_strength=self.repulsionStrength,
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
                    set_op_mix_ratio=0.5,
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
                    output_dens=False).fit_transform(data[abx_channels])

            print('Embedding completed in ' + str(datetime.now() - startTime))
            print()

            np.save(os.path.join(clustering_dir, 'embedding'), embedding)
            data['emb1'] = embedding[:, 0]
            data['emb2'] = embedding[:, 1]

        #######################################################################

        def do_plots():

            # initialize counter to allow for different
            # behavior on the first round of plotting
            iteration = 1

            sns.set_style('whitegrid')

            fig = plt.figure(figsize=(17, 6))
            matplotlib_warnings(fig)

            gs = plt.GridSpec(2, 3, figure=fig, height_ratios=[1, 10])

            # define text box axes
            ax_mcs_tb = fig.add_subplot(gs[0, 0])
            ax_sample_tb = fig.add_subplot(gs[0, 1])

            # define embedding axes
            ax_sc_cluster = fig.add_subplot(gs[1, 0])
            ax_sc_sample = fig.add_subplot(gs[1, 1])
            ax_sc_cond = fig.add_subplot(gs[1, 2])

            # assign default inputs
            range_tuple = default_range_tuple
            sample = default_sample
            selector = None

            def cluster_and_plot():
                nonlocal iteration
                nonlocal selector
                nonlocal range_tuple

                if not isinstance(range_tuple, tuple):
                    range_tuple = (range_tuple, )

                # loop for single-value MCS inputs
                if iteration == 1:
                    if len(range_tuple) == 1:

                        mylist = [range_tuple[0]]

                        for i in mylist:

                            min_cluster_size = i

                            clustering = hdbscan.HDBSCAN(
                                min_cluster_size=min_cluster_size,
                                min_samples=None,
                                metric='euclidean', alpha=1.0, p=None,
                                algorithm='best', leaf_size=40,
                                memory=Memory(
                                    location=None),
                                approx_min_span_tree=True,
                                gen_min_span_tree=False,
                                core_dist_n_jobs=-1,
                                cluster_selection_method='eom',
                                allow_single_cluster=False,
                                prediction_data=False,
                                match_reference_implementation=False).fit(
                                    data[['emb1', 'emb2']]
                                    )
                            data['cluster'] = clustering.labels_

                            print(
                                f'min_cluster_size = {i}',
                                np.unique(clustering.labels_)
                                )
                            print()

                            # PLOT embedding
                            for color_by in ['cluster', 'Sample', 'Condition']:

                                highlight = 'none'

                                if color_by == 'cluster':

                                    # build cmap
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(
                                            data[color_by].unique()),
                                        numCatagories=10,
                                        cmap='tab10',
                                        continuous=False
                                        )

                                    # make black the first color to specify
                                    # unclustered cells (cluster -1)
                                    cmap = ListedColormap(
                                        np.insert(
                                            arr=cmap.colors, obj=0,
                                            values=[0, 0, 0], axis=0)
                                            )

                                    # trim cmap to # unique samples
                                    trim = (
                                        len(cmap.colors) - len(
                                            data[color_by].unique())
                                        )
                                    cmap = ListedColormap(
                                        cmap.colors[:-trim]
                                        )

                                    sample_dict = dict(
                                        zip(
                                            natsorted(
                                                data[color_by].unique()),
                                            list(range(len(data[color_by]
                                                 .unique()))))
                                            )

                                    c = [sample_dict[i] for i
                                         in data[color_by]]

                                    ax_sc_cluster.cla()

                                    cluster_paths = ax_sc_cluster.scatter(
                                        data['emb1'],
                                        data['emb2'],
                                        c=c,
                                        alpha=alpha,
                                        s=point_size,
                                        cmap=cmap,
                                        ec=[
                                            'k' if i == highlight
                                            else 'none' for i in
                                            data[color_by]
                                            ],
                                        linewidth=0.1
                                        )

                                    ax_sc_cluster.set_title('HDBSCAN')
                                    ax_sc_cluster.axis('equal')
                                    ax_sc_cluster.tick_params(labelsize=5)
                                    ax_sc_cluster.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(
                                        natsorted(data[color_by].unique())
                                      ):

                                        hi_markers = cluster_expression(
                                            df=data, markers=abx_channels,
                                            cluster=i, num_proteins=3,
                                            standardize='within',
                                            )

                                        legend_elements.append(
                                            Line2D([0], [0], marker='o',
                                                   color='none',
                                                   label=(
                                                    f'Cluster {i}: '
                                                    f'{hi_markers}'),
                                                   markerfacecolor=(
                                                    cmap.colors[e]),
                                                   markeredgecolor='none',
                                                   lw=0.001, markersize=6))

                                    # cluster_lgd = ax_sc_cluster.legend(
                                    #     handles=legend_elements,
                                    #     prop={'size': 8},
                                    #     bbox_to_anchor=[1.0, 1.0]
                                    #     )

                                elif color_by == 'Sample':

                                    # build cmap
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(
                                            data['Sample'].unique()),
                                        numCatagories=10,
                                        cmap='tab10',
                                        continuous=False
                                        )

                                    sample_dict = dict(
                                        zip(
                                            natsorted(
                                                data['Sample'].unique()),
                                            list(range(len(data['Sample']
                                                 .unique()))))
                                            )

                                    c = [sample_dict[i] for
                                         i in data['Sample']]

                                    ax_sc_sample.cla()

                                    ax_sc_sample.scatter(
                                        data['emb1'],
                                        data['emb2'],
                                        c=c,
                                        cmap=cmap,
                                        alpha=alpha,
                                        s=point_size,
                                        ec=[
                                            'k' if i == highlight
                                            else 'none' for i in
                                            data['Sample']
                                            ],
                                        linewidth=0.1
                                        )
                                    ax_sc_sample.set_title('Sample')

                                    ax_sc_sample.axis('equal')
                                    ax_sc_sample.tick_params(labelsize=5)
                                    ax_sc_sample.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(
                                      natsorted(data['Sample'].unique())):

                                        if i == highlight:
                                            markeredgecolor = 'k'
                                        else:
                                            markeredgecolor = 'none'

                                        sample_to_map = (
                                            data['Sample'][
                                                data['Sample'] == i]
                                            .unique()[0]
                                            )

                                        legend_elements.append(
                                            Line2D([0], [0], marker='o',
                                                   color='none',
                                                   label=i,
                                                   markerfacecolor=(
                                                    cmap.colors[e]
                                                    ),
                                                   markeredgecolor=(
                                                    markeredgecolor),
                                                   lw=0.001,
                                                   markersize=6)
                                                   )

                                    sample_lgd = ax_sc_sample.legend(
                                        handles=legend_elements,
                                        prop={'size': 8}, loc='best')

                                elif color_by == 'Condition':

                                    # build cmap
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(
                                            data['Condition'].unique()),
                                        numCatagories=10,
                                        cmap='tab10',
                                        continuous=False
                                        )

                                    sample_dict = dict(
                                        zip(
                                            natsorted(
                                                data['Condition'].unique()),
                                            list(range(len(data['Condition']
                                                 .unique()))))
                                            )

                                    c = [sample_dict[i] for
                                         i in data['Condition']]

                                    ax_sc_cond.cla()

                                    ax_sc_cond.scatter(
                                        data['emb1'],
                                        data['emb2'],
                                        c=c,
                                        cmap=cmap,
                                        alpha=alpha,
                                        s=point_size,
                                        ec=[
                                            'k' if i == highlight
                                            else 'none' for i in
                                            data['Condition']
                                            ],
                                        linewidth=0.1
                                        )
                                    ax_sc_cond.set_title('Condition')

                                    ax_sc_cond.axis('equal')
                                    ax_sc_cond.tick_params(labelsize=5)
                                    ax_sc_cond.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(
                                        natsorted(data['Condition'].unique())
                                      ):

                                        if i == highlight:
                                            markeredgecolor = 'k'
                                        else:
                                            markeredgecolor = 'none'

                                        sample_to_map = (
                                            data['Condition'][
                                                data['Condition'] == i]
                                            .unique()[0]
                                            )

                                        legend_elements.append(
                                            Line2D([0], [0], marker='o',
                                                   color='none',
                                                   label=i,
                                                   markerfacecolor=(
                                                    cmap.colors[e]
                                                    ),
                                                   markeredgecolor=(
                                                    markeredgecolor),
                                                   lw=0.001,
                                                   markersize=6)
                                                   )

                                    sample_lgd = ax_sc_cond.legend(
                                        handles=legend_elements,
                                        prop={'size': 8}, loc='best')

                            fig.tight_layout()

                            # must call draw() before creating selector,
                            # or alpha setting doesn't work.
                            fig.canvas.draw()
                            print('Done!')
                            print()

                            if selector:
                                selector.disconnect()
                            selector = SelectFromCollection(
                                ax_sc_cluster, cluster_paths
                                )

                            # allow for different behavior
                            # on first iteration of plotting
                            iteration += 1

                            return min_cluster_size

                    else:
                        # loop for a range of input values
                        mylist = list(
                            range(range_tuple[0], range_tuple[1] + 1, 1)
                            )
                        mylist.reverse()  # run higher cluster sizes first

                        for i in mylist:

                            min_cluster_size = i

                            clustering = hdbscan.HDBSCAN(
                                min_cluster_size=min_cluster_size,
                                min_samples=None,
                                metric='euclidean', alpha=1.0, p=None,
                                algorithm='best', leaf_size=40,
                                memory=Memory(
                                    location=None),
                                approx_min_span_tree=True,
                                gen_min_span_tree=False,
                                core_dist_n_jobs=-1,
                                cluster_selection_method='eom',
                                allow_single_cluster=False,
                                prediction_data=False,
                                match_reference_implementation=False).fit(
                                    data[['emb1', 'emb2']]
                                    )
                            data['cluster'] = clustering.labels_

                            print(
                                f'min_cluster_size = {i}',
                                np.unique(clustering.labels_)
                                )
                        print()

                        # acts to terminate loop
                        return min_cluster_size

                if iteration != 1:
                    if len(range_tuple) == 1:

                        mylist = [range_tuple[0]]

                        for i in mylist:

                            min_cluster_size = i

                            clustering = hdbscan.HDBSCAN(
                                min_cluster_size=min_cluster_size,
                                min_samples=None,
                                metric='euclidean', alpha=1.0, p=None,
                                algorithm='best', leaf_size=40,
                                memory=Memory(
                                    location=None),
                                approx_min_span_tree=True,
                                gen_min_span_tree=False,
                                core_dist_n_jobs=-1,
                                cluster_selection_method='eom',
                                allow_single_cluster=False,
                                prediction_data=False,
                                match_reference_implementation=False).fit(
                                    data[['emb1', 'emb2']]
                                    )
                            data['cluster'] = clustering.labels_

                            print(
                                f'min_cluster_size = {i}',
                                np.unique(clustering.labels_)
                                )
                            print()

                            # PLOT embedding
                            for color_by in ['cluster']:

                                highlight = 'none'

                                if color_by == 'cluster':

                                    # build cmap
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(
                                            data[color_by].unique()),
                                        numCatagories=10,
                                        cmap='tab10',
                                        continuous=False
                                        )

                                    # make black the first color to specify
                                    # unclustered cells (cluster -1)
                                    cmap = ListedColormap(
                                        np.insert(
                                            arr=cmap.colors, obj=0,
                                            values=[0, 0, 0], axis=0)
                                            )

                                    # trim cmap to # unique samples
                                    trim = (
                                        len(cmap.colors) - len(
                                            data[color_by].unique())
                                        )
                                    cmap = ListedColormap(
                                        cmap.colors[:-trim]
                                        )

                                    sample_dict = dict(
                                        zip(
                                            natsorted(
                                                data[color_by].unique()),
                                            list(range(len(data[color_by]
                                                 .unique()))))
                                            )

                                    c = [sample_dict[i] for i
                                         in data[color_by]]

                                    ax_sc_cluster.cla()

                                    cluster_paths = ax_sc_cluster.scatter(
                                        data['emb1'],
                                        data['emb2'],
                                        c=c,
                                        alpha=alpha,
                                        s=point_size,
                                        cmap=cmap,
                                        ec=[
                                            'k' if i == highlight
                                            else 'none' for i in
                                            data[color_by]
                                            ],
                                        linewidth=0.1
                                        )

                                    ax_sc_cluster.set_title('HDBSCAN')
                                    ax_sc_cluster.axis('equal')
                                    ax_sc_cluster.tick_params(labelsize=5)
                                    ax_sc_cluster.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(
                                      natsorted(data[color_by].unique())):

                                        hi_markers = cluster_expression(
                                            df=data, markers=abx_channels,
                                            cluster=i, num_proteins=3,
                                            standardize='within',
                                            )

                                        legend_elements.append(
                                            Line2D([0], [0], marker='o',
                                                   color='none',
                                                   label=(
                                                    f'Cluster {i}: '
                                                    f'{hi_markers}'),
                                                   markerfacecolor=(
                                                    cmap.colors[e]),
                                                   markeredgecolor='none',
                                                   lw=0.001, markersize=6))

                                    # cluster_lgd = ax_sc_cluster.legend(
                                    #     handles=legend_elements,
                                    #     prop={'size': 8},
                                    #     bbox_to_anchor=[1.0, 1.0]
                                    #     )

                            fig.tight_layout()

                            # Must call draw() before creating selector,
                            # or alpha setting doesn't work.
                            fig.canvas.draw()

                            if selector:
                                selector.disconnect()
                            selector = SelectFromCollection(
                                ax_sc_cluster, cluster_paths
                                )

                            # acts to terminate loop
                            return min_cluster_size

                    else:
                        # loop for a range of input values
                        mylist = list(
                            range(range_tuple[0], range_tuple[1] + 1, 1)
                            )
                        mylist.reverse()  # run higher cluster sizes first

                        for i in mylist:

                            min_cluster_size = i

                            clustering = hdbscan.HDBSCAN(
                                min_cluster_size=min_cluster_size,
                                min_samples=None,
                                metric='euclidean', alpha=1.0, p=None,
                                algorithm='best', leaf_size=40,
                                memory=Memory(
                                    location=None),
                                approx_min_span_tree=True,
                                gen_min_span_tree=False,
                                core_dist_n_jobs=-1,
                                cluster_selection_method='eom',
                                allow_single_cluster=False,
                                prediction_data=False,
                                match_reference_implementation=False).fit(
                                    data[['emb1', 'emb2']]
                                    )
                            data['cluster'] = clustering.labels_

                            print(
                                f'min_cluster_size = {i}',
                                np.unique(clustering.labels_)
                                )
                        print()

                        # acts to terminate loop
                        return min_cluster_size

            def mcs_tb_submit(text):
                nonlocal range_tuple

                numerical_input = text.split('.save')[0].strip()
                new_range_tuple = tuple(map(
                    int, numerical_input.split('-')))

                # ensure code doesn't refire on
                # navigating away from text box
                if new_range_tuple != range_tuple:
                    range_tuple = new_range_tuple
                    min_cluster_size = cluster_and_plot()

                if '.save' in text:
                    if len(new_range_tuple) == 1:

                        print(
                            f'Saving current minimum cluster size: ' +
                            f'{new_range_tuple[0]}')
                        f = open(os.path.join(
                            clustering_dir, 'MCS.pkl'), 'wb')
                        pickle.dump(new_range_tuple[0], f)
                        f.close()

                        print('Saving current figure...')
                        plt.savefig(
                            os.path.join(
                                clustering_dir,
                                f'{self.embeddingAlgorithm}_'
                                f'{new_range_tuple[0]}.png'),
                            bbox_inches='tight', dpi=1000
                            )
                        plt.close('all')

                    else:
                        print('A single value must be entered to save an' +
                              'optimal minimum cluster size.')

            mcs_text_box = TextBox(
                ax_mcs_tb, "Minimum Cluster Size",
                range_tuple, label_pad=0.1)
            mcs_text_box.on_submit(mcs_tb_submit)

            label = mcs_text_box.ax.get_children()[1]
            label.set_position([0.5, 1.3])
            label.set_horizontalalignment('center')

            def sample_tb_submit(text):
                nonlocal sample

                new_sample = str(text)

                # ensure code doesn't refire on
                # navigating away from text box
                if new_sample != sample:
                    sample = new_sample
                    print(f'Setting sample-to-view to: {sample}')
                    print()

            sample_text_box = TextBox(
                ax_sample_tb, "Sample Name",
                f'{sample}', label_pad=0.1)
            sample_text_box.on_submit(sample_tb_submit)

            label = sample_text_box.ax.get_children()[1]
            label.set_position([0.5, 1.3])
            label.set_horizontalalignment('center')

            # kick things off with cluster_and_plot function
            cluster_and_plot()

            plt.show(block=True)

            # update new default MCS and reclassCutoffs
            # to the last ones entered so window reopens with
            # latest values (not original defaults)
            if len(range_tuple) == 1:
                new_default_range_tuple = range_tuple[0]
                new_default_sample = sample

            return (selector.ind, new_default_range_tuple, new_default_sample)

        #######################################################################
        # define while loop

        # define the point size for cells in the embedding
        point_size = 50000/len(data)

        # interact with plots to identify optimal min cluster size
        while not os.path.isfile(os.path.join(clustering_dir, 'MCS.pkl')):

            (selected_idx,
             default_range_tuple,
             default_sample) = do_plots()

            # if lassoed cells...
            if len(selected_idx) >= 1:

                # show highest expression channels
                data_copy = data.copy()

                # assign lassoed data a dummy cluster variable
                # and get highest expressed markers
                data_copy.loc[
                    data_copy.index.isin(selected_idx),
                    'cluster'] = 1000
                hi_markers = cluster_expression(
                    df=data_copy, markers=abx_channels,
                    cluster=1000, num_proteins=3,
                    standardize='within')

                print(f'Top three markers {hi_markers}')
                print()

                # superimpose centroids of lassoed noisy cells
                # colored by stage removed over channel images
                print('Opening Napari...')
                print()

                napari_warnings()
                if self.showAbChannels:

                    for e, ch in enumerate(reversed(abx_channels)):
                        channel_number = marker_channel_number(
                            markers, ch)

                        # read antibody image
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{default_sample}*.tif'):
                            img = single_channel_pyramid(
                                file_path,
                                channel=channel_number.item() - 1
                                )

                        # initialize Napari viewer with first channel
                        if e == 0:
                            viewer = napari.view_image(
                                img, rgb=False,
                                blending='additive',
                                colormap='green',
                                visible=False,
                                name=ch)
                        else:
                            viewer.add_image(
                                img, rgb=False,
                                blending='additive',
                                colormap='green',
                                visible=False,
                                name=ch)

                # color noisy data points by module used to redact them
                # len(data['cluster'].unique())
                cmap = categorical_cmap(
                    numUniqueSamples=len(data['Condition'].unique()),
                    numCatagories=10,
                    cmap='tab10',
                    continuous=False
                    )

                centroids = data[['Y_centroid', 'X_centroid']][
                        (data.index.isin(selected_idx))
                        & (data['Sample'] == default_sample)]

                viewer.add_points(
                    centroids, name='lassoed cells', visible=False,
                    face_color='lime', edge_width=0.0, size=4.0)

                # read segmentation outlines, add to Napari
                for file_path in glob.glob(
                  f'{self.inDir}/seg/{default_sample}*.ome.tif'):
                    seg = single_channel_pyramid(
                        file_path, channel=0)
                viewer.add_image(
                    seg, rgb=False,
                    blending='additive',
                    colormap='gray',
                    visible=False,
                    name='segmentation')

                # read last DNA, add to Napari
                last_dna_cycle = natsorted(
                    [i for i in data.columns if dna_moniker in i])[-1]
                channel_number = marker_channel_number(
                    markers, last_dna_cycle)
                for file_path in glob.glob(
                  f'{self.inDir}/tif/{default_sample}*.tif'):
                    dna_last = single_channel_pyramid(
                        file_path,
                        channel=channel_number.item() - 1)
                    viewer.add_image(
                        dna_last, rgb=False, blending='additive',
                        opacity=0.5, colormap='gray', visible=False,
                        name=f'{last_dna_cycle}: ' +
                        f'{default_sample}')

                # read first DNA, add to Napari
                for file_path in glob.glob(
                  f'{self.inDir}/tif/{default_sample}*.tif'):
                    dna_first = single_channel_pyramid(
                        file_path, channel=0)
                viewer.add_image(
                    dna_first, rgb=False, blending='additive',
                    opacity=0.5, colormap='gray', visible=True,
                    name=f'{dna1}: ' +
                    f'{default_sample}')

                napari.run()

        #######################################################################
        # apply final MCS and return data from clustering module

        f = open(os.path.join(clustering_dir, 'MCS.pkl'), 'rb')
        min_cluster_size = pickle.load(f)

        clustering = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=None,
            metric='euclidean', alpha=1.0, p=None,
            algorithm='best', leaf_size=40,
            memory=Memory(
                location=None),
            approx_min_span_tree=True,
            gen_min_span_tree=False,
            core_dist_n_jobs=-1,
            cluster_selection_method='eom',
            allow_single_cluster=False,
            prediction_data=False,
            match_reference_implementation=False).fit(
                data[['emb1', 'emb2']]
                )
        data['cluster'] = clustering.labels_

        print()
        print(
            f'FINAL min_cluster_size = {min_cluster_size}',
            np.unique(clustering.labels_)
            )

        print()
        print()
        return data

    @module
    def clustermap(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude)

        # create clustering directory if it hasn't already
        clustering_dir = os.path.join(self.outDir, 'clustering/final')
        if not os.path.exists(clustering_dir):
            os.makedirs(clustering_dir)

        # drop antibody channel exclusions for clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering]

        # drop unclustered cells before plotting clustermap
        clustermap_input = data[data['cluster'] != -1]

        # compute mean antibody signals for clusters
        clustermap_input = (
            clustermap_input[abx_channels + ['cluster']]
            .groupby('cluster').mean())

        sns.set(font_scale=0.8)
        for name, axis in zip(['channels', 'clusters'], [0, 1]):

            g = sns.clustermap(
                clustermap_input, cmap='viridis', standard_scale=axis,
                square=False, yticklabels=1, linewidth=0.1, cbar=True)

            g.fig.suptitle(f'Normalized across {name}', y=0.995, fontsize=10)
            g.fig.set_size_inches(6.0, 6.0)

            plt.savefig(os.path.join(
                    clustering_dir, f'clustermap_norm_{name}.pdf'),
                    bbox_inches='tight')

        plt.show(block=True)

        print()
        print()
        return data

    @module
    def setContrast(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude)

        # create contrast directory if it hasn't already
        contrast_dir = os.path.join(self.outDir, 'contrast')
        if not os.path.exists(contrast_dir):
            os.makedirs(contrast_dir)

        napari_warnings()

        # loop over antibody channels and add them to Napari viewer
        for e, ch in enumerate(reversed(abx_channels)):

            channel_number = marker_channel_number(markers, ch)

            # read antibody image
            for file_path in glob.glob(
              f'{self.inDir}/tif/{self.viewSample}*.tif'):
                img = single_channel_pyramid(
                    file_path, channel=channel_number.item() - 1)

            # initialize Napari viewer with first channel
            if e == 0:
                viewer = napari.view_image(
                    img, rgb=False, blending='additive',
                    colormap='green', visible=False, name=ch)

            else:
                viewer.add_image(
                    img, rgb=False, blending='additive',
                    colormap='green', visible=False, name=ch)

        # read DNA1 channel
        for file_path in glob.glob(
          f'{self.inDir}/tif/{self.viewSample}*.tif'):
            dna = single_channel_pyramid(file_path, channel=0)

        viewer.add_image(
            dna, rgb=False, blending='additive', colormap='gray',
            name=f'{dna1}: {self.viewSample}')

        # apply previously defined contrast limits if they exist
        if os.path.exists(os.path.join(contrast_dir, 'contrast_limits.yml')):

            print()
            print('Applying existing channel contrast settings.')

            contrast_limits = yaml.safe_load(
                open(f'{contrast_dir}/contrast_limits.yml'))

            viewer.layers[f'{dna1}: {self.viewSample}'].contrast_limits = (
                contrast_limits[dna1][0], contrast_limits[dna1][1])

            for ch in reversed(abx_channels):
                viewer.layers[ch].contrast_limits = (
                    contrast_limits[ch][0], contrast_limits[ch][1])

            napari.run()

        else:
            print()
            print('Channel contrast settings have not been defined.')

            napari.run()

            # create channel settings configuration file and
            # update with chosen constrast limits
            contrast_limits = {}
            for ch in [dna1] + abx_channels:
                if ch == dna1:
                    contrast_limits[ch] = (
                        viewer.layers[f'{dna1}: {self.viewSample}']
                        .contrast_limits)
                else:
                    contrast_limits[ch] = viewer.layers[ch].contrast_limits

            with open(f'{contrast_dir}/contrast_limits.yml', 'w') as file:
                yaml.dump(contrast_limits, file)

        print()
        print()
        return data

    @module
    def curateThumbnails(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude)

        # drop antibody channel exclusions for clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering]

        # create thumbnails directory if it hasn't already
        thumbnails_dir = os.path.join(
            self.outDir, 'clustering/final/thumbnails')
        if not os.path.exists(thumbnails_dir):
            os.mkdir(thumbnails_dir)

        # create zarr directory if it hasn't already
        zarr_dir = os.path.join(
            self.outDir, 'clustering/final/thumbnails/zarrs')
        if not os.path.exists(zarr_dir):
            os.mkdir(zarr_dir)

        # create list of tifs sorted from largest to smallest
        # (allows for early test of memory limitation)
        tifs = os.listdir(f'{self.inDir}/tif/')
        tifs.sort(
            key=lambda f: os.stat(os.path.join(
                f'{self.inDir}/tif/', f)).st_size, reverse=True)

        # remove tifs of samples censored from the analysis
        tifs = [i for i in tifs if i.split('.')[0] in data['Sample'].unique()]

        # get image contrast settings
        contrast_dir = os.path.join(self.outDir, 'contrast')
        if os.path.exists(f'{contrast_dir}/contrast_limits.yml'):
            contrast_limits = yaml.safe_load(
                open(f'{contrast_dir}/contrast_limits.yml'))

        # drop unclustered cells from data
        data = data[data['cluster'] != -1]

        #######################################################################

        # read the indices of clusters that have already been run
        if os.path.exists(os.path.join(
          thumbnails_dir, 'completed_clusters.pkl')):
            f = open(os.path.join(
                thumbnails_dir, 'completed_clusters.pkl'), 'rb')
            completed_clusters = pickle.load(f)
            total_clusters = set(data['cluster'].unique())
            clusters_to_run = natsorted(
                total_clusters.difference(set(completed_clusters)))
            print(f'Clusters to run: {len(clusters_to_run)}')
        else:
            # create a list of clusters to run
            completed_clusters = []
            clusters_to_run = natsorted(data['cluster'].unique())
            print(f'Clusters to run: {len(clusters_to_run)}')

        # read names of samples that have already been run
        if os.path.exists(
          os.path.join(thumbnails_dir, 'completed_zarrs.pkl')):
            f = open(os.path.join(
                thumbnails_dir, 'completed_zarrs.pkl'), 'rb')
            completed_zarrs = pickle.load(f)
            total_zarrs = [
                f'c{i}_{j}' for j in tifs
                for i in natsorted(data['cluster'].unique())]
            zarrs_to_run = [
                x for x in total_zarrs if x not in completed_zarrs]
            print(f'Samples to Zarr: {len(zarrs_to_run)}')
            print()
        else:
            # create a list of samples to run
            completed_zarrs = []
            total_zarrs = [
                f'c{i}_{j}' for j in tifs
                for i in natsorted(data['cluster'].unique())]
            zarrs_to_run = total_zarrs
            print(f'Samples to Zarr: {len(zarrs_to_run)}')
            print()

        for cluster in clusters_to_run:

            print(f'Cluster: {cluster}')

            # create dataframe to collect image thumbnails and their metadata
            long_table = pd.DataFrame()

            # get cycle 1 dna, plus top three expressed markers
            markers_to_show = cluster_expression(
                df=data, markers=abx_channels,
                cluster=cluster, num_proteins=3,
                standardize='within')

            markers_to_show = [
                i for i in [dna1] + markers_to_show[1].split(', ')]
            print(f'Markers to show: {markers_to_show}')
            print()

            # create marker color dict
            color_dict = {}
            for i, j in zip(
              markers_to_show,
              [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0),
               (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]):
                color_dict[i] = j

            for sample in tifs:

                print(f"Sample: {sample.split('.')[0]}")
                print(
                    f'Collecting cluster {cluster} ' +
                    f'centroids in {sample}')

                # crop out thumbnail images
                sample_cluster_subset = data[
                    (data['Sample'] == sample.split('.')[0])
                    & (data['cluster'] == cluster)]

                sample_cluster_subset.reset_index(
                    drop=True, inplace=True)

                if self.numThumbnails > len(sample_cluster_subset):

                    dif = self.numThumbnails - len(sample_cluster_subset)

                    extra_rows = pd.DataFrame(
                        data=0,
                        index=list(range(dif)),
                        columns=sample_cluster_subset.columns)
                    sample_cluster_subset = (
                        sample_cluster_subset.append(extra_rows))
                    sample_cluster_subset.reset_index(
                        drop=True, inplace=True)

                else:
                    sample_cluster_subset = (
                        sample_cluster_subset.sample(
                            n=self.numThumbnails, random_state=3))

                centroids = sample_cluster_subset[
                    ['X_centroid', 'Y_centroid']]

                clearRAM(print_usage=False)
                del sample_cluster_subset
                clearRAM(print_usage=False)

                # if a persistent image and seg zarrs for current sample
                # don't already exist, create them
                z1_path = os.path.join(
                    zarr_dir, f'c{cluster}_{sample}.zarr')

                z2_path = os.path.join(
                    zarr_dir, f'c{cluster}_{sample}_seg.zarr')

                if f'c{cluster}_{sample}' in zarrs_to_run:

                    print(
                        'Creating Zarr array for ' +
                        f'cluster {cluster} cells in {sample}')

                    # loop over markers to show to create multi-channel image
                    for e, marker in enumerate(markers_to_show):

                        # if DNA1 channel
                        if e == 0:

                            print(f'Reading {marker} image')

                            # read cycle1 dna channel, convert to float and rgb
                            file_path = os.path.join(
                                f'{self.inDir}/tif/', sample)

                            channel_number = 0
                            tiff = TiffFile(file_path, is_ome=False)
                            dna = zarr.open(
                                tiff.series[0].levels[0][
                                    channel_number].aszarr())

                            dna = img_as_float(dna)

                            dna = gray2rgb(dna)

                        else:
                            print(f'Overlaying {marker} image')

                            # read antibody image, convert to float and rgb
                            file_path = os.path.join(
                                f'{self.inDir}/tif/', sample)

                            channel_number = marker_channel_number(
                                markers, marker)
                            tiff = TiffFile(file_path, is_ome=False)
                            img = zarr.open(
                                tiff.series[0].levels[0][
                                    channel_number.item() - 1].aszarr())

                            img = img_as_float(img)

                            # apply image contrast settings
                            img -= (contrast_limits[marker][0]/65535)
                            img /= (
                                (contrast_limits[marker][1]/65535)
                                - (contrast_limits[marker][0]/65535))

                            img = np.clip(img, 0, 1)

                            img = gray2rgb(img)

                            img = (img * color_dict[marker])

                            # add antibody image to cycle1 dna images
                            dna += img

                            clearRAM(print_usage=False)
                            del img
                            clearRAM(print_usage=False)

                    print('Overlaying cluster centroids image')

                    # create blank array (zeros) for applying a centroids mask
                    centroid_img = np.zeros((dna.shape[0], dna.shape[1]))

                    # specify centroid size (in pixels)
                    centroid_dist = 1

                    # loop over centroids and add them to blank image
                    for example, centroid in enumerate(centroids.iterrows()):

                        ystart_centroid = int(
                            centroid[1]['Y_centroid'] - centroid_dist)
                        ystop_centroid = int(
                            centroid[1]['Y_centroid'] + centroid_dist)

                        xstart_centroid = int(
                            centroid[1]['X_centroid'] - centroid_dist)
                        xstop_centroid = int(
                            centroid[1]['X_centroid'] + centroid_dist)

                        centroid_img[
                            ystart_centroid:ystop_centroid,
                            xstart_centroid:xstop_centroid
                            ] = 1

                    # convert to rgb and colorize
                    centroid_img = gray2rgb(centroid_img)
                    centroid_img = (centroid_img * (1.0, 1.0, 1.0))

                    # add to overlay
                    dna += centroid_img

                    clearRAM(print_usage=False)
                    del centroid_img
                    clearRAM(print_usage=False)

                    print('Reading cell segmentation outlines image')

                    file_path = os.path.join(f'{self.inDir}/seg/', sample)

                    # check for .tif or .ome.tif file extension
                    # (TMA images are still .tif files in mcmicro)
                    if os.path.exists(file_path):
                        seg_img = imread(file_path, key=0)
                    else:
                        file_path = os.path.join(
                            f'{self.inDir}/seg/',
                            f"{sample.split('.tif')[0]}.ome.tif")
                        seg_img = imread(file_path, key=0)

                    seg_img = gray2rgb(seg_img)

                    seg_img = (seg_img * (1.0, 1.0, 1.0))

                    # assign multi-channel array to a persistent zarr array
                    print('Saving image overlay as Zarr array')

                    z1 = zarr.open(
                        z1_path, mode='w',
                        shape=(dna.shape[0], dna.shape[1], dna.shape[2]),
                        chunks=(2500, 2500), dtype='f8')
                    z1[:] = dna

                    # assign segmenation outlines to a separate zarr array
                    print(
                        'Saving cell segmentation outlines as ' +
                        'separate Zarr array')

                    z2 = zarr.open(
                        z2_path, mode='w',
                        shape=(dna.shape[0], dna.shape[1], dna.shape[2]),
                        chunks=(2500, 2500), dtype='f8')
                    z2[:] = seg_img

                    # update completed_zarrs list
                    completed_zarrs.append(f'c{cluster}_{sample}')
                    f = open(os.path.join(
                        thumbnails_dir, 'completed_zarrs.pkl'), 'wb')
                    pickle.dump(completed_zarrs, f)
                    f.close()

                else:
                    print(
                        'Reading Zarr array for ' +
                        f'cluster {cluster} cells in {sample}')
                    z1 = zarr.open(z1_path, mode='r')

                # crop thumbnails
                for example, centroid in enumerate(centroids.iterrows()):

                    if (
                        (centroid[1]['X_centroid'] == 0.0) &
                        (centroid[1]['Y_centroid'] == 0.0)
                      ):

                        blank_img = np.ones(
                            (self.squareWindowDimension,
                             self.squareWindowDimension))

                        long_table = long_table.append(
                            {'sample': sample.split('.')[0],
                             'example': int(example),
                             'image': blank_img},
                            ignore_index=True)

                    else:
                        # specify window x, y ranges
                        ystart_window = int(
                            centroid[1]['Y_centroid']
                            - self.squareWindowDimension)
                        ystop_window = int(
                            centroid[1]['Y_centroid']
                            + self.squareWindowDimension)

                        xstart_window = int(
                            centroid[1]['X_centroid']
                            - self.squareWindowDimension)
                        xstop_window = int(
                            centroid[1]['X_centroid']
                            + self.squareWindowDimension)

                        # for centroids falling within
                        # self.squareWindowDimension pixels of the edge of the
                        # dna image, ensure that the thumbnail image is not
                        # cropped using negative slicing values, as it will
                        # return an empty array and lead to a runtime error
                        # during plotting.
                        window_list = [
                            ystart_window, ystop_window,
                            xstart_window, xstop_window]
                        (ystart_window,
                         ystop_window,
                         xstart_window,
                         xstop_window) = [
                            0 if i < 0 else i for i in window_list]

                        # crop overlay image to window size
                        thumbnail = z1[
                            ystart_window:ystop_window,
                            xstart_window:xstop_window]

                        # add egmentation outlines to thumbnail images
                        if self.segOutlines:
                            z2 = zarr.open(z2_path, mode='r')

                            seg_thumbnail = z2[
                                ystart_window:ystop_window,
                                xstart_window:xstop_window]

                            thumbnail += seg_thumbnail

                        long_table = long_table.append(
                            {'sample': sample.split('.')[0],
                             'example': int(example),
                             'image': thumbnail},
                            ignore_index=True)
                print()

            # known pandas foible: integers are stored as floats
            long_table['example'] = [int(i) for i in long_table['example']]

            # natsort long_table by 'sample' column
            long_table['sample'] = pd.Categorical(
                long_table['sample'], ordered=True,
                categories=natsorted(long_table['sample'].unique()))
            long_table.sort_values('sample', inplace=True)
            long_table['sample'] = long_table['sample'].astype(str)

            # plot cluster facet grid
            fig, ax = plt.subplots()

            g = sns.FacetGrid(
                long_table, row='sample', col='example',
                sharex=False, sharey=False,
                gridspec_kws={'hspace': 0.1, 'wspace': 0.05})

            g.map(
                lambda image, **kwargs: (
                    plt.imshow(np.clip(image.values[0], 0, 1)),
                    plt.grid(False)), 'image')
            # image clipping prevents matplotlib warning

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
                loc='upper left')

            plt.savefig(
                os.path.join(
                    thumbnails_dir,
                    'cluster' + str(cluster) + '_thumbnails.pdf'),
                bbox_inches='tight')
            plt.close('all')

            completed_clusters.append(cluster)

            f = open(os.path.join(
                thumbnails_dir, 'completed_clusters.pkl'), 'wb')
            pickle.dump(completed_clusters, f)
            f.close()
            print()

        print()
        print()
        return data

    @module
    def frequencyStats(data, self, args):

        # prepare input data for computing statistics
        stats_input = data[['Sample', 'Replicate', 'cluster']][
            data['cluster'] >= 0]

        # loop over comma-delimited binary declarations
        for i in range(len(list(self.sampleStatuses.values())[0].split(', '))):

            # get unique declaration categories (should be 2 per test)
            comparison = set(
                [j.split(', ')[i] for j in self.sampleStatuses.values()
                 if '-UNK' not in j.split(', ')[i]])

            if len(comparison) > 1:

                # assign test and control groups
                test = [
                    i for i in comparison if i not in self.controlGroups][0]
                control = [
                    i for i in comparison if i in self.controlGroups][0]

                # create frequency stats directory if it hasn't already
                frequency_dir = os.path.join(
                    self.outDir, 'clustering/final/frequency_stats',
                    f"{test}_v_{control}")
                if not os.path.exists(frequency_dir):
                    os.makedirs(frequency_dir)

                # create single-column dataFrame containing all sample names
                # to pad counts tables with zeros if a celltype
                # is not in a tissue
                pad = pd.DataFrame(
                    natsorted(stats_input['Sample'].unique())).rename(
                        columns={0: 'Sample'})

                cluster_list = []
                ratio_list = []
                dif_list = []
                pval_list = []

                # intialize a dataframe to collect catplot data
                catplot_input = pd.DataFrame()

                # loop over clusters
                for cluster, group in stats_input.groupby('cluster'):

                    print(
                        f'Calculating log2({test}/{control})'
                        f' of mean cell density for cluster {str(cluster)}.')

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

                    # get sample file names (i.e. sampleMetadata keys)
                    # from config.yml based on "Sample" column
                    # (first elements of sampleMetadata vals)
                    def get_key(val):
                        for key, value in self.sampleNames.items():
                            if val == value:
                                return key

                        return "key doesn't exist"
                    file_names = [get_key(i) for i in group['Sample']]

                    # add binary declarations column to group data
                    group['status'] = [
                        self.sampleStatuses[j].split(', ')[i]
                        for j in file_names]

                    # add replicates column to group data
                    group['Replicate'] = [
                        self.sampleReplicates[i]
                        for i in file_names]

                    group['cluster'] = cluster

                    # drop samples for which a declaration cannot be made
                    group = group[~group['status'].str.contains('-UNK')]

                    group.reset_index(drop=True, inplace=True)

                    # get denominator cell count for each sample
                    if self.denominatorCluster is None:
                        group['tissue_count'] = [
                            len(stats_input[stats_input['Sample'] == i]) for
                            i in group['Sample']]
                    else:
                        group['tissue_count'] = [
                            len(stats_input[(stats_input['Sample'] == i) &
                                (stats_input['cluster'] ==
                                 self.denominatorCluster)])
                            for i in group['Sample']]

                    # compute density of cells per sample
                    group['density'] = group['count']/group['tissue_count']

                    # append group data to catplot_input
                    catplot_input = catplot_input.append(group)

                    # isolate test and control group values
                    cnd1_values = group['density'][group['status'] == test]
                    cnd2_values = group['density'][group['status'] == control]

                    # perform Welch's t-test (equal_var=False)
                    stat, pval = ttest_ind(
                        cnd1_values, cnd2_values,
                        axis=0, equal_var=False, nan_policy='propagate')

                    # round resulting values
                    stat = round(stat, 3)
                    pval = round(pval, 3)

                    # compute mean of test and control group values
                    cnd1_mean = np.mean(cnd1_values)
                    cnd2_mean = np.mean(cnd2_values)

                    # compute mean ratio
                    ratio = np.log2(
                        (cnd1_mean + 0.000001)/(cnd2_mean + 0.000001))

                    # compute mean difference
                    dif = cnd1_mean-cnd2_mean

                    cluster_list.append(cluster)
                    ratio_list.append(ratio)
                    dif_list.append(dif)
                    pval_list.append(pval)

                # create stats dataframe
                statistics = pd.DataFrame(
                    list(zip(cluster_list, ratio_list, dif_list, pval_list)),
                    columns=['cluster', 'ratio', 'dif', 'pval']).sort_values(
                        by='cluster')

                # save total stats table
                statistics.to_csv(
                    os.path.join(
                        frequency_dir, 'stats_total.csv'), index=False)

                # compute FDR p-val corrections
                # (uses statsmodels.stats.multitest implementation)
                rejected, p_adjust = fdrcorrection(
                    statistics['pval'].tolist(), alpha=0.05,
                    method='indep', is_sorted=False)

                statistics['qval'] = p_adjust

                if self.FDRCorrection:
                    stat = 'qval'
                else:
                    stat = 'pval'

                # isolate statistically significant stat values
                significant = statistics[
                    statistics[stat] <= 0.05].sort_values(by=stat)

                # save significant stats table
                significant.to_csv(
                    os.path.join(
                        frequency_dir, 'stats_sig.csv'), index=False)

                # plot
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
                    zip(natsorted(catplot_input['Sample'].unique()),
                        cmap.colors))

                catplot_input.sort_values(
                    by=['cluster', 'status', 'density'],
                    ascending=[True, False, True], inplace=True)

                catplot_input['cluster'] = (
                    catplot_input['cluster'].astype(str) + f'; {stat} = ' +
                    catplot_input[stat].astype(str)
                    )

                catplot_input['cluster'] = [
                    i.split(f'; {stat} = ns')[0]
                    for i in catplot_input['cluster']]

                sns.set(font_scale=0.4)
                g = sns.catplot(
                    x='status', y='density',
                    hue=catplot_input['Sample'], col='cluster', col_wrap=6,
                    data=catplot_input, kind='bar', palette=sample_color_dict,
                    height=2, aspect=0.8, sharex=True, sharey=False,
                    edgecolor='k', linewidth=0.1, legend=False)

                g.set(ylim=(0.0, None))

                file_names = [
                    get_key(i) for i in
                    natsorted(catplot_input['Sample'].unique())]

                sample_conds = [
                    self.sampleConditions[i]
                    for i in file_names]

                sample_abbrs = [
                    self.sampleConditionAbbrs[i]
                    for i in file_names]

                cond_abbr = [
                    f'{i}-{j}' for i, j in zip(sample_conds, sample_abbrs)]

                handles_dict = dict(zip(
                    natsorted(catplot_input['Sample'].unique()), cond_abbr))

                legend_handles = []
                for k, v in handles_dict.items():
                    legend_handles.append(
                        Line2D([0], [0], marker='o', color='none',
                               label=v, markerfacecolor=sample_color_dict[k],
                               markeredgecolor='k', markeredgewidth=0.2,
                               markersize=5.0))

                plt.legend(
                    handles=legend_handles,
                    prop={'size': 5.0},
                    bbox_to_anchor=[1.03, 1.0])

                plt.savefig(
                    os.path.join(frequency_dir, 'catplot.pdf'),
                    bbox_inches='tight')
                plt.close('all')

                print()

            else:
                print(
                    'Only one binary declaration ' +
                    f'class represented for {list(comparison)[0]}. ' +
                    'No statistics will be computed.')
                print()

        print()
        print()
        return data
