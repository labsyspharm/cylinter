import logging
import functools

import sys
import os
import re
import glob
import yaml
import csv as csv_module
import math
import pickle
from ast import literal_eval

import gc
import hdbscan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from magicgui import magicgui
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from qtpy.QtCore import QTimer

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
from matplotlib import animation

from PIL import Image, ImageDraw

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
# from numcodecs import Blosc
from datetime import datetime
from joblib import Memory
from scipy.stats import ttest_ind

from subprocess import run

from .utils import (
    dataset_files, log_banner, log_multiline, reorganize_dfcolumns,
    SelectFromCollection, read_dataframe, save_dataframe, read_markers,
    marker_channel_number, categorical_cmap, cluster_expression, clearRAM,
    single_channel_pyramid, matplotlib_warnings, napari_warnings,
    fdrcorrection, open_file, dict_to_csv, csv_to_dict
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
                 samplesForROISelection=None,

                 # crossCycleCorrelation -
                 yAxisGating=None,

                 # pruneOutliers -
                 hexbins=None,
                 hexbinGridSize=None,

                 # metaQC -
                 metaQC=None,
                 default_mcs=200,
                 default_reclass_tuple='0.75, 0.75',
                 embeddingAlgorithmQC=None,
                 channelExclusionsClusteringQC=None,
                 samplesToRemoveClusteringQC=None,
                 fracForEmbeddingQC=None,
                 dimensionEmbeddingQC=None,
                 topMarkersQC=None,
                 metricQC=None,
                 perplexityQC=None,
                 earlyExaggerationQC=None,
                 learningRateTSNEQC=None,

                 randomStateQC=None,
                 nNeighborsQC=None,
                 learningRateUMAPQC=None,
                 minDistQC=None,
                 repulsionStrengthQC=None,

                 # PCA module —
                 channelExclusionsPCA=None,
                 samplesToRemovePCA=None,
                 dimensionPCA=None,
                 pointSize=None,
                 labelPoints=None,
                 distanceCutoff=None,
                 conditionsToSilhouette=None,

                 # clustering module —
                 embeddingAlgorithm=None,
                 channelExclusionsClustering=None,
                 samplesToRemoveClustering=None,
                 normalizeTissueCounts=None,
                 fracForEmbedding=None,
                 dimensionEmbedding=None,
                 topMarkers=None,
                 colormapChannel=None,
                 perplexity=None,
                 earlyExaggeration=None,
                 learningRateTSNE=None,
                 metric=None,
                 randomStateTSNE=None,
                 nNeighbors=None,
                 learningRateUMAP=None,
                 minDist=None,
                 repulsionStrength=None,
                 randomStateUMAP=None,

                 # frequencyStats —
                 controlGroups=None,
                 denominatorCluster=None,
                 FDRCorrection=None,

                 # curateThumbnails —
                 numThumbnails=None,
                 topMarkersThumbnails=None,
                 windowSize=None,
                 segOutlines=None,
                 ):

        # assert(SOMETHING)  # placeholder

        self.inDir = inDir
        self.outDir = outDir
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
        self.samplesForROISelection = samplesForROISelection

        self.yAxisGating = yAxisGating

        self.hexbins = hexbins
        self.hexbinGridSize = hexbinGridSize

        self.metaQC = metaQC
        self.default_mcsQC = default_mcs
        self.default_reclass_tuple = default_reclass_tuple
        self.embeddingAlgorithmQC = embeddingAlgorithmQC
        self.channelExclusionsClusteringQC = channelExclusionsClusteringQC
        self.samplesToRemoveClusteringQC = samplesToRemoveClusteringQC
        self.fracForEmbeddingQC = fracForEmbeddingQC
        self.dimensionEmbeddingQC = dimensionEmbeddingQC
        self.topMarkersQC = topMarkersQC
        self.metricQC = metricQC
        self.perplexityQC = perplexityQC
        self.earlyExaggerationQC = earlyExaggerationQC
        self.learningRateTSNEQC = learningRateTSNEQC
        self.randomStateQC = randomStateQC
        self.nNeighborsQC = nNeighborsQC
        self.learningRateUMAPQC = learningRateUMAPQC
        self.minDistQC = minDistQC
        self.repulsionStrengthQC = repulsionStrengthQC

        self.channelExclusionsPCA = channelExclusionsPCA
        self.samplesToRemovePCA = samplesToRemovePCA
        self.dimensionPCA = dimensionPCA
        self.pointSize = pointSize
        self.labelPoints = labelPoints
        self.distanceCutoff = distanceCutoff
        self.conditionsToSilhouette = conditionsToSilhouette

        self.embeddingAlgorithm = embeddingAlgorithm
        self.channelExclusionsClustering = channelExclusionsClustering
        self.samplesToRemoveClustering = samplesToRemoveClustering
        self.normalizeTissueCounts = normalizeTissueCounts
        self.fracForEmbedding = fracForEmbedding
        self.dimensionEmbedding = dimensionEmbedding
        self.topMarkers = topMarkers
        self.colormapChannel = colormapChannel
        self.perplexity = perplexity
        self.earlyExaggeration = earlyExaggeration
        self.learningRateTSNE = learningRateTSNE
        self.metric = metric
        self.randomStateTSNE = randomStateTSNE
        self.nNeighbors = nNeighbors
        self.learningRateUMAP = learningRateUMAP
        self.minDist = minDist
        self.repulsionStrength = repulsionStrength
        self.randomStateUMAP = randomStateUMAP

        self.controlGroups = controlGroups
        self.denominatorCluster = denominatorCluster
        self.FDRCorrection = FDRCorrection

        self.numThumbnails = numThumbnails
        self.topMarkersThumbnails = topMarkersThumbnails
        self.windowSize = windowSize
        self.segOutlines = segOutlines

    @module
    def aggregateData(data, self, args):

        print()
        files = natsorted(dataset_files(f'{self.inDir}/csv'))

        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=None
            )

        df_list = []
        raw_sample_names_dict = {}
        channel_setlist = []
        for file in files:
            if not file.startswith('.'):

                file_name = file.split('.csv')[0]

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
#                     mask_dict = {
#                         'Hoechst0': 'nucleiRingMask',
#                         'Hoechst1': 'nucleiRingMask',
#                         'Hoechst2': 'nucleiRingMask',
#                         'anti_CD3': 'cytoRingMask',
#                         'anti_CD45RO': 'cytoRingMask',
#                         'Hoechst3': 'nucleiRingMask',
#                         'Keratin_570': 'cellRingMask',
#                         'aSMA_660': 'cellRingMask',
#                         'Hoechst4': 'nucleiRingMask',
#                         'CD4_488': 'cytoRingMask',
#                         'CD45_PE': 'cytoRingMask',
#                         'PD1_647': 'cytoRingMask',
#                         'Hoechst5': 'nucleiRingMask',
#                         'CD20_488': 'cytoRingMask',
#                         'CD68_555': 'cellRingMask',
#                         'CD8a_660': 'cytoRingMask',
#                         'Hoechst6': 'nucleiRingMask',
#                         'CD163_488': 'cellRingMask',
#                         'FOXP3_570': 'nucleiRingMask',
#                         'PDL1_647': 'cytoRingMask',
#                         'Hoechst7': 'nucleiRingMask',
#                         'Ecad_488': 'cellRingMask',
#                         'Vimentin_555': 'cellRingMask',
#                         'CDX2_647': 'cellRingMask',
#                         'Hoechst8': 'nucleiRingMask',
#                         'LaminABC_488': 'nucleiRingMask',
#                         'Desmin_555': 'cellRingMask',
#                         'CD31_647': 'nucleiRingMask',
#                         'Hoechst9': 'nucleiRingMask',
#                         'PCNA_488': 'nucleiRingMask',
#                         'CollagenIV_647': 'cellRingMask'}

#                     mask_object_cols = (
#                         ['CellID', 'X_centroid', 'Y_centroid', 'Area',
#                          'MajorAxisLength', 'MinorAxisLength',
#                          'Eccentricity', 'Solidity', 'Extent',
#                          'Orientation'] +
#                         [f'{i}_{mask_dict[i]}' for i
#                          in markers['marker_name']])

                    # select boilerplate columns
                    cols = (
                        ['CellID', 'X_centroid', 'Y_centroid', 'Area',
                         'MajorAxisLength', 'MinorAxisLength',
                         'Eccentricity', 'Solidity', 'Extent',
                         'Orientation'] +
                        [i for i in markers['marker_name'] if i in csv.columns]
                         )

                    csv = csv[cols]

                    # add sample column
                    csv['Sample'] = sample_name

                    # add condition column
                    csv['Condition'] = self.sampleConditionAbbrs[file_name]

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
            print()
            print(
                f'Features {markers_to_drop} were not probed in all' +
                ' samples and will be dropped from the analysis.'
                )
        data = data[channels_set].copy()

        # sort by Sample and CellID to be tidy
        data.sort_values(by=['Sample', 'CellID'], inplace=True)

        # assign global index
        data.reset_index(drop=True, inplace=True)

        # ensure MCMICRO-generated columns come first and
        # are in the same order as csv input
        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def selectROIs(data, self, args):

        print()
        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # create ROIs directory if it doesn't already exist
        roi_dir = os.path.join(self.outDir, 'ROIs')
        if not os.path.exists(roi_dir):
            os.makedirs(roi_dir)

        # load polygon dictionary if it exists
        if os.path.exists(os.path.join(roi_dir, 'polygon_dict.pkl')):
            f = open(os.path.join(roi_dir, 'polygon_dict.pkl'), 'rb')
            polygon_dict = pickle.load(f)
        else:
            # if polygon dictionary does not exist, create it
            polygon_dict = {}

        napari_warnings()

        # loop over samples
        for sample_name in self.samplesForROISelection:

            if sample_name in data['Sample'].unique():

                print(f'Working on sample {sample_name}')

                try:
                    # make a list of exisiting polygon vertices
                    polygon_dict[sample_name]

                    shapes = [polygon_dict[sample_name][i][0] for
                              i in range(0, len(polygon_dict[sample_name]))]
                    polygons = [polygon_dict[sample_name][i][1] for
                                i in range(0, len(polygon_dict[sample_name]))]

                except KeyError:
                    # create empty list to store polygon vertices
                    polygons = []

                if self.showAbChannels:
                    for e, ch in enumerate(reversed(abx_channels)):
                        channel_number = marker_channel_number(markers, ch)

                        # read antibody image
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{sample_name}.*tif'):

                            img = single_channel_pyramid(
                                file_path, channel=channel_number.item() - 1
                                )

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
                file_path = f'{self.inDir}/seg/{sample_name}.*tif'
                seg = single_channel_pyramid(
                    glob.glob(file_path)[0], channel=0)

                viewer.add_image(
                    seg, rgb=False, blending='additive',
                    opacity=0.5, colormap='gray', visible=False,
                    name='segmentation')

                # read DNA1 channel
                for file_path in glob.glob(
                  f'{self.inDir}/tif/{sample_name}.*tif'):
                    dna = single_channel_pyramid(file_path, channel=0)

                viewer.add_image(
                    dna, rgb=False, blending='additive',
                    colormap='gray', visible=True,
                    name=f'{dna1}: {sample_name}')

                if polygons:
                    selection_layer = viewer.add_shapes(
                        data=polygons,
                        shape_type=shapes,
                        ndim=2,
                        face_color=[1.0, 1.0, 1.0, 0.2],
                        edge_color=[0.0, 0.66, 1.0, 1.0],
                        edge_width=10.0,
                        name='ROI(s)')
                else:
                    selection_layer = viewer.add_shapes(
                        data=None,
                        shape_type='polygon',
                        ndim=2,
                        face_color=[1.0, 1.0, 1.0, 0.2],
                        edge_color=[0.0, 0.66, 1.0, 1.0],
                        edge_width=10.0,
                        name='ROI(s)')

                # clear polygon vertices from polygons list for re-entry below
                updated_polygons = []

                napari.run()

                # store lists vertices per sample as a dictionary
                for shape_type, roi in zip(
                  selection_layer.shape_type, selection_layer.data):
                    updated_polygons.append((shape_type, roi))

                polygon_dict[sample_name] = updated_polygons

                f = open(os.path.join(roi_dir, 'polygon_dict.pkl'), 'wb')
                pickle.dump(polygon_dict, f)
                f.close()

        print()

        # create txt files per sample with cell IDs to drop
        idxs_to_drop = {}
        for sample_name in natsorted(data['Sample'].unique()):

            # test if sample key in polygon_dict
            try:
                polygon_dict[sample_name]

                # if polygon list is not empty
                if polygon_dict[sample_name]:

                    print(
                        'Generating polygon mask(s) ' +
                        f'for sample: {sample_name}')

                    sample_data = data[['X_centroid', 'Y_centroid', 'CellID']][
                        data['Sample'] == sample_name].astype(int)

                    sample_data['tuple'] = list(
                        zip(sample_data['X_centroid'],
                            sample_data['Y_centroid'])
                        )

                    for file_path in glob.glob(
                      f'{self.inDir}/tif/{sample_name}.*tif'):
                        dna = imread(file_path, key=0)

                    columns, rows = np.meshgrid(
                        np.arange(dna.shape[1]),
                        np.arange(dna.shape[0])
                        )
                    columns, rows = columns.flatten(), rows.flatten()
                    pixel_coords = np.vstack((rows, columns)).T

                    # create pillow image to convert into boolean mask
                    img = Image.new('L', (dna.shape[1], dna.shape[0]))

                    polygons = []
                    for shape_type, verts in polygon_dict[sample_name]:

                        selection_verts = np.round(verts).astype(int)

                        if shape_type == 'ellipse':
                            col_min = selection_verts[0][1]
                            col_max = selection_verts[1][1]
                            row_min = selection_verts[0][0]
                            row_max = selection_verts[2][0]

                            ellipse = [(col_min, row_min), (col_max, row_max)]

                            # update pillow image with ellipse
                            ImageDraw.Draw(img).ellipse(
                                ellipse, outline=1, fill=1
                                )

                        else:
                            polygon = list(tuple(
                                zip(selection_verts[:, 1],
                                    selection_verts[:, 0])
                                ))

                            # update pillow image with polygon
                            ImageDraw.Draw(img).polygon(
                                polygon, outline=1, fill=1
                                )

                    # convert pillow image into boolean numpy array
                    mask = np.array(img, dtype=bool)

                    # use numpy fancy indexing to get centroids
                    # where boolean mask is True
                    xs, ys = zip(*sample_data['tuple'])

                    inter = mask[ys, xs]

                    # update sample_data with boolean calls per cell
                    sample_data['inter'] = inter

                    if self.delintMode is True:
                        idxs_to_drop[sample_name] = list(
                            sample_data['CellID'][sample_data['inter']]
                            )

                    else:
                        idxs_to_drop[sample_name] = list(
                            sample_data['CellID'][~sample_data['inter']]
                            )

                else:
                    print(f'No ROIs selected for sample: {sample_name}')
                    idxs_to_drop[sample_name] = []

            except KeyError:
                print(f'No ROIs selected for sample: {sample_name}')
                idxs_to_drop[sample_name] = []

        print()

        # drop cells from samples
        for sample_name, cell_ids in idxs_to_drop.items():
            if cell_ids:
                print(f'Dropping cells from sample: {sample_name}')
                global_idxs_to_drop = data[
                    (data['Sample'] == sample_name)
                    & (data['CellID'].isin(set(cell_ids)))].index
                data.drop(global_idxs_to_drop, inplace=True)
            else:
                pass
        print()

        # save images of tissue with selected data points
        image_dir = os.path.join(roi_dir, 'images')
        if not os.path.exists(image_dir):
            os.mkdir(image_dir)

        for sample_name in natsorted(data['Sample'].unique()):

            print(f'Plotting ROI selections for sample: {sample_name}')

            for file_path in glob.glob(
              f'{self.inDir}/tif/{sample_name}.*tif'):
                dna = imread(file_path, key=0)

            fig, ax = plt.subplots()
            ax.imshow(dna, cmap='gray')
            ax.grid(False)
            coords = data[['X_centroid', 'Y_centroid', 'Area']][
                data['Sample'] == sample_name]
            sp = ax.scatter(
                coords['X_centroid'], coords['Y_centroid'],
                s=0.35, lw=0.0,
                c=coords['Area'], cmap='viridis')
            plt.title(
                f'Sample {sample_name}. ' +
                'Selected cells colored by segmentation area')
            plt.colorbar(sp)
            plt.savefig(
                os.path.join(
                    image_dir, f'{sample_name}.png'), dpi=1000)
            plt.close('all')

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

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
            data=data
            )

        # create intensity directory if it doesn't already exist
        intensity_dir = os.path.join(self.outDir, 'intensity')
        if not os.path.exists(intensity_dir):
            os.makedirs(intensity_dir)

        # set histogram bin size
        num_bins = 125

        # set histogram type
        histtype = 'stepfilled'

        # pick up where samples-loop left off
        if os.path.exists(
          os.path.join(intensity_dir, 'idxs_to_drop.csv')):

            # read stored dict
            idxs_to_drop = csv_to_dict(
                os.path.join(intensity_dir, 'idxs_to_drop.csv'))

            # dictionary values from strings to lists
            idxs_to_drop = {
                k: literal_eval(v) for k, v in idxs_to_drop.items()}

            # get names of samples already run
            previously_run_samples = idxs_to_drop.keys()

            # get names of samples remaining
            samples_remaining = (
                len(data['Sample'].unique())
                - len(previously_run_samples))

            print()
            print(f'Samples to threshold: {samples_remaining}')

            # drop samples previously run
            df = data[~data['Sample'].isin(previously_run_samples)].copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique()))
            df.sort_values('Sample', inplace=True)

            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        else:
            # get names of samples remaining (total samples in this case)
            samples_remaining = len(data['Sample'].unique())

            print()
            print(f'Samples to threshold: {samples_remaining}')

            # initialize dictionary of sample indices to drop
            idxs_to_drop = {}

            # drop samples previously run (total samples in this case)
            df = data.copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique()))
            df.sort_values('Sample', inplace=True)
            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        # loop over all (or remaining) samples
        for name, group in df.groupby('Sample', sort=False):

            print()
            print(f'Sample: {name}')

            # read segmentation outlines
            file_path = f'{self.inDir}/seg/{name}.*tif'

            # create image pyramid
            seg = single_channel_pyramid(
                glob.glob(file_path)[0], channel=0)

            # add DNA image to viewer
            viewer = napari.view_image(
                seg, rgb=False, visible=False, colormap='gray',
                opacity=0.5, name='segmentation'
                )

            # read DNA1 channel
            for file_path in glob.glob(
              f'{self.inDir}/tif/{name}.*tif'):
                # create image pyramid
                dna = single_channel_pyramid(file_path, channel=0)

            # add segmentation image to viewer
            viewer.add_image(
                dna, rgb=False, blending='additive',
                name=f'{dna1}: {name}'
                )

            # generate Qt widget to dock in Napari viewer
            hist_widget = QtWidgets.QWidget()

            # generate a blank figure canvas
            canvas = FigureCanvas(Figure())

            # construct vertical box layout object
            # to line up widgets vertically
            hist_layout = QtWidgets.QVBoxLayout(hist_widget)

            # add navigation tool bar and figure canvas to widget
            hist_layout.addWidget(NavigationToolbar(canvas, hist_widget))
            hist_layout.addWidget(canvas)

            # set plot style
            sns.set_style('whitegrid')

            # get figure object from canvas
            fig = canvas.figure

            # adjust plot on canvas to taste
            fig.subplots_adjust(left=0.25, bottom=0.25)

            # set plot title
            fig.suptitle(f'Sample={name} mean DNA intensity', size=10)

            # get axis object from canvas
            ax = canvas.figure.subplots()

            # plot histogram on canvas axis
            n, bins, patches = ax.hist(
                group[dna1], bins=num_bins,
                density=False, color='grey', ec='none',
                alpha=0.75, histtype=histtype,
                range=None, label='before')

            # set y-axis label
            ax.set_ylabel('count')

            # add sliders to plot
            axcolor = 'lightgoldenrodyellow'
            axLowerCutoff = fig.add_axes(
                [0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
            axUpperCutoff = fig.add_axes(
                [0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

            # specify data range
            rnge = [bins.min(), bins.max()]

            # add slider functionality
            sLower = Slider(
                axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
                valinit=0.00, valstep=(rnge[1]/100000))
            sLower.label.set_fontsize(11)
            sLower.label.set_color('b')
            sUpper = Slider(
                axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
                valinit=0.00, valstep=(rnge[1]/100000))
            sUpper.label.set_fontsize(11)
            sUpper.label.set_color('r')

            # specify function for updating sliders
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

            # update sliders when moved
            sLower.on_changed(update)
            sUpper.on_changed(update)

            # add button to show selected centroids in Napari viewer
            button_ax = fig.add_axes([0.65, 0.025, 0.25, 0.06])
            button = Button(
                button_ax, 'Plot Points', color=axcolor, hovercolor='0.975')
            button.label.set_fontsize(11)

            def apply_cutoffs(event):

                # get current cutoffs
                lowerCutoff, upperCutoff = update(val=None)

                # apply lower and upper cutoffs
                group_update = group[
                    (group[dna1] > lowerCutoff) &
                    (group[dna1] < upperCutoff)]

                # isolate x, y coordinates of selected centroids
                centroids = group_update[['Y_centroid', 'X_centroid']]

                # isolate cycle 1 DNA intensity values and assign
                # as quantitative point properties
                dna_intensity = group_update[dna1].values
                point_properties = {'dna_intensity': dna_intensity}

                # remove existing centroids and
                # plot new centroid selection in Napari window
                if not centroids.empty:
                    if len(viewer.layers) == 3:
                        viewer.layers.pop(2)
                    viewer.add_points(
                        centroids, name='Intensity',
                        properties=point_properties,
                        face_color='dna_intensity',
                        face_colormap='viridis',
                        edge_width=0.0, size=4.0)

            # add button functionality
            button.on_clicked(apply_cutoffs)

            # dock plot widget in Napari window
            viewer.window.add_dock_widget(
                hist_widget, name='DNA intensity histogram', area='right')

            hist_widget.setSizePolicy(
                QtWidgets.QSizePolicy.Fixed,
                QtWidgets.QSizePolicy.Preferred,
            )

            # show Napari window
            napari.run()

            # get current cutoffs
            lowerCutoff, upperCutoff = update(val=None)

            # plot DNA intensity histogram BEFORE filtering
            fig, ax = plt.subplots()
            plt.hist(
                group[dna1], bins=bins,
                density=False, color='b', ec='none',
                alpha=0.5, histtype=histtype,
                range=None, label='before')

            if lowerCutoff == upperCutoff:
                lowerCutoff = group[dna1].min()
                upperCutoff = group[dna1].max()
                print('No cutoffs applied, selecting all data points.')
            else:
                print(
                    f'Applied Lower Cutoff: {round(lowerCutoff, 2)}')
                print(
                    f'Applied Upper Cutoff: {round(upperCutoff, 2)}')

            # apply lower and upper cutoffs
            group_update = group[
                (group[dna1] > lowerCutoff) &
                (group[dna1] < upperCutoff)]

            # plot DNA intensity histogram AFTER filtering
            plt.hist(
                group_update[dna1], bins=bins, color='r', ec='none',
                alpha=0.5, histtype=histtype, range=None,
                label='after')
            plt.xlabel('mean DNA intensity')
            plt.ylabel('count')
            plt.title(f'Sample={name} mean DNA intensity)', size=10)

            legend_elements = []
            legend_elements.append(
                Line2D([0], [0], marker='o', color='none',
                       label='excluded data',
                       markerfacecolor='b', alpha=0.5,
                       markeredgecolor='none', lw=0.001,
                       markersize=8))
            legend_elements.append(
                Line2D([0], [0], marker='o', color='none',
                       label='included data',
                       markerfacecolor='r', alpha=0.5,
                       markeredgecolor='none', lw=0.001,
                       markersize=8))
            plt.legend(
                handles=legend_elements, prop={'size': 10},
                loc='best')

            plt.tight_layout()
            plt.savefig(
                os.path.join(intensity_dir, f'{name}.pdf'))
            plt.close('all')

            # isolate sample data to drop
            data_to_drop = group.copy()[
                (group[dna1] < lowerCutoff) |
                (group[dna1] > upperCutoff)]

            # create unique IDs for cells to drop in current sample
            data_to_drop['handle'] = (
                data_to_drop['CellID'].map(str) + '_' +
                data_to_drop['Sample'])

            # add sample indices to drop to idxs_to_drop dictionary
            idxs_to_drop[name] = [i for i in data_to_drop['handle']]

            # save updated idxs_to_drop dictionary as a csv
            dict_to_csv(
                dict=idxs_to_drop,
                path=os.path.join(intensity_dir, 'idxs_to_drop.csv'))

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

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def areaFilter(data, self, args):

        napari_warnings()

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # create area directory if it doesn't already exist
        area_dir = os.path.join(self.outDir, 'area')
        if not os.path.exists(area_dir):
            os.makedirs(area_dir)

        # set histogram bin size
        num_bins = 90

        # set histogram type
        histtype = 'stepfilled'

        # pick up where samples-loop left off
        if os.path.exists(
          os.path.join(area_dir, 'idxs_to_drop.csv')):

            # read stored dict
            idxs_to_drop = csv_to_dict(
                os.path.join(area_dir, 'idxs_to_drop.csv'))

            # dictionary values from strings to lists
            idxs_to_drop = {
                k: literal_eval(v) for k, v in idxs_to_drop.items()}

            # get names of samples already run
            previously_run_samples = idxs_to_drop.keys()

            # get names of samples remaining
            samples_remaining = (
                len(data['Sample'].unique())
                - len(previously_run_samples))

            print()
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

            print()
            print(f'Samples to threshold: {samples_remaining}')

            # initialize dictionary of sample indices to drop
            idxs_to_drop = {}

            # drop samples previously run (total samples in this case)
            df = data.copy()

            # natsort df by 'Sample' column
            df['Sample'] = pd.Categorical(
                df['Sample'], ordered=True,
                categories=natsorted(df['Sample'].unique()))
            df.sort_values('Sample', inplace=True)
            # convert 'Sample' column dtype back to string
            df['Sample'] = df['Sample'].astype(str)

        # loop over remaining samples
        for name, group in df.groupby('Sample', sort=False):

            print()
            print(f'Sample: {name}')

            # read segmentation outlines
            file_path = f'{self.inDir}/seg/{name}.*tif'

            # create image pyramid
            seg = single_channel_pyramid(
                glob.glob(file_path)[0], channel=0)

            # add DNA image to viewer
            viewer = napari.view_image(
                seg, rgb=False, visible=False, colormap='gray',
                opacity=0.5, name='segmentation'
                )

            # read DNA1 channel
            for file_path in glob.glob(
              f'{self.inDir}/tif/{name}.*tif'):
                # create image pyramid
                dna = single_channel_pyramid(file_path, channel=0)

            # add segmentation image to viewer
            viewer.add_image(
                dna, rgb=False, blending='additive',
                name=f'{dna1}: {name}'
                )

            # generate Qt widget to dock in Napari viewer
            hist_widget = QtWidgets.QWidget()

            # generate a blank figure canvas
            canvas = FigureCanvas(Figure())

            # construct vertical box layout object
            # to line up widgets vertically
            hist_layout = QtWidgets.QVBoxLayout(hist_widget)

            # add navigation tool bar and figure canvas to widget
            hist_layout.addWidget(NavigationToolbar(canvas, hist_widget))
            hist_layout.addWidget(canvas)

            # set plot style
            sns.set_style('whitegrid')

            # get figure object from canvas
            fig = canvas.figure

            # adjust plot on canvas to taste
            fig.subplots_adjust(left=0.25, bottom=0.25)

            # set plot title
            fig.suptitle(
                f'Sample={name} segmentation area', size=10)

            # get axis object from canvas
            ax = canvas.figure.subplots()

            # plot histogram on canvas axis
            n, bins, patches = ax.hist(
                group['Area'], bins=num_bins,
                density=False, color='grey', ec='none',
                alpha=0.75, histtype=histtype,
                range=None, label='before')

            # set y-axis label
            ax.set_ylabel('count')

            # add sliders to plot
            axcolor = 'lightgoldenrodyellow'
            axLowerCutoff = fig.add_axes(
                [0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
            axUpperCutoff = fig.add_axes(
                [0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

            # specify data range
            rnge = [bins.min(), bins.max()]

            # add slider functionality
            sLower = Slider(
                axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
                valinit=0.00, valstep=(rnge[1]/100000))
            sLower.label.set_fontsize(11)
            sLower.label.set_color('b')
            sUpper = Slider(
                axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
                valinit=0.00, valstep=(rnge[1]/100000))
            sUpper.label.set_fontsize(11)
            sUpper.label.set_color('r')

            # specify function for updating sliders
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

            # update sliders when moved
            sLower.on_changed(update)
            sUpper.on_changed(update)

            # add button to show selected centroids in Napari viewer
            button_ax = fig.add_axes([0.65, 0.025, 0.25, 0.06])
            button = Button(
                button_ax, 'Plot Points', color=axcolor, hovercolor='0.975')
            button.label.set_fontsize(11)

            def apply_cutoffs(event):

                # get current cutoffs
                lowerCutoff, upperCutoff = update(val=None)

                # apply lower and upper cutoffs
                group_update = group[
                    (group['Area'] > lowerCutoff) &
                    (group['Area'] < upperCutoff)]

                # isolate x, y coordinates of selected centroids
                centroids = group_update[['Y_centroid', 'X_centroid']]

                # isolate Area values and assign
                # as quantitative point properties
                cell_area = group_update['Area'].values
                point_properties = {'cell_area': cell_area}

                # remove existing centroids and
                # plot new centroid selection in Napari window
                if not centroids.empty:
                    if len(viewer.layers) == 3:
                        viewer.layers.pop(2)
                    viewer.add_points(
                        centroids, name='Area',
                        properties=point_properties,
                        face_color='cell_area',
                        face_colormap='viridis',
                        edge_width=0.0, size=4.0)

            # add button functionality
            button.on_clicked(apply_cutoffs)

            # dock plot widget in Napari window
            viewer.window.add_dock_widget(
                hist_widget, name='segmentation area histogram', area='right')

            hist_widget.setSizePolicy(
                QtWidgets.QSizePolicy.Fixed,
                QtWidgets.QSizePolicy.Preferred,
            )

            # show Napari window
            napari.run()

            # get current cutoffs
            lowerCutoff, upperCutoff = update(val=None)

            # plot DNA area histogram BEFORE filtering
            fig, ax = plt.subplots()
            plt.hist(
                group['Area'], bins=bins,
                density=False, color='b', ec='none',
                alpha=0.5, histtype=histtype,
                range=None, label='before')

            if lowerCutoff == upperCutoff:
                lowerCutoff = group['Area'].min()
                upperCutoff = group['Area'].max()
                print('No cutoffs applied, selecting all data points.')
            else:
                print(
                    f'Applied Lower Cutoff: {round(lowerCutoff, 2)}')
                print(
                    f'Applied Upper Cutoff: {round(upperCutoff, 2)}')

            # apply lower and upper cutoffs
            group_update = group[
                (group['Area'] > lowerCutoff) &
                (group['Area'] < upperCutoff)]

            # plot DNA area histogram AFTER filtering
            plt.hist(
                group_update['Area'], bins=bins, color='r', ec='none',
                alpha=0.5, histtype=histtype, range=None, label='after')
            plt.xlabel('mean DNA area')
            plt.ylabel('count')
            plt.title(f'Sample={name} cell area (pixel count))', size=10)

            legend_elements = []
            legend_elements.append(
                Line2D([0], [0], marker='o', color='none',
                       label='excluded data',
                       markerfacecolor='b', alpha=0.5,
                       markeredgecolor='none', lw=0.001,
                       markersize=8))
            legend_elements.append(
                Line2D([0], [0], marker='o', color='none',
                       label='included data',
                       markerfacecolor='r', alpha=0.5,
                       markeredgecolor='none', lw=0.001,
                       markersize=8))
            plt.legend(
                handles=legend_elements, prop={'size': 10},
                loc='best')

            plt.tight_layout()
            plt.savefig(os.path.join(area_dir, f'{name}.pdf'))
            plt.close('all')

            # isolate sample data to drop
            data_to_drop = group.copy()[
                (group['Area'] < lowerCutoff) |
                (group['Area'] > upperCutoff)]

            # create unique IDs for cells to drop in current sample
            data_to_drop['handle'] = (
                data_to_drop['CellID'].map(str) + '_' +
                data_to_drop['Sample'])

            # add sample indices to drop to idxs_to_drop dictionary
            idxs_to_drop[name] = [i for i in data_to_drop['handle']]

            # save updated idxs_to_drop dictionary as a csv
            dict_to_csv(
                dict=idxs_to_drop,
                path=os.path.join(area_dir, 'idxs_to_drop.csv'))

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

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def cycleCorrelation(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # create cycles directory if it doesn't already exist
        cycles_dir = os.path.join(self.outDir, 'cycles')
        if not os.path.exists(cycles_dir):
            os.makedirs(cycles_dir)

        # get ordered list of DNA cycles
        dna_cycles = natsorted(
            data.columns[data.columns.str.contains(dna_moniker)])

        # compute log(cycle 1/n) ratios
        ratios = pd.DataFrame(
            [np.log10((data[dna1] + 0.00001) /
             (data[i] + 0.00001)) for i in dna_cycles]).T

        # computing ratios changes columns headers,
        # create new ratio column headers
        unnamed_headers = [
            i for i in ratios.columns if i.startswith('Unnamed')]
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
                  value_name='log10(ratio)'))

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
                col='cycle', sharey=False)

            g = g.map(
                plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
                ec='none', range=None, bins=200, density=True)

            plt.savefig(
                os.path.join(cycles_dir, 'cycle_correlation(logRatio).pdf'))
            plt.close('all')

        if self.yAxisGating is True:

            # grab selected y-axis count cutoff if one was entered
            if os.path.exists(os.path.join(cycles_dir, 'count_cutoff.txt')):

                with open(
                  os.path.join(cycles_dir, 'count_cutoff.txt'), 'r') as f:
                    count_cutoff = f.readlines()
                    count_cutoff = float(count_cutoff[0])

            else:
                filename = os.path.join(
                    cycles_dir, 'cycle_correlation(logRatio).pdf')
                open_file(filename)

                def submit(text):

                    # save selected y-axis count cutoff as pickle
                    with open(
                      os.path.join(cycles_dir, 'count_cutoff.txt'), 'w') as f:
                        f.write(text)

                    count_cutoff = float(text)

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
                        col='cycle', sharey=False)
                    g = g.map(
                        plt.hist, 'log10(ratio)', color='r',
                        histtype='stepfilled', ec='none',
                        range=None,
                        bins=200, density=True)

                    for ax in g.axes.ravel():
                        ax.axhline(y=count_cutoff, c='k', linewidth=0.5)

                    plt.savefig(
                        os.path.join(
                            cycles_dir, 'cycle_correlation(logRatio).pdf')
                            )

                    filename = os.path.join(
                        cycles_dir, 'cycle_correlation(logRatio).pdf')
                    open_file(filename)

                    plt.show(block=False)
                    plt.close('all')

                plt.rcParams['figure.figsize'] = (6, 2)
                axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
                text_box = TextBox(
                    axbox, 'countCutoff', initial='',
                    color='0.95',
                    hovercolor='1.0',
                    label_pad=0.05)
                text_box.label.set_size(12)

                text_box.on_submit(submit)

                plt.show(block=True)

                # grab selected y-axis count cutoff if one was entered
                if os.path.exists(
                  os.path.join(cycles_dir, 'count_cutoff.txt')):
                    with open(
                      os.path.join(cycles_dir, 'count_cutoff.txt'), 'r') as f:
                        count_cutoff = f.readlines()
                        count_cutoff = float(count_cutoff[0])
                else:
                    # save count cutoff as pickle as 0.0
                    count_cutoff = '0.0'
                    with open(
                      os.path.join(cycles_dir, 'count_cutoff.txt'), 'w') as f:
                        f.write(count_cutoff)

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
                            bins=200, density=True)
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
                                count_indices[0], count_indices[0].max() + 1)]

                        if len(bin_values) > 1:
                            min_bin_val = min(bin_values)
                            max_bin_val = max(bin_values)

                            # get indices in log(ratio) series outside
                            # min_bin_val and max_bin_val
                            idxs = list(
                                cycle_data['index'][
                                    (cycle_data['log10(ratio)'] <
                                     min_bin_val)
                                    |
                                    (cycle_data['log10(ratio)'] >
                                     max_bin_val)])

                            # append indices of uncorrelated
                            # log(ratios) to idx_list
                            indices_to_drop.update(set(idxs))
            print()

            # filter dataframe by selecting indices NOT in the
            # indices_to_drop list
            print('Y-axis gating: Dropping unstable cells from samples...')
            df = data.loc[~data.index.isin(indices_to_drop)]
            plt.close('all')

        elif self.yAxisGating is False:

            napari_warnings()

            # pick up where samples loop left off
            if os.path.exists(os.path.join(cycles_dir, 'idxs_to_drop.csv')):

                # csv file might contain very huge fields, therefore increase
                # the field_size_limit, then read stored dict
                csv_module.field_size_limit(sys.maxsize)

                idxs_to_drop = csv_to_dict(
                    os.path.join(cycles_dir, 'idxs_to_drop.csv'))

                # dictionary values from strings to lists
                idxs_to_drop = {
                    k: literal_eval(v) for k, v in idxs_to_drop.items()}

                samples_to_threshold = (
                    len(data['Sample'].unique())
                    - len(idxs_to_drop.keys()))

                print()
                print(f'Samples to threshold: {samples_to_threshold}')

            else:
                # initialize a dictionary to append indices to drop
                samples_to_threshold = len(data['Sample'].unique())

                print()
                print(f'Samples to threshold: {samples_to_threshold }')

                idxs_to_drop = {}

            # drop samples previously run
            ratios_melt = ratios_melt[~ratios_melt['sample'].isin(
                    [i for i in idxs_to_drop.keys()])]

            if samples_to_threshold > 0:

                # show plot of all histograms
                # filename = os.path.join(
                #     cycles_dir, 'cycle_correlation(logRatio).pdf')
                # open_file(filename)

                for name, group in natsorted(ratios_melt.groupby(['sample'])):

                    print()
                    print(f'Sample: {name}')

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

                    # save dataframe of current cycle
                    cycle_data.to_parquet(
                        os.path.join(cycles_dir, 'cycle_data.parquet'))

                    # get channel number from markers.csv
                    channel_number = marker_channel_number(
                        markers, cycle_num)

                    # read segmentation outlines
                    file_path = f'{self.inDir}/seg/{name}.*tif'

                    # create image pyramid
                    seg = single_channel_pyramid(
                        glob.glob(file_path)[0], channel=0)

                    # add segmentation outlines image to viewer
                    viewer = napari.view_image(
                        seg, rgb=False, blending='additive',
                        opacity=0.5, colormap='gray', visible=False,
                        name='segmentation')

                    # read last DNA channel
                    for file_path in glob.glob(
                      f'{self.inDir}/tif/{name}.*tif'):
                        dna_last = single_channel_pyramid(
                            file_path, channel=channel_number.item() - 1)

                    # add last DNA image to viewer
                    viewer.add_image(
                        dna_last, rgb=False,
                        blending='additive',
                        colormap='magenta',
                        name=f'{cycle_num}')

                    # read first DNA channel
                    for file_path in glob.glob(
                      f'{self.inDir}/tif/{name}.*tif'):

                        # create image pyramid
                        dna_first = single_channel_pyramid(
                            file_path, channel=0)

                    # add first DNA image to viewer
                    viewer.add_image(
                        dna_first, rgb=False, blending='additive',
                        colormap='green',
                        name=f'{dna1}')

                    # generate Qt widget to dock in Napari viewer
                    widget = QtWidgets.QWidget()

                    # generate a blank figure canvas
                    canvas = FigureCanvas(Figure())

                    # construct vertical box layout object
                    # to line up widgets vertically
                    layout = QtWidgets.QVBoxLayout(widget)

                    # add navigation tool bar and figure canvas to widget
                    layout.addWidget(NavigationToolbar(canvas, widget))
                    layout.addWidget(canvas)

                    num_bins = 125
                    histtype = 'stepfilled'

                    # set plot style
                    sns.set_style('whitegrid')

                    # get figure object from canvas
                    fig = canvas.figure

                    # adjust plot on canvas to taste
                    fig.subplots_adjust(left=0.25, bottom=0.25)

                    # set plot title
                    log10 = '$log_{10}$'
                    fig.suptitle(
                        f'Sample={name}  {log10}({dna1}/{cycle_num})', size=10)

                    # get axis object from canvas
                    ax = canvas.figure.subplots()

                    # plot log(cycle 1/n) histogram for current sample
                    counts, bins, patches = ax.hist(
                        cycle_data['log10(ratio)'], bins=num_bins,
                        density=False, color='grey', ec='none',
                        alpha=0.75, histtype=histtype,
                        range=None, label='before')

                    ax.set_ylabel('count', size=13)
                    ax.tick_params(axis='x', labelsize=10)
                    ax.tick_params(axis='y', labelsize=10)

                    # add sliders to plot
                    axcolor = 'lightgoldenrodyellow'
                    axLowerCutoff = fig.add_axes(
                        [0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
                    axUpperCutoff = fig.add_axes(
                        [0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

                    # specify data range
                    rnge = [bins.min(), bins.max()]

                    # add slider functionality
                    sLower = Slider(
                        axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
                        valinit=0.00, valstep=(rnge[1]/100000))
                    sLower.label.set_fontsize(11)
                    sLower.label.set_color('b')
                    sUpper = Slider(
                        axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
                        valinit=0.00, valstep=(rnge[1]/100000))
                    sUpper.label.set_fontsize(11)
                    sUpper.label.set_color('r')

                    # specify function for updating sliders
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

                    # update sliders when moved
                    sLower.on_changed(update)
                    sUpper.on_changed(update)

                    # add button to show selected centroids in Napari viewer
                    button_ax = fig.add_axes([0.65, 0.025, 0.25, 0.06])
                    button = Button(
                        button_ax, 'Plot Points',
                        color=axcolor, hovercolor='0.975')
                    button.label.set_fontsize(11)

                    def apply_cutoffs(event):

                        # get current cutoffs
                        lowerCutoff, upperCutoff = update(val=None)

                        # read group
                        cycle_data = pd.read_parquet(
                            os.path.join(cycles_dir, 'cycle_data.parquet'))

                        # get indices in log(ratio) series outside
                        # lower and upper cutoffs
                        idxs = list(
                            cycle_data['index'][
                                (cycle_data['log10(ratio)']
                                 < lowerCutoff) |
                                (cycle_data['log10(ratio)']
                                 > upperCutoff)])

                        # filter group data by selecting
                        # indices NOT in idxs
                        sample_data = data[data['Sample'] == name]
                        drop_df = sample_data.index.isin(idxs)
                        centroids = sample_data[
                            ['Y_centroid', 'X_centroid']][~drop_df]

                        # remove existing centroids and
                        # plot new centroid selection in Napari window
                        if not centroids.empty:
                            if len(viewer.layers) == 4:
                                viewer.layers.pop(3)
                            viewer.add_points(
                                centroids,
                                name='Selected Cells',
                                properties=None,
                                face_color='yellow',
                                edge_color='k',
                                edge_width=0.0, size=7.0)

                    # add button functionality
                    button.on_clicked(apply_cutoffs)

                    # dock plot widget in Napari window
                    viewer.window.add_dock_widget(
                        widget, name=f'first/last DNA histogram',
                        area='right')

                    # show Napari window
                    napari.run()

                    # get final lower and upper cutoffs
                    lowerCutoff, upperCutoff = update(val=None)

                    # read dataframe of current cycle and apply final cutoffs
                    cycle_data = pd.read_parquet(
                        os.path.join(cycles_dir, 'cycle_data.parquet'))

                    # take all data if sliders not moved
                    if lowerCutoff == upperCutoff:
                        lowerCutoff = cycle_data['log10(ratio)'].min()
                        upperCutoff = cycle_data['log10(ratio)'].max()
                        print(
                            'No cutoffs applied, selecting all data points.')
                    else:
                        print(
                            f'Applied Lower Cutoff: {round(lowerCutoff, 2)}')
                        print(
                            f'Applied Upper Cutoff: {round(upperCutoff, 2)}')

                    # otherwise, get indices outside lower and upper cutoffs
                    # and append to idxs_to_drop dictionary
                    idxs = list(
                        cycle_data['index'][
                            (cycle_data['log10(ratio)'] < lowerCutoff)
                            |
                            (cycle_data['log10(ratio)'] > upperCutoff)])
                    idxs_to_drop[name] = idxs

                    # save updated idxs_to_drop dictionary as a csv
                    dict_to_csv(
                        dict=idxs_to_drop,
                        path=os.path.join(cycles_dir, 'idxs_to_drop.csv'))
            print()

            # filter data with indices to frop from all samples
            print('X-axis gating: Dropping unstable cells from samples.')
            indices_to_drop = set()
            for k, v in idxs_to_drop.items():
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
            .melt(id_vars=['Sample', 'index'], var_name='cycle'))

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
            ['Sample', 'cycle', 'index'])

        # plot dna intensity correlation per cycle
        fig, ax = plt.subplots(figsize=(5, 5))
        g = sns.FacetGrid(
            facet_per_cycle_melt, col='cycle', col_wrap=5,
            sharex=True, sharey=False)

        g.map(
            lambda y, color: plt.scatter(
                facet_per_cycle_melt['value'].loc[
                    facet_per_cycle_melt['cycle']
                    == dna1], y, s=0.05, alpha=0.1, linewidth=None,
                marker='o', c='r'), 'value')

        plt.savefig(
            os.path.join(
                cycles_dir, 'cycle_correlation(perCycle).png'), dpi=600)
        plt.close('all')

        # plot dna intensity correlation per cycle (color by sample)
        fig, ax = plt.subplots(figsize=(5, 5))

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(
                facet_per_cycle_melt['Sample'].unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False)

        sample_color_dict = dict(zip(
            natsorted(facet_per_cycle_melt['Sample'].unique()),
            cmap.colors))

        g = sns.FacetGrid(
            facet_per_cycle_melt, col='cycle', hue='Sample',
            col_wrap=5, sharex=True, sharey=True)

        g.map(
            lambda sam, y, color, **kwargs: plt.scatter(
                facet_per_cycle_melt.loc[
                    (facet_per_cycle_melt['Sample'] ==
                     sam.unique()[0])
                    & (facet_per_cycle_melt['cycle'] == dna1),
                    'value'], y,
                c=np.reshape(sample_color_dict[sam.unique()[0]], (-1, 3)),
                s=0.05, linewidth=None, marker='o', **kwargs),
            'Sample', 'value')

        plt.legend(markerscale=10, bbox_to_anchor=(1.1, 1.05))

        plt.savefig(
            os.path.join(
                cycles_dir, 'cycle_correlation(perSample).png'), dpi=600,
            bbox_inches='tight')
        plt.close('all')
        print()

        # remove last sample groupby dataframe
        if os.path.exists(os.path.join(cycles_dir, 'cycle_data.parquet')):
            os.remove(os.path.join(cycles_dir, 'cycle_data.parquet'))

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        return df

    @module
    def logTransform(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        abx_channels_mod = data[abx_channels].copy()
        abx_channels_mod += 0.00000000001
        abx_channels_mod = np.log10(abx_channels_mod)
        data.loc[:, abx_channels] = abx_channels_mod

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def pruneOutliers(data, self, args):

        napari_warnings()

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # create pruning directory if it doesn't already exist
        pruning_dir = os.path.join(self.outDir, 'pruning')
        if not os.path.exists(pruning_dir):
            os.mkdir(pruning_dir)

        # if percentile cutoffs have already been assigned for some antibodies
        # open pruning dict, get list of remaining markers to prune
        # sorted by markers.csv, and read partially-pruned data
        if os.path.exists(os.path.join(pruning_dir, 'pruning_dict.csv')):

            # read stored dict
            pruning_dict = csv_to_dict(
                os.path.join(pruning_dir, 'pruning_dict.csv'))

            # dictionary tuple values are converted to strings when saved,
            # convert back to tuples
            pruning_dict = {
                k: tuple(map(
                    float, v.lstrip('(').rstrip(')').split(', ')))
                for (k, v) in pruning_dict.items()}

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

            print()
            print(f'Immunomarker channels to prune: {len(markers_to_prune)}')

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

            print()
            print(f'Immunomarker channels to prune: {len(markers_to_prune)}')

        # view raw signal intensities for all samples
        for ab in markers_to_prune:

            print()
            print(f'Working on channel: {ab}')

            # initialize Napari window without an image
            # before first cutoffs are applied
            viewer = napari.Viewer()

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
                hist_facet['Condition'])

            # set plot style
            sns.set_style('white')

            # plot raw facets
            g_raw = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=4,
                height=1, aspect=1.0, sharex=True, sharey=False)

            # use hexbins for plotting
            if self.hexbins:
                g_raw.map(
                    plt.hexbin, 'signal', 'Area',
                    gridsize=self.hexbinGridSize,
                    linewidths=0.02, color='dimgrey')

            # use scatter plots for plotting
            else:
                g_raw.map(
                    plt.scatter, 'signal', 'Area', s=0.05,
                    linewidths=0.0, color='k')

                # suppresses matplotlib tight_layout warning
                g_raw.fig.set_figheight(2.5)
                g_raw.fig.set_figwidth(7)

            g_raw.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=5.0, pad=4.0)

            for ax in g_raw.axes.flatten():
                ax.tick_params(
                    axis='both', which='major',
                    labelsize=5.0, pad=-2)

                ax.xaxis.label.set_size(5.0)
                ax.yaxis.label.set_size(5.0)

                if self.hexbins:
                    ax.spines['left'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                else:
                    ax.spines['left'].set_linewidth(0.1)
                    ax.spines['bottom'].set_linewidth(0.1)

            # generate Qt widget
            raw_widget = QtWidgets.QWidget()

            # add FacetGrid object to a blank canvas
            raw_canvas = FigureCanvas(g_raw.fig)

            # construct vertical box layout object
            # to line up widgets vertically
            raw_layout = QtWidgets.QVBoxLayout(raw_widget)

            # add navigation tool bar and figure canvas to widget
            raw_layout.addWidget(NavigationToolbar(raw_canvas, raw_widget))
            raw_layout.addWidget(raw_canvas)

            plt.subplots_adjust(
                left=0.03, bottom=0.13, right=0.99,
                top=0.92, hspace=0.8, wspace=0.3)

            plt.savefig(
                os.path.join(pruning_dir, f'{ab}_raw.png'), dpi=300,
                bbox_inches='tight')

            # generate Qt widget
            pruned_widget = QtWidgets.QWidget()

            # construct vertical box layout object
            pruned_layout = QtWidgets.QVBoxLayout(pruned_widget)

            @magicgui(
                layout='vertical',
                call_button='Apply Cutoffs',
                lower_cutoff={'widget_type': 'FloatSlider', 'max': 100.0,
                              'label': 'Lower Percentile'},
                upper_cutoff={'widget_type': 'FloatSlider', 'max': 100.0,
                              'label': 'Upper Percentile'},
                )
            def select_parameters(
              lower_cutoff: float = 0.0,
              upper_cutoff: float = 100.0,
            ):
                # close exisiting plots
                plt.close('all')

                # add entered cutoffs to pruning_dict
                pruning_dict[ab] = (
                    lower_cutoff,
                    upper_cutoff)

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
                            sample_channel_data, lower_cutoff)]
                    indices_to_drop.extend(low_drop_idxs)

                    high_drop_idxs = sample_channel_data.index[
                        sample_channel_data > np.percentile(
                            sample_channel_data, upper_cutoff)]
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
                g_pruned = sns.FacetGrid(
                    hist_facet, col='for_plot', col_wrap=4,
                    height=1.0, aspect=1.0, sharex=True, sharey=False)

                # use hexbins for plotting
                if self.hexbins:
                    g_pruned.map(
                        plt.hexbin, 'signal', 'Area',
                        gridsize=self.hexbinGridSize,
                        linewidths=0.02, color='dimgrey')

                # use scatter plots for plotting
                else:
                    g_pruned.map(
                        plt.scatter, 'signal', 'Area', s=0.05,
                        linewidths=0.0, color='k')

                    # suppresses matplotlib tight_layout warning
                    g_pruned.fig.set_figheight(2.5)
                    g_pruned.fig.set_figwidth(7)

                g_pruned.set_titles(
                    col_template="{col_name}", fontweight='bold',
                    size=5.0, pad=4.0)

                for ax in g_pruned.axes.flatten():
                    ax.tick_params(
                        axis='both', which='major',
                        labelsize=5.0, pad=-2)

                    ax.xaxis.label.set_size(5.0)
                    ax.yaxis.label.set_size(5.0)

                    if self.hexbins:
                        ax.spines['left'].set_visible(False)
                        ax.spines['bottom'].set_visible(False)
                    else:
                        ax.spines['left'].set_linewidth(0.1)
                        ax.spines['bottom'].set_linewidth(0.1)

                plt.subplots_adjust(
                    left=0.03, bottom=0.13, right=0.99,
                    top=0.92, hspace=0.8, wspace=0.3)

                count = pruned_layout.count()
                # print('layout count:', count)
                # print('widget children:', pruned_widget.children())

                # remove old widgets from widget layout
                for i in range(count - 1, -1, -1):
                    item = pruned_layout.itemAt(i)
                    widget = item.widget()
                    # print('    item:', item)
                    # print('        widget:', widget)
                    if widget:
                        widget.setParent(None)

                # add updated widgets to widget layout
                pruned_canvas = FigureCanvas(g_pruned.fig)
                pruned_layout.addWidget(
                    NavigationToolbar(pruned_canvas, pruned_widget))
                pruned_layout.addWidget(pruned_canvas)

                ###############################################################
                @magicgui(
                    layout='horizontal',
                    call_button='View Outliers',
                    sample_name={'label': 'Sample Name'},
                    )
                def sample_selector(
                  sample_name: str,
                  ):

                    return sample_name

                ###############################################################

                ###############################################################
                @sample_selector.called.connect
                def my_callback(value: str):

                    channel_number = marker_channel_number(markers, ab)

                    if value in data_copy1['Sample'].unique():

                        print()
                        print(f'Sample selection: {value}')

                        # read DNA1 channel
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{value}.*tif'):
                            dna = single_channel_pyramid(file_path, channel=0)

                        # read antibody channel
                        for file_path in glob.glob(
                          f'{self.inDir}/tif/{value}.*tif'):
                            channel = single_channel_pyramid(
                                file_path, channel=channel_number.item() - 1)

                        # read segmentation outlines
                        file_path = f'{self.inDir}/seg/{value}.*tif'
                        seg = single_channel_pyramid(
                            glob.glob(file_path)[0], channel=0)

                        # remove existing layers
                        for i in reversed(range(len(viewer.layers))):
                            viewer.layers.pop(i)

                        # decorate Napari viewer
                        viewer.add_image(
                            dna, rgb=False, opacity=0.5, name=dna1)

                        viewer.add_image(
                            channel, rgb=False, blending='additive',
                            colormap='green', visible=False, name=ab)

                        viewer.add_image(
                            seg, rgb=False, blending='additive', opacity=0.5,
                            colormap='gray', visible=False,
                            name='segmentation'
                            )

                        # grab centroids of low signal intensity outliers
                        low_centroids = data_copy1[
                            ['Y_centroid', 'X_centroid']][
                            (data_copy1.index.isin(total_low_idxs)) &
                            (data_copy1['Sample'] == value)]

                        # grab centroids of high signal intensity outliers
                        high_centroids = data_copy1[
                            ['Y_centroid', 'X_centroid']][
                            (data_copy1.index.isin(total_high_idxs)) &
                            (data_copy1['Sample'] == value)]

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
                    else:
                        print()
                        print('Invalid sample name entered.')
                        pass

                sample_selector.native.setSizePolicy(
                    QtWidgets.QSizePolicy.Maximum,
                    QtWidgets.QSizePolicy.Maximum,
                )

                pruned_layout.addWidget(sample_selector.native)

            viewer.window.add_dock_widget(
                select_parameters, name='select percentile cutoffs',
                area='right')

            select_parameters.native.setSizePolicy(
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Maximum,
            )

            raw_widget.setSizePolicy(
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Fixed,
            )

            pruned_widget.setSizePolicy(
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Fixed,
            )

            raw_dock = viewer.window.add_dock_widget(
                raw_widget, name=f'{ab} raw', area='right')

            pruned_dock = viewer.window.add_dock_widget(
                pruned_widget, name=f'{ab} pruned/rescaled', area='right')

            napari.run()

            # take all data if cutoff window is closed without entering values
            if ab not in pruning_dict.keys():
                print()
                print(f'{ab}:')
                print('No cutoffs applied, selecting all data points.')
                pruning_dict[ab] = (0.0, 100.0)
            else:
                print()
                print(f'{ab}:')
                print(
                    'Applied lower percentile cutoff:', pruning_dict[ab][0])
                print('Applied upper percentile cutoff:', pruning_dict[ab][1])

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

            # save updated pruning_dict as a csv
            dict_to_csv(
                dict=pruning_dict,
                path=os.path.join(pruning_dir, 'pruning_dict.csv'))

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

            for sample in natsorted(data['Sample'].unique()):

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

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def metaQC(data, self, args):

        print()
        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # drop antibody channel exclusions for metaQC clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClusteringQC]

        # create metaQC directory if metaQC is to be performed and
        # the directory doesn't already exist
        reclass_dir = os.path.join(
            self.outDir, 'metaQC')
        if not os.path.exists(reclass_dir):
            os.makedirs(reclass_dir)

        # specify the names of modules in the pipeline that perform
        # data redaction prior to the metQC module
        modules = ['aggregateData', 'selectROIs', 'intensityFilter',
                   'areaFilter', 'cycleCorrelation', 'pruneOutliers']

        # build a dictionary of data returned by each module (clean)
        module_dict = {}
        for module_idx, module in enumerate(modules):
            data = pd.read_parquet(
                os.path.join(
                    self.outDir, f'checkpoints/{module}.parquet'))
            module_dict[module_idx] = [module, data]

        #######################################################################
        # build QCData: QCData is a combination of clean and noisy data

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
                        noisyData = pd.concat(
                            [noisyData, module_dict[module_idx][2]], axis=0
                            )
            else:
                # if postive ROI selection, exlcude cells not gated on
                noisyData = pd.DataFrame()
                for module_idx, module in enumerate(modules):
                    if len(module_dict[module_idx]) == 3:
                        if not module == 'selectROIs':
                            noisyData = pd.concat(
                                [noisyData, module_dict[module_idx][2]],
                                axis=0
                                )

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
                with open(
                  os.path.join(reclass_dir, 'chunk_index.txt'), 'r') as f:
                    chunk_index = f.readlines()
                    chunk_index = int(chunk_index[0])

                # read current chunk data if it exists
                if os.path.exists(os.path.join(reclass_dir, 'chunk.pkl')):
                    f = open(os.path.join(
                        reclass_dir, 'chunk.pkl'), 'rb')
                    chunk = pickle.load(f)
            else:
                # if QCData.pkl doesn't exist, append noisyData
                # to cleanDataRaw, row-wise
                QCData = pd.concat([cleanDataRaw, noisyData], axis=0)

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
                chunk_index = '0'
                with open(
                  os.path.join(reclass_dir, 'chunk_index.txt'), 'w') as f:
                    f.write(chunk_index)
                chunk_index = int(chunk_index)

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
            # the chunk.pkl which has clustering labels if it exits.
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
            ###################################################################

            # loop over QCData chunks
            for chunk in chunks[chunk_index:]:

                print(
                    f'Clustering: {chunk_index + 1} of ' +
                    f'{len(chunks)} data chunks')

                # make directory for current chunk if it hasn't already
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
                        str(datetime.now() - startTime))

                    np.save(
                        os.path.join(chunk_dir, 'embedding'),
                        embedding)

                    chunk['emb1'] = embedding[:, 0]
                    chunk['emb2'] = embedding[:, 1]

                # define the point size for cells in the embedding
                point_size = 50000/len(chunk)

                def reclassify_chunk(chunk, clean_cutoff, noisy_cutoff):

                    clean = pd.DataFrame()
                    noisy = pd.DataFrame()
                    for name, cluster in chunk.groupby(
                      f'cluster_{self.dimensionEmbeddingQC}d'):
                        if name != -1:
                            # if a cluster contains
                            # >= n% clean data,
                            # reclassify all clustering
                            # cells as noisy
                            if (
                              (len(cluster[cluster[
                                'QC_status'] == 'clean'])
                               / len(cluster)) >=
                              clean_cutoff):
                                clean = pd.concat([clean, cluster], axis=0)
                            # elif a cluster contains
                            # >= n% noisy data,
                            # reclassify all clustering
                            # cells as clean
                            elif (
                              (len(cluster[cluster[
                                'QC_status'] == 'noisy'])
                               / len(cluster)) >=
                              noisy_cutoff):
                                noisy = pd.concat([noisy, cluster], axis=0)
                            else:
                                noisy = pd.concat(
                                    [noisy,
                                     cluster[cluster['QC_status'] == 'noisy']],
                                    axis=0
                                    )
                                clean = pd.concat(
                                    [clean,
                                     cluster[cluster['QC_status'] == 'clean']],
                                    axis=0
                                    )

                    # consider -1 cells from clean data
                    # to be "noisy"
                    clean_outliers = chunk[
                        (chunk[f'cluster_{self.dimensionEmbeddingQC}d'] == -1)
                        &
                        (chunk['QC_status'] == 'clean')
                        ].copy()
                    noisy = pd.concat([noisy, clean_outliers], axis=0)

                    # consider -1 cells from noisy data
                    # to be "noisy"
                    noisy_outliers = chunk[
                        (chunk[f'cluster_{self.dimensionEmbeddingQC}d'] == -1)
                        &
                        (chunk['QC_status'] == 'noisy')
                        ].copy()
                    noisy = pd.concat([noisy, noisy_outliers], axis=0)

                    chunk['Reclass'] = 'init'
                    chunk.loc[
                        chunk.index.isin(clean.index),
                        'Reclass'] = 'clean'
                    chunk.loc[
                        chunk.index.isin(noisy.index),
                        'Reclass'] = 'noisy'

                    return chunk, clean, noisy

                # interact with plots to identify optimal min cluster size
                while not os.path.isfile(os.path.join(reclass_dir, 'MCS.txt')):

                    # initial Napari viewer without images
                    viewer = napari.Viewer()

                    # generate Qt widget
                    cluster_widget = QtWidgets.QWidget()

                    # generate vertical widget layout
                    cluster_layout = QtWidgets.QVBoxLayout(cluster_widget)

                    cluster_widget.setSizePolicy(
                        QtWidgets.QSizePolicy.Minimum,
                        QtWidgets.QSizePolicy.Maximum,
                    )

                    ###########################################################
                    @magicgui(
                        layout='horizontal',
                        call_button='Cluster and Plot',
                        MCS={'label': 'Min Cluster Size (MCS)'},
                        )
                    def cluster_and_plot(MCS: int = 200.0):

                        # placeholder for lasso selection
                        selector = None

                        sns.set_style('whitegrid')

                        fig = plt.figure(figsize=(8, 7))
                        matplotlib_warnings(fig)

                        gs = plt.GridSpec(2, 4, figure=fig)

                        # define axes
                        ax_cluster = fig.add_subplot(gs[0, 0])
                        ax_status = fig.add_subplot(gs[0, 1])
                        ax_reclass = fig.add_subplot(gs[0, 2])
                        ax_sample = fig.add_subplot(gs[0, 3])

                        ax_cluster_lbs = fig.add_subplot(gs[1, 0])
                        ax_status_lbs = fig.add_subplot(gs[1, 1])
                        ax_reclass_lbs = fig.add_subplot(gs[1, 2])
                        ax_sample_lbs = fig.add_subplot(gs[1, 3])

                        plt.subplots_adjust(
                            left=0.01, right=0.99, bottom=0.0,
                            top=0.9, wspace=0.0, hspace=0.0)

                        clustering = hdbscan.HDBSCAN(
                            min_cluster_size=MCS,
                            min_samples=None,
                            metric='euclidean', alpha=1.0, p=None,
                            algorithm='best', leaf_size=40,
                            memory=Memory(location=None),
                            approx_min_span_tree=True,
                            gen_min_span_tree=False,
                            core_dist_n_jobs=-1,
                            cluster_selection_method='eom',
                            allow_single_cluster=False,
                            prediction_data=False,
                            match_reference_implementation=False).fit(
                                chunk[['emb1', 'emb2']])

                        chunk[f'cluster_{self.dimensionEmbeddingQC}d'] = (
                            clustering.labels_
                            )

                        # scatter point selection tool assumes a
                        # sorted index, but the index of QCdata is
                        # shuffled to acheive a mix of clean and
                        # noisy data per chunk
                        chunk.sort_index(inplace=True)

                        print()
                        print(
                            f'min_cluster_size = {MCS}',
                            np.unique(clustering.labels_))

                        # PLOT embedding
                        for color_by in [
                          f'cluster_{self.dimensionEmbeddingQC}d',
                          'QC_status', 'Reclass', 'Sample']:

                            highlight = 'none'

                            if color_by == f'cluster_{self.dimensionEmbeddingQC}d':

                                # build cmap
                                cmap = categorical_cmap(
                                    numUniqueSamples=len(
                                        chunk[color_by].unique()),
                                    numCatagories=10, cmap='tab10',
                                    continuous=False)

                                # make black the first color to specify
                                # unclustered cells (cluster -1)
                                if -1 in chunk[color_by].unique():

                                    cmap = ListedColormap(
                                        np.insert(
                                            arr=cmap.colors, obj=0,
                                            values=[0, 0, 0], axis=0))

                                    # trim cmap to # unique samples
                                    trim = (
                                        len(cmap.colors) - len(
                                            chunk[color_by].unique()))
                                    cmap = ListedColormap(
                                        cmap.colors[:-trim])

                                sample_dict = dict(
                                    zip(natsorted(
                                            chunk[color_by].unique()),
                                        list(range(len(chunk[color_by]
                                             .unique())))))

                                c = [sample_dict[i] for i
                                     in chunk[color_by]]

                                ax_cluster.cla()

                                cluster_paths = ax_cluster.scatter(
                                    chunk['emb1'], chunk['emb2'], c=c,
                                    alpha=1.0, s=point_size,
                                    cmap=cmap, ec='k', linewidth=0.0)

                                ax_cluster.set_title(
                                    'HDBSCAN (lasso)', fontsize=10)
                                ax_cluster.axis('equal')
                                ax_cluster.axes.xaxis.set_visible(False)
                                ax_cluster.axes.yaxis.set_visible(False)
                                ax_cluster.grid(False)

                                legend_elements = []
                                for e, i in enumerate(
                                  natsorted(chunk[color_by].unique())):

                                    norm_ax, hi_markers = cluster_expression(
                                        df=chunk, markers=abx_channels,
                                        cluster=i, num_proteins=3,
                                        clus_dim=self.dimensionEmbeddingQC,
                                        norm_ax=self.topMarkersQC
                                        )

                                    legend_elements.append(
                                        Line2D([0], [0], marker='o',
                                               color='none',
                                               label=(
                                                f'Cluster {i}: '
                                                f'{hi_markers} {norm_ax}'),
                                               markerfacecolor=(
                                                cmap.colors[e]),
                                               markeredgecolor='none',
                                               lw=0.001, markersize=2))

                                cluster_lgd = ax_cluster_lbs.legend(
                                    handles=legend_elements,
                                    prop={'size': 5}, loc='upper left',
                                    frameon=False)

                                ax_cluster_lbs.axis('off')

                            elif color_by == 'QC_status':

                                # build cmap
                                cmap = ListedColormap(
                                    np.array([[0.91, 0.29, 0.235],
                                              [0.18, 0.16, 0.15]]))

                                sample_dict = dict(
                                    zip(natsorted(
                                        chunk['QC_status'].unique()),
                                        list(range(len(chunk['QC_status']
                                             .unique())))))

                                c = [sample_dict[i] for
                                     i in chunk['QC_status']]

                                ax_status.cla()

                                ax_status.scatter(
                                    chunk['emb1'], chunk['emb2'], c=c,
                                    cmap=cmap, alpha=1.0, s=point_size,
                                    ec='k', linewidth=0.0)

                                ax_status.set_title(
                                    'QC Status', fontsize=10)
                                ax_status.axis('equal')
                                ax_status.axes.xaxis.set_visible(False)
                                ax_status.axes.yaxis.set_visible(False)
                                ax_status.grid(False)

                                legend_elements = []
                                for e, i in enumerate(
                                    natsorted(
                                        chunk['QC_status'].unique())):

                                    if i == highlight:
                                        markeredgecolor = 'k'
                                    else:
                                        markeredgecolor = 'none'

                                    sample_to_map = (
                                        chunk['Sample'][
                                            chunk['QC_status'] == i]
                                        .unique()[0])

                                    legend_elements.append(
                                        Line2D([0], [0], marker='o',
                                               color='none',
                                               label=i,
                                               markerfacecolor=(
                                                cmap.colors[e]),
                                               markeredgecolor=(
                                                markeredgecolor),
                                               lw=0.001,
                                               markersize=2))

                                qc_lgd = ax_status_lbs.legend(
                                    handles=legend_elements,
                                    prop={'size': 5}, loc='upper left',
                                    frameon=False)

                                ax_status_lbs.axis('off')

                            elif color_by == 'Reclass':

                                # build cmap
                                cmap = ListedColormap(
                                    np.array([[0.91, 0.29, 0.235],
                                              [0.18, 0.16, 0.15]]))

                                sample_dict = dict(
                                    zip(natsorted(
                                            chunk['QC_status'].unique()),
                                        list(range(len(chunk['QC_status']
                                             .unique())))))

                                c = [sample_dict[i] for
                                     i in chunk['QC_status']]

                                ax_reclass.cla()

                                ax_reclass.scatter(
                                    chunk['emb1'], chunk['emb2'], c=c,
                                    cmap=cmap, alpha=1.0, s=point_size,
                                    ec='k', linewidth=0.0)

                                ax_reclass.set_title(
                                    'Reclassification', fontsize=10)
                                ax_reclass.axis('equal')
                                ax_reclass.axes.xaxis.set_visible(False)
                                ax_reclass.axes.yaxis.set_visible(False)
                                ax_reclass.grid(False)

                                legend_elements = []
                                for e, i in enumerate(
                                  natsorted(chunk['QC_status'].unique())):

                                    if i == highlight:
                                        markeredgecolor = 'k'
                                    else:
                                        markeredgecolor = 'none'

                                    sample_to_map = (
                                        chunk['Sample'][
                                            chunk['QC_status'] == i]
                                        .unique()[0])

                                    legend_elements.append(
                                        Line2D([0], [0], marker='o',
                                               color='none',
                                               label=i,
                                               markerfacecolor=(
                                                cmap.colors[e]),
                                               markeredgecolor=(
                                                markeredgecolor),
                                               lw=0.001,
                                               markersize=2))

                                reclass_lgd = ax_reclass_lbs.legend(
                                    handles=legend_elements,
                                    prop={'size': 5}, loc='upper left',
                                    frameon=False)

                                ax_reclass_lbs.axis('off')

                            elif color_by == 'Sample':

                                # build cmap
                                cmap = categorical_cmap(
                                    numUniqueSamples=len(
                                        chunk['Sample'].unique()),
                                    numCatagories=10,
                                    cmap='tab10',
                                    continuous=False)

                                sample_dict = dict(
                                    zip(natsorted(
                                            chunk['Sample'].unique()),
                                        list(range(len(chunk['Sample']
                                             .unique())))))

                                c = [sample_dict[i] for
                                     i in chunk['Sample']]

                                ax_sample.cla()

                                ax_sample.scatter(
                                    chunk['emb1'], chunk['emb2'], c=c,
                                    cmap=cmap, alpha=1.0, s=point_size,
                                    ec='k', linewidth=0.0)

                                ax_sample.set_title('Sample', fontsize=10)
                                ax_sample.axis('equal')
                                ax_sample.axes.xaxis.set_visible(False)
                                ax_sample.axes.yaxis.set_visible(False)
                                ax_sample.grid(False)

                                legend_elements = []
                                for e, i in enumerate(
                                  natsorted(chunk['Sample'].unique())):

                                    if i == highlight:
                                        markeredgecolor = 'k'
                                    else:
                                        markeredgecolor = 'none'

                                    sample_to_map = (
                                        chunk['Sample'][
                                            chunk['Sample'] == i]
                                        .unique()[0])

                                    legend_elements.append(
                                        Line2D([0], [0], marker='o',
                                               color='none',
                                               label=i,
                                               markerfacecolor=(
                                                cmap.colors[e]),
                                               markeredgecolor=(
                                                markeredgecolor),
                                               lw=0.001,
                                               markersize=2))

                                sample_lgd = ax_sample_lbs.legend(
                                    handles=legend_elements,
                                    prop={'size': 5}, loc='upper left',
                                    frameon=False)

                                ax_sample_lbs.axis('off')

                        count = cluster_layout.count()
                        # print('layout count:', count)
                        # print('widget children:', pruned_widget.children())

                        # remove old widgets from widget layout
                        for i in range(count - 1, -1, -1):
                            item = cluster_layout.itemAt(i)
                            widget = item.widget()
                            # print('    item:', item)
                            # print('        widget:', widget)
                            if widget:
                                widget.setParent(None)

                        # add updated widgets to widget layout
                        cluster_canvas = FigureCanvas(fig)
                        cluster_layout.addWidget(
                            NavigationToolbar(cluster_canvas, cluster_widget))
                        cluster_layout.addWidget(cluster_canvas)

                        # must call draw() before creating selector,
                        # or alpha setting doesn't work.
                        fig.canvas.draw()

                        if selector:
                            selector.disconnect()
                        selector = SelectFromCollection(
                            ax_cluster, cluster_paths)

                        #######################################################
                        @magicgui(
                            layout='horizontal',
                            call_button='View Lassoed Points',
                            sample_name={'label': 'Sample Name'},
                            )
                        def sample_selector(
                          sample_name: str,
                          ):

                            return sample_name
                        #######################################################

                        #######################################################
                        @sample_selector.called.connect
                        def my_callback(value: str):

                            print()
                            print(f'Sample selection: {value}')

                            # if cells lassoed
                            if not (
                              selector.ind is not None and
                              len(selector.ind) >= 1
                              ):

                                print()
                                print(
                                    'Cells must be lassoed in HDBSCAN plot '
                                    'before sample inspection.')
                                pass

                            else:
                                # if valid sample name entered
                                if value in chunk['Sample'].unique():

                                    # show highest expression channels
                                    lasso = chunk.loc[
                                        chunk.index.isin(selector.ind)
                                        ].copy()

                                    # assign lassoed data a dummy
                                    # cluster variable and get highest
                                    # expressed markers across channels
                                    lasso[
                                        f'cluster_{self.dimensionEmbeddingQC}d'
                                        ] = 1000

                                    norm_ax, hi_markers = cluster_expression(
                                        df=lasso, markers=abx_channels,
                                        cluster=1000, num_proteins=3,
                                        clus_dim=self.dimensionEmbeddingQC,
                                        norm_ax='channels'
                                        )

                                    print()
                                    print(
                                        'Top three expressed '
                                        f'markers {hi_markers} channels'
                                        )

                                    # superimpose centroids of lassoed noisy
                                    # cells colored by stage removed over
                                    # channel images
                                    print()
                                    print(
                                        f'Opening sample {value} in Napari...')

                                    # remove previous samples' image layers
                                    for i in reversed(
                                      range(0, len(viewer.layers))):
                                        viewer.layers.pop(i)

                                    napari_warnings()
                                    if self.showAbChannels:

                                        for ch in reversed(abx_channels):
                                            channel_number = (
                                                marker_channel_number(
                                                    markers, ch))

                                            # read antibody image
                                            for file_path in glob.glob(
                                              f'{self.inDir}/tif/' +
                                              f'{value}.*tif'):
                                                img = single_channel_pyramid(
                                                    file_path,
                                                    channel=(
                                                        channel_number
                                                        .item()-1)
                                                        )

                                            viewer.add_image(
                                                img, rgb=False,
                                                blending='additive',
                                                colormap='green',
                                                visible=False,
                                                name=ch)

                                    # color noisy data points by
                                    # module used to redact them
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(modules[1:]),
                                        numCatagories=10,
                                        cmap='tab10',
                                        continuous=False)

                                    # reverse module order so they appear in
                                    # correct order in Napari
                                    QC_color_dict = dict(
                                        zip(modules[1:], cmap.colors))

                                    for module, color in reversed(
                                      QC_color_dict.items()):

                                        centroids = chunk[
                                            ['Y_centroid', 'X_centroid']][
                                                (chunk.index.isin(
                                                 selector.ind))
                                                & (chunk['Sample'] ==
                                                   value)
                                                & (chunk['QC_status'] ==
                                                   'noisy')
                                                & (chunk['Filter'] ==
                                                   module)]

                                        viewer.add_points(
                                            centroids, name=module,
                                            visible=False, face_color=color,
                                            edge_width=0.0, size=4.0)

                                    # read segmentation outlines, add to Napari
                                    for file_path in glob.glob(
                                      f'{self.inDir}/seg/' +
                                      f'{value}.*tif'):
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
                                        [i for i in chunk.columns
                                         if dna_moniker in i])[-1]
                                    channel_number = marker_channel_number(
                                        markers, last_dna_cycle)
                                    for file_path in glob.glob(
                                      f'{self.inDir}/tif/' +
                                      f'{value}.*tif'):
                                        dna_last = single_channel_pyramid(
                                            file_path,
                                            channel=channel_number.item() - 1)
                                        viewer.add_image(
                                            dna_last, rgb=False,
                                            blending='additive',
                                            opacity=0.5, colormap='gray',
                                            visible=False,
                                            name=f'{last_dna_cycle}: ' +
                                            f'{value}')

                                    # read first DNA, add to Napari
                                    for file_path in glob.glob(
                                      f'{self.inDir}/tif/' +
                                      f'{value}.*tif'):
                                        dna_first = single_channel_pyramid(
                                            file_path, channel=0)
                                    viewer.add_image(
                                        dna_first, rgb=False,
                                        blending='additive',
                                        opacity=0.5, colormap='gray',
                                        visible=True, name=f'{dna1}: ' +
                                        f'{value}')

                                else:
                                    print()
                                    print('Invalid sample name entered.')
                                    pass

                        #######################################################

                        sample_selector.native.setSizePolicy(
                            QtWidgets.QSizePolicy.Maximum,
                            QtWidgets.QSizePolicy.Maximum,
                        )

                        #######################################################
                        @magicgui(
                            layout='horizontal',
                            call_button='Reclass Cutoffs',
                            cleanReclass={
                                'label': 'Clean Reclass', 'max': 1.0},
                            noisyReclass={
                                'label': 'Noisy Reclass', 'max': 1.0},
                            )
                        def reclass_selector(
                            cleanReclass: float = 1.0,
                            noisyReclass: float = 1.0,
                          ):

                            return cleanReclass, noisyReclass, chunk
                        #######################################################

                        #######################################################
                        @reclass_selector.called.connect
                        def my_callback(value: str):

                            print()
                            print(
                                'Lower reclassification ' +
                                f'selection: {value[0]}')
                            print(
                                'Upper reclassification ' +
                                f'selection: {value[1]}')

                            chunk, clean, noisy = reclassify_chunk(
                                chunk=value[2],
                                clean_cutoff=value[0],
                                noisy_cutoff=value[1])

                            # build cmap
                            cmap = ListedColormap(
                                np.array([[0.91, 0.29, 0.235],
                                          [0.18, 0.16, 0.15]]))

                            sample_dict = dict(
                                zip(natsorted(
                                        chunk['Reclass'].unique()),
                                    list(range(len(chunk['Reclass']
                                         .unique())))))

                            c = [sample_dict[i] for
                                 i in chunk['Reclass']]

                            ax_reclass.cla()

                            ax_reclass.scatter(
                                chunk['emb1'], chunk['emb2'], c=c,
                                cmap=cmap, alpha=1.0, s=point_size,
                                ec='k', linewidth=0.0)

                            ax_reclass.set_title(
                                'Reclassification', fontsize=10)
                            ax_reclass.axis('equal')
                            ax_reclass.axes.xaxis.set_visible(False)
                            ax_reclass.axes.yaxis.set_visible(False)
                            ax_reclass.grid(False)

                            fig.canvas.draw()

                        #######################################################

                        reclass_selector.native.setSizePolicy(
                            QtWidgets.QSizePolicy.Maximum,
                            QtWidgets.QSizePolicy.Maximum,
                        )

                        #######################################################

                        #######################################################
                        @magicgui(
                            layout='horizontal',
                            call_button='Save'
                            )
                        def save_selector():

                            current_MCS = cluster_and_plot[0].value
                            current_cleanReclass = reclass_selector[0].value
                            current_noisyReclass = reclass_selector[1].value

                            print()
                            print(
                                'Saving current min ' +
                                f'cluster size: {current_MCS}')
                            with open(
                               os.path.join(reclass_dir, 'MCS.txt'), 'w') as f:
                                f.write(str(current_MCS))

                            print()
                            print('Saving current clean reclassification ' +
                                  f'cutoff: {current_cleanReclass}')
                            print('Saving current noisy reclassification ' +
                                  f'cutoff: {current_noisyReclass}')
                            with open(os.path.join(
                               reclass_dir, 'RECLASS_TUPLE.txt'), 'w') as f:
                                f.write(
                                    str((current_cleanReclass,
                                         current_noisyReclass)))

                            print()
                            print('Saving cluster plot')
                            if selector:
                                selector.disconnect()
                            fig.savefig(os.path.join(
                                reclass_dir,
                                f'{self.embeddingAlgorithm}_'
                                f'{current_MCS}.png'),
                                bbox_inches='tight', dpi=1000)
                            plt.close('all')

                            QTimer().singleShot(0, viewer.close)

                        #######################################################

                        save_selector.native.setSizePolicy(
                            QtWidgets.QSizePolicy.Maximum,
                            QtWidgets.QSizePolicy.Maximum,
                        )

                        #######################################################

                        cluster_layout.addWidget(sample_selector.native)

                        cluster_layout.addWidget(reclass_selector.native)

                        cluster_layout.addWidget(save_selector.native)

                        return MCS

                        #######################################################

                    cluster_and_plot.native.setSizePolicy(
                        QtWidgets.QSizePolicy.Maximum,
                        QtWidgets.QSizePolicy.Maximum,
                    )

                    #######################################################
                    @magicgui(
                        layout='horizontal',
                        call_button='Sweep Range',
                        lowerMCS={'label': 'Lower MCS'},
                        upperMCS={'label': 'Upper MCS'},
                        )
                    def sweep_MCS(
                      lowerMCS: int = 200.0,
                      upperMCS: int = 200.0,
                      ):

                        rnge = list(
                            range(lowerMCS, upperMCS + 1, 1))

                        print()
                        for i in rnge:

                            clustering = hdbscan.HDBSCAN(
                                min_cluster_size=i,
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
                                    chunk[['emb1', 'emb2']])

                            chunk[f'cluster_{self.dimensionEmbeddingQC}d'] = (
                                clustering.labels_
                                )

                            print(
                                f'min_cluster_size = {i}',
                                np.unique(clustering.labels_))
                    ##########################################################

                    sweep_MCS.native.setSizePolicy(
                        QtWidgets.QSizePolicy.Maximum,
                        QtWidgets.QSizePolicy.Maximum,
                    )

                    viewer.window.add_dock_widget(
                        cluster_and_plot, name='plot single MCS',
                        area='right')

                    viewer.window.add_dock_widget(
                        sweep_MCS, name='sweep MCS range',
                        area='right')

                    cluster_dock = viewer.window.add_dock_widget(
                        cluster_widget, name='clustering result', area='right')

                    napari.run()

                ###############################################################
                # once optimal MCS has been saved
                if os.path.exists(os.path.join(reclass_dir, 'MCS.txt')):

                    with open(
                      os.path.join(reclass_dir, 'MCS.txt'), 'r') as f:
                        final_mcs_entry = f.readlines()
                        final_mcs_entry = int(final_mcs_entry[0])

                    with open(
                      os.path.join(
                       reclass_dir, 'RECLASS_TUPLE.txt'), 'r') as f:
                        final_reclass_entry = f.readlines()
                        final_reclass_entry = (
                            final_reclass_entry[0]
                            .lstrip('([')
                            .rstrip('])')
                            .split(', '))
                        final_reclass_entry = [
                            float(i) for i in final_reclass_entry]

                    ###########################################################
                    # cluster chunk using selected MCS
                    # (not applicable to first chunk, which gets
                    # clustered during plotting above)

                    if f'cluster_{self.dimensionEmbeddingQC}d' not in chunk.columns:

                        print()
                        print(
                            'Applying saved minimum cluster size: ' +
                            f'{final_mcs_entry}')

                        clustering = hdbscan.HDBSCAN(
                            min_cluster_size=final_mcs_entry,
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
                        chunk[f'cluster_{self.dimensionEmbeddingQC}d'] = (
                            clustering.labels_
                            )

                    ###########################################################
                    # add clean and noisy data (based on final reclass_tuple)
                    # to reclass_storage_dict

                    print()
                    print(
                        'Applying saved clean reclassification ' +
                        f'cutoff: {final_reclass_entry[0]}')
                    print(
                        'Applying saved noisy reclassification ' +
                        f'cutoff: {final_reclass_entry[1]}')

                    chunk, clean, noisy = reclassify_chunk(
                        chunk=chunk,
                        clean_cutoff=final_reclass_entry[0],
                        noisy_cutoff=final_reclass_entry[1])

                    # update clean and noisy storage dataframes and save
                    reclass_storage_dict['clean'] = pd.concat(
                        [reclass_storage_dict['clean'], clean], axis=0
                        )
                    reclass_storage_dict['noisy'] = pd.concat(
                        [reclass_storage_dict['noisy'], noisy], axis=0
                        )

                    f = open(os.path.join(
                        reclass_dir, 'reclass_storage_dict.pkl'), 'wb')
                    pickle.dump(reclass_storage_dict, f)
                    f.close()

                    print()
                    print('Reclassified clean tally: ' +
                          f"{len(reclass_storage_dict['clean'])}")
                    print('Reclassified noisy tally: ' +
                          f"{len(reclass_storage_dict['noisy'])}")
                    print()

                    ###########################################################
                    # get clustermap

                    # exit program if all cells are considered ambiguous by the
                    # clustering algorithm (likely too few cells per chunk)
                    if chunk[
                        f'cluster_{self.dimensionEmbeddingQC}d'].eq(-1).all():
                        print(
                            f'WARNING: All cells in chunk {chunk_index + 1} ' +
                            'were deemed ambiguous by clustering algorithm ' +
                            '(i.e. assigned to cluster -1), exiting program. '
                            + 'Try using a larger batch size.')
                        exit()

                    clustermap_input = chunk[
                        chunk[f'cluster_{self.dimensionEmbeddingQC}d'] != -1]

                    cluster_heatmap_input = (
                        clustermap_input[
                            abx_channels +
                            [f'cluster_{self.dimensionEmbeddingQC}d']]
                        .groupby(f'cluster_{self.dimensionEmbeddingQC}d')
                        .mean()
                        )

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

                    with open(
                      os.path.join(reclass_dir, 'chunk_index.txt'), 'w') as f:
                        f.write(str(chunk_index))

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

            data = pd.concat([dropped, replaced], axis=0)

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

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

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
                data=data
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

                                neighbors_df = scatter_input.loc[
                                    list(neighbors)
                                    ]

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

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def clustering(data, self, args):

        print()
        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # drop antibody channel exclusions for clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering]

        # create N-dimensional clustering subdirectory if it hasn't already
        dim_dir = os.path.join(
            self.outDir, 'clustering', f'{self.dimensionEmbedding}d'
            )
        if not os.path.exists(dim_dir):
            os.makedirs(dim_dir)

        # recapitulate df index at the point of embedding
        data = data[~data['Sample'].isin(self.samplesToRemoveClustering)]

        # pick a random seed for reproducibility of
        # data sampling via self.fracForEmbedding
        random_state = 5

        if self.normalizeTissueCounts:
            print('Performing weighted sampling of cells from tissues...')
            print(
                'Check that resulting cell counts are similar across samples.')
            print('If not, try embedding a smaller fraction of data ' +
                  'with the fracForEmbedding configuration setting.')
            print()

            # calculate per tissue random samples weighted by cell count
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
        if os.path.exists(os.path.join(dim_dir, 'embedding.npy')):

            embedding = np.load(os.path.join(dim_dir, 'embedding.npy'))

            if embedding.shape[1] == 2:
                data['emb1'] = embedding[:, 0]
                data['emb2'] = embedding[:, 1]
                clustering_input = data[['emb1', 'emb2']]

            elif embedding.shape[1] == 3:
                data['emb1'] = embedding[:, 0]
                data['emb2'] = embedding[:, 1]
                data['emb3'] = embedding[:, 2]
                clustering_input = data[['emb1', 'emb2', 'emb3']]

        else:
            # exit program is dimensionEmbedding configuration is not 2 or 3
            if self.dimensionEmbedding not in [2, 3]:
                print()
                print('Embedding dimension must be set to 2 or 3.')
                exit()

            startTime = datetime.now()

            if self.embeddingAlgorithm == 'TSNE':
                print('Computing TSNE embedding.')
                embedding = TSNE(
                    n_components=self.dimensionEmbedding,
                    perplexity=self.perplexity,
                    early_exaggeration=self.earlyExaggeration,
                    learning_rate=self.learningRateTSNE,
                    metric=self.metric,
                    random_state=self.randomStateTSNE,
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
                    random_state=self.randomStateUMAP,
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

            np.save(
                os.path.join(dim_dir, 'embedding'), embedding
                )

            if embedding.shape[1] == 2:
                data['emb1'] = embedding[:, 0]
                data['emb2'] = embedding[:, 1]
                clustering_input = data[['emb1', 'emb2']]

            elif embedding.shape[1] == 3:
                data['emb1'] = embedding[:, 0]
                data['emb2'] = embedding[:, 1]
                data['emb3'] = embedding[:, 2]
                clustering_input = data[['emb1', 'emb2', 'emb3']]

        #######################################################################

        # define the point size for cells in the embedding
        point_size = 30000/len(data)

        # interact with plots to identify optimal minimum cluster size
        while not os.path.isfile(os.path.join(dim_dir, 'MCS.txt')):

            # initial Napari viewer without images
            viewer = napari.Viewer()

            # generate Qt widget
            cluster_widget = QtWidgets.QWidget()

            # generate vertical widget layout
            cluster_layout = QtWidgets.QVBoxLayout(cluster_widget)

            cluster_widget.setSizePolicy(
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Maximum,
                )

            ###################################################################
            @magicgui(
                layout='horizontal',
                call_button='Cluster and Plot',
                MCS={'label': 'Min Cluster Size (MCS)'},
                )
            def cluster_and_plot(MCS: int = 200.0):

                clustering = hdbscan.HDBSCAN(
                    min_cluster_size=MCS,
                    min_samples=None,
                    metric='euclidean', alpha=1.0, p=None,
                    algorithm='best', leaf_size=40,
                    memory=Memory(location=None),
                    approx_min_span_tree=True,
                    gen_min_span_tree=False,
                    core_dist_n_jobs=-1,
                    cluster_selection_method='eom',
                    allow_single_cluster=False,
                    prediction_data=False,
                    match_reference_implementation=False).fit(clustering_input)

                data[f'cluster_{self.dimensionEmbedding}d'] = (
                    clustering.labels_
                    )

                print()
                print(
                    f'min_cluster_size = {MCS}',
                    np.unique(clustering.labels_))

                if embedding.shape[1] == 2:

                    # placeholder for lasso selection
                    selector = None

                    sns.set_style('whitegrid')

                    fig = plt.figure(figsize=(8, 7))
                    matplotlib_warnings(fig)

                    gs = plt.GridSpec(2, 3, figure=fig)

                    # define axes
                    ax_cluster = fig.add_subplot(gs[0, 0])
                    ax_channel = fig.add_subplot(gs[0, 1])
                    ax_sample = fig.add_subplot(gs[0, 2])

                    ax_cluster_lbs = fig.add_subplot(gs[1, 0])
                    ax_channel_lbs = fig.add_subplot(gs[1, 1])
                    ax_sample_lbs = fig.add_subplot(gs[1, 2])

                    plt.subplots_adjust(
                        left=0.01, right=0.99, bottom=0.0,
                        top=0.95, wspace=0.0, hspace=0.0)

                    # PLOT embedding
                    for color_by in [
                      f'cluster_{self.dimensionEmbedding}d',
                      'channel', 'Sample'
                      ]:

                        highlight = 'none'

                        if color_by == f'cluster_{self.dimensionEmbedding}d':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(data[color_by].unique()),
                                numCatagories=10, cmap='tab10',
                                continuous=False)

                            # make black the first color to specify
                            # unclustered cells (cluster -1)
                            if -1 in data[color_by].unique():

                                cmap = ListedColormap(
                                    np.insert(arr=cmap.colors, obj=0,
                                              values=[0, 0, 0], axis=0))

                                # trim cmap to # unique samples
                                trim = (
                                    len(cmap.colors) - len(
                                        data[color_by].unique()))
                                cmap = ListedColormap(
                                    cmap.colors[:-trim])

                            sample_dict = dict(
                                zip(natsorted(data[color_by].unique()),
                                    list(range(len(data[color_by]
                                         .unique())))))

                            c = [sample_dict[i] for i
                                 in data[color_by]]

                            ax_cluster.cla()

                            cluster_paths = ax_cluster.scatter(
                                data['emb1'], data['emb2'], c=c, alpha=1.0,
                                s=point_size, cmap=cmap, ec='k', linewidth=0.0)

                            ax_cluster.set_title(
                                'HDBSCAN (lasso)', fontsize=10)
                            ax_cluster.set_aspect('equal')
                            ax_cluster.axes.xaxis.set_visible(False)
                            ax_cluster.axes.yaxis.set_visible(False)
                            ax_cluster.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                              natsorted(data[color_by].unique())):

                                norm_ax, hi_markers = cluster_expression(
                                    df=data, markers=abx_channels,
                                    cluster=i, num_proteins=3,
                                    clus_dim=self.dimensionEmbedding,
                                    norm_ax=self.topMarkers
                                    )

                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none',
                                           label=(
                                            f'Cluster {i}: '
                                            f'{hi_markers} {norm_ax}'),
                                           markerfacecolor=(
                                            cmap.colors[e]),
                                           markeredgecolor='none',
                                           lw=0.001, markersize=3))

                            cluster_lgd = ax_cluster_lbs.legend(
                                handles=legend_elements, prop={'size': 5},
                                loc='upper left', frameon=False)

                            ax_cluster_lbs.axis('off')

                        elif color_by == 'channel':

                            ax_channel.cla()

                            if self.colormapChannel is None:
                                c = 'gray'
                                title = 'No channel selected'
                            elif self.colormapChannel in abx_channels:
                                c = data[self.colormapChannel]
                                title = self.colormapChannel
                            else:
                                print()
                                print(
                                    'Selected channel for colormap ' +
                                    'not in markers.csv'
                                    )
                                sys.exit()

                            ax_channel.scatter(
                                data['emb1'], data['emb2'], cmap='viridis',
                                c=c, alpha=1.0, s=point_size, ec='k',
                                linewidth=0.0
                                )

                            ax_channel.set_title(title, fontsize=10)
                            ax_channel.set_aspect('equal')
                            ax_channel.axes.xaxis.set_visible(False)
                            ax_channel.axes.yaxis.set_visible(False)
                            ax_channel.grid(False)

                            ax_channel_lbs.axis('off')

                        elif color_by == 'Sample':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(data['Sample'].unique()),
                                numCatagories=10, cmap='tab10',
                                continuous=False)

                            sample_dict = dict(
                                zip(natsorted(data['Sample'].unique()),
                                    list(range(len(data['Sample']
                                         .unique())))))

                            c = [sample_dict[i] for
                                 i in data['Sample']]

                            ax_sample.cla()

                            ax_sample.scatter(
                                data['emb1'], data['emb2'], c=c, cmap=cmap,
                                alpha=1.0, s=point_size, ec='k', linewidth=0.0)

                            ax_sample.set_title('Sample', fontsize=10)
                            ax_sample.set_aspect('equal')
                            ax_sample.axes.xaxis.set_visible(False)
                            ax_sample.axes.yaxis.set_visible(False)
                            ax_sample.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                              natsorted(data['Sample'].unique())):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                sample_to_map = (
                                    data['Sample'][data['Sample'] == i]
                                    .unique()[0])

                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none', label=i,
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001, markersize=3))

                            sample_lgd = ax_sample_lbs.legend(
                                handles=legend_elements, prop={'size': 5},
                                loc='upper left', frameon=False)

                            ax_sample_lbs.axis('off')

                elif embedding.shape[1] == 3:
                    # initialize figure and add to FigureCanvas
                    # before rendering plot in 3D
                    sns.set_style('whitegrid')
                    fig = plt.figure(figsize=(7, 7))
                    matplotlib_warnings(fig)

                count = cluster_layout.count()
                # print('layout count:', count)
                # print('widget children:', pruned_widget.children())

                # remove old widgets from widget layout
                for i in range(count - 1, -1, -1):
                    item = cluster_layout.itemAt(i)
                    widget = item.widget()
                    # print('    item:', item)
                    # print('        widget:', widget)
                    if widget:
                        widget.setParent(None)

                # add updated widgets to widget layout
                cluster_canvas = FigureCanvas(fig)

                if embedding.shape[1] == 3:

                    gs = plt.GridSpec(2, 3, figure=fig)

                    # define axes
                    ax_cluster = fig.add_subplot(gs[0, 0], projection='3d')
                    ax_channel = fig.add_subplot(gs[0, 1], projection='3d')
                    ax_sample = fig.add_subplot(gs[0, 2], projection='3d')

                    ax_cluster_lbs = fig.add_subplot(gs[1, 0])
                    ax_channel_lbs = fig.add_subplot(gs[1, 1])
                    ax_sample_lbs = fig.add_subplot(gs[1, 2])

                    plt.subplots_adjust(
                        left=0.0, right=1.00, bottom=0.0,
                        top=0.94, wspace=0.0, hspace=0.0)

                    for color_by in [
                      f'cluster_{self.dimensionEmbedding}d',
                      'channel', 'Sample'
                      ]:

                        highlight = 'none'

                        if color_by == f'cluster_{self.dimensionEmbedding}d':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(data[color_by].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            # make black the first color to specify
                            # cluster outliers (i.e. cluster -1 cells)
                            cmap = ListedColormap(
                                np.insert(
                                    arr=cmap.colors, obj=0,
                                    values=[0.0, 0.0, 0.0], axis=0)
                                    )

                            # trim qualitative cmap to number of unique samples
                            trim = (
                                len(cmap.colors) - len(data[color_by].unique())
                                )
                            cmap = ListedColormap(
                                cmap.colors[:-trim]
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(data[color_by].unique()),
                                    list(range(len(data[color_by].unique()))))
                                    )

                            c = [sample_dict[i] for i in data[color_by]]

                            ax_cluster.scatter(
                                data['emb1'], data['emb2'], data['emb3'],
                                c=c, cmap=cmap, alpha=1.0, s=point_size,
                                ec='k', linewidth=0.0)

                            ax_cluster.set_title(
                                'HDBSCAN', fontsize=10)
                            ax_cluster.tick_params(labelsize=5)
                            ax_cluster.xaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})
                            ax_cluster.yaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})
                            ax_cluster.zaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})

                            legend_elements = []
                            for e, i in enumerate(
                              natsorted(data[color_by].unique())):

                                norm_ax, hi_markers = cluster_expression(
                                    df=data, markers=abx_channels,
                                    cluster=i, num_proteins=3,
                                    clus_dim=self.dimensionEmbedding,
                                    norm_ax=self.topMarkers
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
                                           lw=0.001, markersize=3))

                            cluster_lgd = ax_cluster_lbs.legend(
                                handles=legend_elements, prop={'size': 5},
                                loc='upper left', frameon=False)
                            ax_cluster_lbs.axis('off')

                        elif color_by == 'channel':

                            if self.colormapChannel is None:
                                c = 'gray'
                                title = 'No channel selected'
                            elif self.colormapChannel in abx_channels:
                                c = data[self.colormapChannel]
                                title = self.colormapChannel
                            else:
                                print()
                                print(
                                    'Selected channel for colormap ' +
                                    'not in markers.csv'
                                    )
                                sys.exit()

                            ax_channel.scatter(
                                data['emb1'], data['emb2'], data['emb3'],
                                cmap='viridis', c=c, alpha=1.0, s=point_size,
                                ec='k', linewidth=0.0)

                            ax_channel.set_title(title, fontsize=10)
                            ax_channel.tick_params(labelsize=5)
                            ax_channel.xaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})
                            ax_channel.yaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})
                            ax_channel.zaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})

                            ax_channel_lbs.axis('off')

                        elif color_by == 'Sample':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(data['Sample'].unique()),
                                numCatagories=10, cmap='tab10',
                                continuous=False)

                            sample_dict = dict(
                                zip(natsorted(data['Sample'].unique()),
                                    list(range(len(data['Sample']
                                         .unique())))))

                            c = [sample_dict[i] for
                                 i in data['Sample']]

                            ax_sample.scatter(
                                data['emb1'], data['emb2'], data['emb3'],
                                c=c, cmap=cmap, alpha=1.0, s=point_size,
                                ec='k', linewidth=0.0)

                            ax_sample.set_title(
                                'Sample', fontsize=10)
                            ax_sample.tick_params(labelsize=5)
                            ax_sample.xaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})
                            ax_sample.yaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})
                            ax_sample.zaxis._axinfo['grid'].update(
                                {'linewidth': 0.5})

                            legend_elements = []
                            for e, i in enumerate(
                              natsorted(data['Sample'].unique())):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                sample_to_map = (
                                    data['Sample'][data['Sample'] == i]
                                    .unique()[0])

                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none', label=i,
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001, markersize=3))

                            sample_lgd = ax_sample_lbs.legend(
                                handles=legend_elements, prop={'size': 5},
                                loc='upper left', frameon=False)

                            ax_sample_lbs.axis('off')

                cluster_layout.addWidget(
                    NavigationToolbar(cluster_canvas, cluster_widget))
                cluster_layout.addWidget(cluster_canvas)

                if embedding.shape[1] == 2:
                    # must call draw() before creating selector,
                    # or alpha setting doesn't work.
                    fig.canvas.draw()

                    if selector:
                        selector.disconnect()
                    selector = SelectFromCollection(
                        ax_cluster, cluster_paths)
                ###############################################################

                ###############################################################
                @magicgui(
                    layout='horizontal',
                    call_button='View Lassoed Points',
                    sample_name={'label': 'Sample Name'},
                    )
                def sample_selector(
                  sample_name: str,
                  ):

                    return sample_name
                ###############################################################

                ###############################################################
                @sample_selector.called.connect
                def my_callback(value: str):

                    print()
                    print(f'Sample selection: {value}')

                    # if cells lassoed
                    if (
                      selector.ind is None or
                      len(selector.ind) == 0
                      ):

                        print()
                        print(
                            'Cells must be lassoed in HDBSCAN plot '
                            'before sample inspection.')
                        pass

                    else:
                        # if valid sample name entered
                        if value in data['Sample'].unique():

                            # show highest expression channels
                            lasso = data.loc[
                                data.index.isin(selector.ind)
                                ].copy()

                            # assign lassoed data a dummy
                            # cluster variable and get highest
                            # expressed markers across channels
                            lasso[
                                f'cluster_{self.dimensionEmbedding}d'
                                ] = 1000

                            norm_ax, hi_markers = cluster_expression(
                                df=lasso, markers=abx_channels,
                                cluster=1000, num_proteins=3,
                                clus_dim=self.dimensionEmbedding,
                                norm_ax='channels'
                                )

                            print()
                            print(
                                'Top three expressed '
                                f'markers {hi_markers} channels'
                                )

                            # superimpose centroids of lassoed noisy cells
                            # colored by stage removed over channel images
                            print()
                            print(f'Opening sample {value} in Napari...')

                            # remove previous samples' image layers
                            for i in reversed(range(0, len(viewer.layers))):
                                viewer.layers.pop(i)

                            napari_warnings()
                            if self.showAbChannels:

                                for ch in reversed(abx_channels):
                                    channel_number = marker_channel_number(
                                        markers, ch)

                                    # read antibody image
                                    for file_path in glob.glob(
                                      f'{self.inDir}/tif/' +
                                      f'{value}.*tif'):
                                        img = single_channel_pyramid(
                                            file_path,
                                            channel=channel_number.item() - 1)

                                    viewer.add_image(
                                        img, rgb=False,
                                        blending='additive',
                                        colormap='green',
                                        visible=False,
                                        name=ch)

                            centroids = data[['Y_centroid', 'X_centroid']][
                                    (data.index.isin(selector.ind))
                                    & (data['Sample'] == value)]

                            viewer.add_points(
                                centroids, name='lassoed cells', visible=True,
                                face_color='lime', edge_width=0.0, size=4.0)

                            # read segmentation outlines, add to Napari
                            for file_path in glob.glob(
                              f'{self.inDir}/seg/{value}.*tif'):
                                seg = single_channel_pyramid(
                                    file_path, channel=0)
                            viewer.add_image(
                                seg, rgb=False,
                                blending='additive',
                                colormap='gray',
                                visible=False,
                                name='segmentation')

                            # read first DNA, add to Napari
                            for file_path in glob.glob(
                              f'{self.inDir}/tif/{value}.*tif'):
                                dna_first = single_channel_pyramid(
                                    file_path, channel=0)
                                viewer.add_image(
                                    dna_first, rgb=False, blending='additive',
                                    opacity=0.5, colormap='gray', visible=True,
                                    name=f'{dna1}: ' +
                                    f'{value}')

                        else:
                            print()
                            print('Invalid sample name entered.')
                            pass

                ###############################################################

                sample_selector.native.setSizePolicy(
                    QtWidgets.QSizePolicy.Maximum,
                    QtWidgets.QSizePolicy.Maximum,
                )

                ###############################################################
                @magicgui(
                    layout='horizontal',
                    call_button='Save'
                    )
                def save_selector():

                    current_MCS = cluster_and_plot[0].value

                    print()
                    print(
                        f'Saving current min cluster size: {current_MCS}')
                    with open(
                       os.path.join(dim_dir, 'MCS.txt'), 'w') as f:
                        f.write(str(current_MCS))

                    if embedding.shape[1] == 2:
                        print()
                        print('Saving 2D cluster plot')

                        if selector:
                            selector.disconnect()

                        fig.savefig(os.path.join(
                            dim_dir,
                            f'{self.embeddingAlgorithm}_'
                            f'{current_MCS}.png'),
                            bbox_inches='tight', dpi=1000)
                        plt.close('all')

                    elif embedding.shape[1] == 3:
                        print()
                        print('Saving 3D cluster plot')

                        def animate(i):
                            ax_cluster.view_init(elev=10., azim=i)
                            ax_channel.view_init(elev=10., azim=i)
                            ax_sample.view_init(elev=10., azim=i)
                            return fig,

                        anim = animation.FuncAnimation(
                            fig, animate, frames=360, interval=20, blit=True)

                        anim.save(os.path.join(
                            dim_dir,
                            f'{self.embeddingAlgorithm}_'
                            f'{current_MCS}.mp4'), dpi=500, fps=30,
                            extra_args=['-vcodec', 'libx264'],
                            savefig_kwargs={'bbox_inches': 'tight'})
                        plt.close('all')

                    QTimer().singleShot(0, viewer.close)

                ###############################################################

                save_selector.native.setSizePolicy(
                    QtWidgets.QSizePolicy.Maximum,
                    QtWidgets.QSizePolicy.Maximum,
                )

                ###############################################################
                if self.dimensionEmbedding == 2:
                    cluster_layout.addWidget(sample_selector.native)

                cluster_layout.addWidget(save_selector.native)

                return MCS

            ###################################################################

            cluster_and_plot.native.setSizePolicy(
                QtWidgets.QSizePolicy.Maximum,
                QtWidgets.QSizePolicy.Maximum,
            )

            ###################################################################
            @magicgui(
                layout='horizontal',
                call_button='Sweep Range',
                lowerMCS={'label': 'Lower MCS'},
                upperMCS={'label': 'Upper MCS'},
                )
            def sweep_MCS(
              lowerMCS: int = 200.0,
              upperMCS: int = 200.0,
              ):

                rnge = list(
                    range(lowerMCS, upperMCS + 1, 1))

                print()
                for i in rnge:

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=i,
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
                            clustering_input)

                    data[f'cluster_{self.dimensionEmbedding}d'] = (
                        clustering.labels_
                        )

                    print(
                        f'min_cluster_size = {i}',
                        np.unique(clustering.labels_))
            ###################################################################

            sweep_MCS.native.setSizePolicy(
                QtWidgets.QSizePolicy.Maximum,
                QtWidgets.QSizePolicy.Maximum,
            )

            mcs_dock = viewer.window.add_dock_widget(
                cluster_and_plot, name='plot single MCS', area='right')

            sweep_dock = viewer.window.add_dock_widget(
                sweep_MCS, name='sweep MCS range', area='right')

            cluster_dock = viewer.window.add_dock_widget(
                cluster_widget, name='clustering result', area='right')

            napari.run()

        #######################################################################
        # apply final MCS and return data from clustering module

        with open(
          os.path.join(dim_dir, 'MCS.txt'), 'r') as f:
            final_mcs_entry = f.readlines()
            final_mcs_entry = int(final_mcs_entry[0])

        print()
        print(
            'Applying saved minimum cluster size: ' +
            f'{final_mcs_entry}')

        clustering = hdbscan.HDBSCAN(
            min_cluster_size=final_mcs_entry,
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
                clustering_input)
        data[f'cluster_{self.dimensionEmbedding}d'] = clustering.labels_

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def clustermap(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # create clustering dimension directory if it hasn't already
        dim_dir = os.path.join(
            self.outDir, 'clustering', f'{self.dimensionEmbedding}d'
            )
        if not os.path.exists(dim_dir):
            os.makedirs(dim_dir)

        # drop antibody channel exclusions for clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering]

        # drop unclustered cells before plotting clustermap
        clustermap_input = data[
            data[f'cluster_{self.dimensionEmbedding}d'] != -1]

        # compute mean antibody signals for clusters
        clustermap_input = (
            clustermap_input[abx_channels +
                             [f'cluster_{self.dimensionEmbedding}d']]
            .groupby(f'cluster_{self.dimensionEmbedding}d').mean())

        sns.set(font_scale=0.8)
        for name, axis in zip(['channels', 'clusters'], [0, 1]):

            g = sns.clustermap(
                clustermap_input, cmap='viridis', standard_scale=axis,
                square=False, yticklabels=1, linewidth=0.1, cbar=True)

            g.fig.suptitle(f'Normalized across {name}', y=0.995, fontsize=10)
            g.fig.set_size_inches(6.0, 6.0)

            plt.savefig(os.path.join(
                    dim_dir, f'clustermap_norm_{name}.pdf'),
                    bbox_inches='tight')

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data

    @module
    def setContrast(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # create contrast directory if it hasn't already
        contrast_dir = os.path.join(self.outDir, 'contrast')
        if not os.path.exists(contrast_dir):
            os.makedirs(contrast_dir)

        napari_warnings()

        if self.viewSample in data['Sample'].unique():
            # loop over antibody channels and add them to Napari viewer
            for e, ch in enumerate(reversed(abx_channels)):

                channel_number = marker_channel_number(markers, ch)

                # read antibody image
                for file_path in glob.glob(
                  f'{self.inDir}/tif/{self.viewSample}.*tif'):
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
              f'{self.inDir}/tif/{self.viewSample}.*tif'):
                dna = single_channel_pyramid(file_path, channel=0)

            viewer.add_image(
                dna, rgb=False, blending='additive', colormap='gray',
                name=f'{dna1}: {self.viewSample}')

            # apply previously defined contrast limits if they exist
            if os.path.exists(
              os.path.join(contrast_dir, 'contrast_limits.yml')):

                print()
                print('Reading existing contrast settings.')

                contrast_limits = yaml.safe_load(
                    open(f'{contrast_dir}/contrast_limits.yml'))

                viewer.layers[f'{dna1}: {self.viewSample}'].contrast_limits = (
                    contrast_limits[dna1][0], contrast_limits[dna1][1])

                for ch in reversed(abx_channels):
                    viewer.layers[ch].contrast_limits = (
                        contrast_limits[ch][0], contrast_limits[ch][1])

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

            data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

            print()
            print()
            return data

        else:
            print()
            print(
                'Sample name for image contrast adjustment not valid. '
                'Please assign a valid sample name to the "viewSample" '
                'setting in config.yml.')
            exit()

    @module
    def frequencyStats(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # prepare input data for computing statistics
        stats_input = data[
            ['Sample', 'Replicate', f'cluster_{self.dimensionEmbedding}d']][
            data[f'cluster_{self.dimensionEmbedding}d'] >= 0]

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
                    self.outDir, 'clustering', f'{self.dimensionEmbedding}d',
                    'frequency_stats', f'{test}_v_{control}'
                    )
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
                for cluster, group in stats_input.groupby(
                  f'cluster_{self.dimensionEmbedding}d'
                  ):

                    print(
                        f'Calculating log2({test}/{control})'
                        f' of mean cell density for cluster {str(cluster)}.')

                    group = (
                        group.groupby(
                            ['Sample', 'Replicate',
                             f'cluster_{self.dimensionEmbedding}d'])
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

                    group[f'cluster_{self.dimensionEmbedding}d'] = cluster

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
                                (stats_input[
                                    f'cluster_{self.dimensionEmbedding}d'] ==
                                 self.denominatorCluster)])
                            for i in group['Sample']]

                    # compute density of cells per sample
                    group['density'] = group['count']/group['tissue_count']

                    # append group data to catplot_input
                    catplot_input = pd.concat([catplot_input, group], axis=0)

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
                    columns=[f'cluster_{self.dimensionEmbedding}d',
                             'ratio', 'dif', 'pval']).sort_values(
                        by=f'cluster_{self.dimensionEmbedding}d')

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
                  significant[f'cluster_{self.dimensionEmbedding}d'],
                  significant[stat],
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
                     'ns' if i not in
                     significant[
                        f'cluster_{self.dimensionEmbedding}d'].unique() else
                     significant[stat][
                        significant[
                            f'cluster_{self.dimensionEmbedding}d']
                        == i].values[0]
                     for i in catplot_input[
                        f'cluster_{self.dimensionEmbedding}d']]

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
                    by=[f'cluster_{self.dimensionEmbedding}d',
                        'status', 'density'],
                    ascending=[True, False, True], inplace=True)

                catplot_input[f'cluster_{self.dimensionEmbedding}d'] = (
                    catplot_input[
                        f'cluster_{self.dimensionEmbedding}d'].astype(str) +
                    f'; {stat} = ' + catplot_input[stat].astype(str)
                    )

                catplot_input[f'cluster_{self.dimensionEmbedding}d'] = [
                    i.split(f'; {stat} = ns')[0] for i in
                    catplot_input[f'cluster_{self.dimensionEmbedding}d']]

                sns.set(font_scale=0.4)
                g = sns.catplot(
                    x='status', y='density',
                    hue=catplot_input['Sample'],
                    col=f'cluster_{self.dimensionEmbedding}d', col_wrap=6,
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

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        return data

    @module
    def curateThumbnails(data, self, args):

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=os.path.join(self.inDir, 'markers.csv'),
            markers_to_exclude=self.markersToExclude,
            data=data
            )

        # drop antibody channel exclusions for clustering
        abx_channels = [
            i for i in abx_channels
            if i not in self.channelExclusionsClustering]

        # drop unclustered cells from data
        data = data[data[f'cluster_{self.dimensionEmbedding}d'] != -1]

        # create thumbnails directory
        thumbnails_dir = os.path.join(
            self.outDir, 'clustering', f'{self.dimensionEmbedding}d',
            'thumbnails'
            )
        if not os.path.exists(thumbnails_dir):
            os.makedirs(thumbnails_dir)

        # create zarr directory
        zarr_dir = os.path.join(thumbnails_dir, 'zarrs')
        if not os.path.exists(zarr_dir):
            os.makedirs(zarr_dir)

        # read image contrast settings
        contrast_dir = os.path.join(self.outDir, 'contrast')
        if os.path.exists(f'{contrast_dir}/contrast_limits.yml'):
            contrast_limits = yaml.safe_load(
                open(f'{contrast_dir}/contrast_limits.yml'))

        #######################################################################

        # read the indices of clusters that have already been run
        if os.path.exists(os.path.join(
          thumbnails_dir, 'completed_clusters.txt')):

            with open(
              os.path.join(
               thumbnails_dir, 'completed_clusters.txt'), 'r') as f:
                completed_clusters = f.readlines()
                completed_clusters = (
                    completed_clusters[0].lstrip('[').rstrip(']').split(', ')
                    )
                completed_clusters = set([int(i) for i in completed_clusters])

            total_clusters = set(
                data[f'cluster_{self.dimensionEmbedding}d'].unique()
                )

            clusters_to_run = natsorted(
                total_clusters.difference(completed_clusters)
                )

            # convert completed_clusters from a set to a list to append
            # to while looping over clusters
            completed_clusters = list(completed_clusters)

            print(f'Clusters to run: {len(clusters_to_run)}')
            print()
        else:
            # create a list of clusters to run
            completed_clusters = []
            clusters_to_run = natsorted(
                data[f'cluster_{self.dimensionEmbedding}d'].unique()
                )
            print(f'Clusters to run: {len(clusters_to_run)}')
            print()

        #######################################################################

        for cluster in clusters_to_run:

            # create dataframe to collect thumbnail images and their metadata
            long_table = pd.DataFrame()

            # identify top expressed markers
            norm_ax, hi_markers = cluster_expression(
                df=data, markers=abx_channels,
                cluster=cluster, num_proteins=3,
                clus_dim=self.dimensionEmbedding,
                norm_ax=self.topMarkersThumbnails
                )

            # combine DNA1 with top expressed markers
            markers_to_show = [dna1] + hi_markers

            # get marker channel numbers from markers.csv
            channel_nums = []
            for marker in markers_to_show:
                channel_number = marker_channel_number(markers, marker)
                channel_nums.append(str(channel_number.item()))

            # create marker LUT
            color_dict = {}
            for i, j in zip(
              markers_to_show,
              [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0),
               (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]):
                color_dict[i] = j

            for sample in [i for i in data['Sample'].unique()]:

                # isolate data for current cluster and sample
                cellcutter_input = data[
                    (data[f'cluster_{self.dimensionEmbedding}d']
                     == cluster) &
                    (data['Sample'] == sample)
                    ]

                # randomly select example cells
                if len(cellcutter_input) > self.numThumbnails:

                    cellcutter_input = cellcutter_input.sample(
                        n=self.numThumbnails, random_state=1
                        )

                # generate cellcutter zarr file if it doesn't already exist
                if not os.path.exists(
                  os.path.join(
                    zarr_dir, f'clus{cluster}_{sample}'
                    f'_win{self.windowSize}.zarr'
                    )
                  ):

                    # write cellcutter_input to disk
                    cellcutter_input.to_csv(
                        os.path.join(thumbnails_dir, 'csv_data.csv'),
                        index=False
                        )

                    print(
                        f'Cutting cluster {cluster} cells '
                        f'showing {markers_to_show} '
                        f'normalized across {norm_ax} '
                        f'from sample {sample}...'
                        )
                    print()

                    # run cellcutter on sample image
                    run(
                        ["cut_cells", "--window-size", f"{self.windowSize}",
                         "--cells-per-chunk", "200", "--cache-size", "57711",
                         f"{self.inDir}/tif/{sample}.ome.tif",
                         f"{self.inDir}/mask/{sample}.tif",
                         f"{thumbnails_dir}/csv_data.csv",
                         f"{zarr_dir}/clus{cluster}_{sample}"
                         f"_win{self.windowSize}.zarr",
                         "--channels"] + channel_nums
                         )
                    print()

                # read multi-channel zarr file created by cellcutter
                z_path_img = os.path.join(
                    zarr_dir, f'clus{cluster}_{sample}'
                    f'_win{self.windowSize}.zarr'
                    )
                z_img = zarr.open(z_path_img, mode='r')

                if self.segOutlines:

                    if not os.path.exists(
                      os.path.join(
                        zarr_dir, f'clus{cluster}_{sample}'
                        f'_win{self.windowSize}_seg.zarr'
                        )
                      ):

                        # write cellcutter_input to disk
                        cellcutter_input.to_csv(
                            os.path.join(thumbnails_dir, 'csv_data.csv'),
                            index=False
                            )

                        # run cellcutter on segmentation outlines image
                        run(
                            ["cut_cells", "--window-size",
                             f"{self.windowSize}",
                             "--cells-per-chunk", "200",
                             "--cache-size", "57711",
                             f"{self.inDir}/seg/{sample}.ome.tif",
                             f"{self.inDir}/mask/{sample}.tif",
                             f"{thumbnails_dir}/csv_data.csv",
                             f"{zarr_dir}/clus{cluster}_{sample}"
                             f"_win{self.windowSize}_seg.zarr",
                             "--channels", "1"]
                             )
                        print()

                    # read segmentation outlines zarr file
                    # created by cellcutter
                    z_path_seg = os.path.join(
                        zarr_dir, f'clus{cluster}_{sample}'
                        f'_win{self.windowSize}_seg.zarr'
                        )
                    z_seg = zarr.open(z_path_seg, mode='r')

                if os.path.exists(
                  os.path.join(thumbnails_dir, 'csv_data.csv')
                  ):
                    # remove cellcutter_input file after cells have been cut
                    os.remove(os.path.join(thumbnails_dir, 'csv_data.csv'))

                # create composite thumbnail images
                for cell in range(z_img.shape[1]):

                    # create blank image with same x, y dimensions as thumbnail
                    blank_img = np.zeros((z_img.shape[2], z_img.shape[3]))

                    # add centroid point at the center of the image
                    blank_img[
                        int(z_img.shape[2]/2):int(z_img.shape[2]/2)+1,
                        int(z_img.shape[3]/2):int(z_img.shape[3]/2)+1
                        ] = 1

                    # convert blank image to rgb and colorize
                    blank_img = gray2rgb(blank_img)

                    # loop over markers to show
                    for ch, marker in enumerate(markers_to_show):

                        # slice marker channel of current thumbnail
                        slice = img_as_float(z_img[ch, cell, :, :])

                        # apply image contrast settings
                        slice -= (contrast_limits[marker][0]/65535)
                        slice /= (
                            (contrast_limits[marker][1]/65535)
                            - (contrast_limits[marker][0]/65535))
                        slice = np.clip(slice, 0, 1)

                        # convert channel slice to RGB, colorize,
                        # and add to blank image
                        slice = gray2rgb(slice)
                        slice = (slice * color_dict[marker])
                        blank_img += slice

                    if self.segOutlines:

                        # get segmentation thumbnails
                        seg = img_as_float(z_seg[0, cell, :, :])

                        # convert segmentation thumbnail to RGB
                        # and add to blank image
                        seg = gray2rgb(seg)
                        blank_img += seg

                    # append merged thumbnail image to long_table
                    append_df = pd.DataFrame.from_dict(
                        {'sample': sample,
                         'example': int(cell+1),
                         'image': blank_img}, orient='index'
                         ).T
                    long_table = pd.concat(
                        [long_table, append_df], axis=0, ignore_index=True
                        )

            # plot facet grid of thumbnails for current cluster
            if not long_table.empty:
                fig, ax = plt.subplots()

                g = sns.FacetGrid(
                    long_table, row='sample', col='example',
                    sharex=False, sharey=False,
                    gridspec_kws={'hspace': 0.1, 'wspace': 0.1})

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
                    fontweight='bold', size=6)

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

            # update completed clusters list
            completed_clusters.append(cluster)

            # overwrite completed_clusters.txt file
            with open(
              os.path.join(
               thumbnails_dir, 'completed_clusters.txt'), 'w') as f:
                f.write(str(completed_clusters))
            print()

        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        return data
