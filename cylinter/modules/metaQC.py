import os
import sys
import pickle
import logging

import numpy as np
import pandas as pd

import math
from natsort import natsorted
from sklearn.preprocessing import MinMaxScaler

from datetime import datetime

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from qtpy.QtCore import QTimer
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

from umap import UMAP
from sklearn.manifold import TSNE
import hdbscan
from joblib import Memory

from magicgui import magicgui

import napari

from ..utils import (
    input_check, read_markers, matplotlib_warnings, categorical_cmap, SelectFromCollection, 
    cluster_expression, napari_notification, marker_channel_number, single_channel_pyramid,
    get_filepath, reorganize_dfcolumns
)

logger = logging.getLogger(__name__)


def metaQC(data, self, args):

    print()
    
    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
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
                    module_dict[module_idx + 1][1].index)].copy()

            # if noisy data exists, add a QC stage column
            if not noise.empty:
                noise.loc[:, 'Filter'] = (
                    module_dict[module_idx + 1][0])

            # append to module_dict
            module_dict[module_idx + 1].append(noise)

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
            logger.info(
                'No data were filtered during prior QC steps. '
                'Returning unfiltered data without reclassification.'
            )
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
            with open(os.path.join(reclass_dir, 'chunk_index.txt'), 'r') as f:
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
            with open(os.path.join(reclass_dir, 'chunk_index.txt'), 'w') as f:
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
        if len(QCData) < (batch_size) * 2:
            num_chunks = 1
            chunks = np.array_split(QCData, num_chunks)
        else:
            num_chunks = math.ceil(len(QCData) / batch_size)
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
        if os.path.exists(os.path.join(reclass_dir, 'reclass_storage_dict.pkl')):
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

            logger.info(f'Clustering: {chunk_index + 1} of {len(chunks)} data chunks')

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

                embedding = np.load(os.path.join(chunk_dir, 'embedding.npy'))
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
                    logger.info('Computing TSNE embedding...')
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
                    logger.info('Computing UMAP embedding...')
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

                logger.info(f'Embedding completed in {str(datetime.now() - startTime)}')

                np.save(os.path.join(chunk_dir, 'embedding'), embedding)

                chunk['emb1'] = embedding[:, 0]
                chunk['emb2'] = embedding[:, 1]

            # define the point size for cells in the embedding
            point_size = 50000 / len(chunk)

            def reclassify_chunk(chunk, clean_cutoff, noisy_cutoff):

                clean = pd.DataFrame()
                noisy = pd.DataFrame()
                for name, cluster in chunk.groupby(f'cluster_{self.dimensionEmbeddingQC}d'):
                    if name != -1:
                        
                        # if a cluster contains >= n% clean data,
                        # reclassify all clustering cells as noisy
                        if (
                           (len(cluster[cluster[
                            'QC_status'] == 'clean']) / len(cluster)) >=
                           clean_cutoff):
                            clean = pd.concat([clean, cluster], axis=0)
                        
                        # elif a cluster contains >= n% noisy data,
                        # reclassify all clustering cells as clean
                        elif (
                             (len(cluster[cluster[
                              'QC_status'] == 'noisy']) / len(cluster)) >=
                                noisy_cutoff):
                            noisy = pd.concat([noisy, cluster], axis=0)
                        else:
                            noisy = pd.concat(
                                [noisy, cluster[cluster['QC_status'] == 'noisy']], axis=0
                            )
                            clean = pd.concat(
                                [clean, cluster[cluster['QC_status'] == 'clean']], axis=0
                            )

                # consider -1 cells from clean data
                # to be "noisy"
                clean_outliers = chunk[
                    (chunk[f'cluster_{self.dimensionEmbeddingQC}d'] == -1)
                    &
                    (chunk['QC_status'] == 'clean')].copy()
                noisy = pd.concat([noisy, clean_outliers], axis=0)

                # consider -1 cells from noisy data
                # to be "noisy"
                noisy_outliers = chunk[
                    (chunk[f'cluster_{self.dimensionEmbeddingQC}d'] == -1)
                    &
                    (chunk['QC_status'] == 'noisy')].copy()
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
                viewer = napari.Viewer(title='CyLinter')

                # generate Qt widget
                cluster_widget = QtWidgets.QWidget()

                # generate vertical widget layout
                cluster_layout = QtWidgets.QVBoxLayout(cluster_widget)

                cluster_widget.setSizePolicy(
                    QtWidgets.QSizePolicy.Fixed,
                    QtWidgets.QSizePolicy.Maximum,
                )

                ###########################################################
                @magicgui(
                    layout='horizontal',
                    call_button='Cluster and Plot',
                    MCS={'label': 'Min Cluster Size (MCS)', 'step': 1},
                )
                def cluster_and_plot(MCS: int = 200.0):

                    # placeholder for lasso selection
                    selector = None

                    sns.set_style('whitegrid')

                    fig = plt.figure(figsize=(3.5, 8))
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

                    chunk[f'cluster_{self.dimensionEmbeddingQC}d'] = clustering.labels_

                    # scatter point selection tool assumes a
                    # sorted index, but the index of QCdata is
                    # shuffled to acheive a mix of clean and
                    # noisy data per chunk
                    chunk.sort_index(inplace=True)

                    print()
                    logger.info(f'min_cluster_size = {MCS} {np.unique(clustering.labels_)}')
                   
                    # PLOT embedding
                    for color_by in [
                        f'cluster_{self.dimensionEmbeddingQC}d',
                            'QC_status', 'Reclass', 'annotation']:

                        highlight = 'none'

                        if color_by == f'cluster_{self.dimensionEmbeddingQC}d':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(chunk[color_by].unique()),
                                numCatagories=10, cmap='tab10', continuous=False
                            )

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
                                zip(natsorted(chunk[color_by].unique()),
                                    list(range(len(chunk[color_by].unique())))))

                            c = [sample_dict[i] for i in chunk[color_by]]

                            ax_cluster.cla()

                            cluster_paths = ax_cluster.scatter(
                                chunk['emb1'], chunk['emb2'], c=c, alpha=1.0, s=point_size,
                                cmap=cmap, ec='k', linewidth=0.0
                            )

                            ax_cluster.set_title(
                                'HDBSCAN (lasso)', fontsize=7)
                            ax_cluster.axis('equal')
                            ax_cluster.axes.xaxis.set_visible(False)
                            ax_cluster.axes.yaxis.set_visible(False)
                            ax_cluster.grid(False)

                            legend_elements = []
                            for e, i in enumerate(natsorted(chunk[color_by].unique())):

                                norm_ax, hi_markers = cluster_expression(
                                    df=chunk, markers=abx_channels,
                                    cluster=i, num_proteins=3,
                                    clus_dim=self.dimensionEmbeddingQC,
                                    norm_ax=self.topMarkersQC
                                )

                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none',
                                           label=(f'Cluster {i}: 'f'{hi_markers} by {norm_ax}'),
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor='none',
                                           lw=0.001, markersize=2)
                                )

                            ax_cluster_lbs.legend(
                                handles=legend_elements, prop={'size': 3}, loc='upper left',
                                frameon=False
                            )

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

                            c = [sample_dict[i] for i in chunk['QC_status']]

                            ax_status.cla()

                            ax_status.scatter(
                                chunk['emb1'], chunk['emb2'], c=c, cmap=cmap, alpha=1.0,
                                s=point_size, ec='k', linewidth=0.0
                            )

                            ax_status.set_title('QC Status', fontsize=7)
                            ax_status.axis('equal')
                            ax_status.axes.xaxis.set_visible(False)
                            ax_status.axes.yaxis.set_visible(False)
                            ax_status.grid(False)

                            legend_elements = []
                            for e, i in enumerate(natsorted(chunk['QC_status'].unique())):

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
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001,
                                           markersize=2)
                                )

                            ax_status_lbs.legend(
                                handles=legend_elements, prop={'size': 5}, loc='upper left',
                                frameon=False
                            )

                            ax_status_lbs.axis('off')

                        elif color_by == 'Reclass':

                            # build cmap
                            cmap = ListedColormap(
                                np.array([[0.91, 0.29, 0.235],
                                          [0.18, 0.16, 0.15]]))

                            sample_dict = dict(
                                zip(natsorted(chunk['QC_status'].unique()),
                                    list(range(len(chunk['QC_status'].unique())))))

                            c = [sample_dict[i] for i in chunk['QC_status']]

                            ax_reclass.cla()

                            ax_reclass.scatter(
                                chunk['emb1'], chunk['emb2'], c=c, cmap=cmap, alpha=1.0,
                                s=point_size, ec='k', linewidth=0.0
                            )

                            ax_reclass.set_title('Reclassification', fontsize=7)
                            ax_reclass.axis('equal')
                            ax_reclass.axes.xaxis.set_visible(False)
                            ax_reclass.axes.yaxis.set_visible(False)
                            ax_reclass.grid(False)

                            legend_elements = []
                            for e, i in enumerate(natsorted(chunk['QC_status'].unique())):

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
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001,
                                           markersize=2)
                                )

                            ax_reclass_lbs.legend(
                                handles=legend_elements, prop={'size': 5}, loc='upper left',
                                frameon=False
                            )

                            ax_reclass_lbs.axis('off')

                        elif color_by == 'annotation':

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(
                                    chunk[self.colormapAnnotationQC].unique()),
                                numCatagories=10, cmap='tab10', continuous=False
                            )

                            sample_dict = dict(
                                zip(natsorted(chunk[self.colormapAnnotationQC].unique()),
                                    list(range(len(chunk[self.colormapAnnotationQC].unique())))))

                            c = [sample_dict[i] for i in chunk[self.colormapAnnotationQC]]

                            ax_sample.cla()

                            ax_sample.scatter(
                                chunk['emb1'], chunk['emb2'], c=c, cmap=cmap, alpha=1.0,
                                s=point_size, ec='k', linewidth=0.0
                            )

                            ax_sample.set_title(self.colormapAnnotationQC, fontsize=7)
                            ax_sample.axis('equal')
                            ax_sample.axes.xaxis.set_visible(False)
                            ax_sample.axes.yaxis.set_visible(False)
                            ax_sample.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                    natsorted(chunk[self.colormapAnnotationQC].unique())):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none',
                                           label=i,
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001, markersize=2)
                                )

                            ax_sample_lbs.legend(
                                handles=legend_elements, prop={'size': 3}, loc='upper left',
                                frameon=False
                            )

                            ax_sample_lbs.axis('off')

                    count = cluster_layout.count()
                    # logger.info('widget children:', pruned_widget.children())

                    # remove old widgets from widget layout
                    for i in range(count - 1, -1, -1):
                        item = cluster_layout.itemAt(i)
                        widget = item.widget()
                        # logger.info('    item:', item)
                        # logger.info('        widget:', widget)
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
                    def sample_selector(sample_name: str):

                        return sample_name
                    #######################################################

                    #######################################################
                    @sample_selector.called.connect
                    def sample_selector_callback(value: str):

                        print()
                        logger.info(f'Sample selection: {value}')

                        # if cells lassoed
                        if not (selector.ind is not None and len(selector.ind) >= 1):

                            print()
                            logger.info(
                                'Cells must be lassoed in HDBSCAN plot before sample inspection.'
                            )
                            pass

                        else:
                            # if valid sample name entered
                            if value in chunk['Sample'].unique():

                                # show highest expression channels
                                lasso = chunk.loc[chunk.index.isin(selector.ind)].copy()

                                # assign lassoed data a dummy
                                # cluster variable and get highest
                                # expressed markers across channels
                                lasso[f'cluster_{self.dimensionEmbeddingQC}d'] = 1000

                                norm_ax, hi_markers = cluster_expression(
                                    df=lasso, markers=abx_channels, cluster=1000,
                                    num_proteins=3, clus_dim=self.dimensionEmbeddingQC,
                                    norm_ax='channels'
                                )

                                print()
                                logger.info(f'Top three expressed markers {hi_markers} channels')

                                # superimpose centroids of lassoed noisy cells colored 
                                # by stage removed over channel images
                                print()
                                logger.info(f'Opening sample {value} in Napari...')

                                # remove previous samples' image layers
                                for i in reversed(range(0, len(viewer.layers))):
                                    viewer.layers.pop(i)

                                if self.showAbChannels:

                                    for ch in reversed(abx_channels):
                                        channel_number = (
                                            marker_channel_number(
                                                markers, ch))

                                        # read antibody image
                                        file_path = get_filepath(self, check, value, 'TIF')
                                        img, min, max = single_channel_pyramid(
                                            file_path, channel=channel_number
                                        )
                                        viewer.add_image(
                                            img, rgb=False, blending='additive',
                                            colormap='green', visible=False, name=ch,
                                            contrast_limits=(min, max)
                                        )

                                # color noisy data points by
                                # module used to redact them
                                cmap = categorical_cmap(
                                    numUniqueSamples=len(modules[1:]), numCatagories=10,
                                    cmap='tab10', continuous=False
                                )

                                # reverse module order so they appear in
                                # correct order in Napari
                                QC_color_dict = dict(zip(modules[1:], cmap.colors))

                                for module, color in reversed(QC_color_dict.items()):

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
                                file_path = get_filepath(self, check, value, 'SEG')
                                seg, min, max = single_channel_pyramid(file_path, channel=0)
                                viewer.add_image(
                                    seg, rgb=False, blending='additive', colormap='gray',
                                    visible=False, name='segmentation', 
                                    contrast_limits=(min, max)
                                )

                                # read last DNA, add to Napari
                                last_dna_cycle = natsorted(
                                    [i for i in chunk.columns
                                     if dna_moniker in i])[-1]
                                channel_number = marker_channel_number(
                                    markers, last_dna_cycle)
                                
                                file_path = get_filepath(self, check, value, 'TIF')
                                dna_last, min, max = single_channel_pyramid(
                                    file_path, channel=channel_number
                                )
                                viewer.add_image(
                                    dna_last, rgb=False, blending='additive',
                                    opacity=0.5, colormap='gray', visible=False,
                                    name=f'{last_dna_cycle}: {value}',
                                    contrast_limits=(min, max)
                                )

                                # read first DNA, add to Napari
                                file_path = get_filepath(self, check, value, 'TIF')
                                dna_first, min, max = single_channel_pyramid(file_path, channel=0)
                                viewer.add_image(
                                    dna_first, rgb=False, blending='additive',
                                    opacity=0.5, colormap='gray', visible=True,
                                    name=f'{dna1}: {value}', contrast_limits=(min, max)
                                )

                            else:
                                print()
                                logger.info('Invalid sample name entered.')
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
                    def reclass_selector_callback(value: str):

                        print()
                        napari_notification(f'Lower reclassification selection: {value[0]}')
                        napari_notification(f'Upper reclassification selection: {value[1]}')

                        chunk, clean, noisy = reclassify_chunk(
                            chunk=value[2], clean_cutoff=value[0], noisy_cutoff=value[1]
                        )

                        # build cmap
                        cmap = ListedColormap(
                            np.array([[0.91, 0.29, 0.235],
                                      [0.18, 0.16, 0.15]]))

                        sample_dict = dict(
                            zip(natsorted(chunk['Reclass'].unique()),
                                list(range(len(chunk['Reclass'].unique())))))

                        c = [sample_dict[i] for i in chunk['Reclass']]

                        ax_reclass.cla()

                        ax_reclass.scatter(
                            chunk['emb1'], chunk['emb2'], c=c, cmap=cmap, alpha=1.0,
                            s=point_size, ec='k', linewidth=0.0
                        )

                        ax_reclass.set_title('Reclassification', fontsize=7)
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
                    @magicgui(layout='horizontal', call_button='Save')
                    def save_selector():

                        current_MCS = cluster_and_plot[0].value
                        current_cleanReclass = reclass_selector[0].value
                        current_noisyReclass = reclass_selector[1].value

                        print()
                        logger.info(f'Saving current min cluster size: {current_MCS}')
                        with open(os.path.join(reclass_dir, 'MCS.txt'), 'w') as f:
                            f.write(str(current_MCS))

                        print()
                        logger.info(
                            'Saving current clean reclassification ' +
                            f'cutoff: {current_cleanReclass}'
                        )
                        logger.info(
                            'Saving current noisy reclassification ' +
                            f'cutoff: {current_noisyReclass}'
                        )
                        with open(os.path.join(reclass_dir, 'RECLASS_TUPLE.txt'), 'w') as f:
                            f.write(str((current_cleanReclass, current_noisyReclass)))

                        print()
                        logger.info('Saving cluster plot')
                        if selector:
                            selector.disconnect()
                        fig.savefig(os.path.join(
                            reclass_dir,
                            f'{self.embeddingAlgorithm}_'
                            f'{current_MCS}.png'),
                            bbox_inches='tight', dpi=1000)

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
                    QtWidgets.QSizePolicy.Fixed,
                    QtWidgets.QSizePolicy.Maximum,
                )

                #######################################################
                @magicgui(
                    layout='horizontal',
                    call_button='Sweep Range',
                    lowerMCS={'label': 'Lower MCS', 'step': 1},
                    upperMCS={'label': 'Upper MCS', 'step': 1},
                )
                def sweep_MCS(lowerMCS: int = 200.0, upperMCS: int = 200.0,):

                    rnge = list(range(lowerMCS, upperMCS + 1, 1))

                    print()
                    for i in rnge:

                        clustering = hdbscan.HDBSCAN(
                            min_cluster_size=i,
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

                        chunk[f'cluster_{self.dimensionEmbeddingQC}d'] = clustering.labels_

                        logger.info(f'min_cluster_size = {i} {np.unique(clustering.labels_)}')
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

                viewer.window.add_dock_widget(
                    cluster_widget, name='clustering result', area='right')

                viewer.scale_bar.visible = True
                viewer.scale_bar.unit = 'um'

                napari.run()

            ###############################################################
            # once optimal MCS has been saved
            if os.path.exists(os.path.join(reclass_dir, 'MCS.txt')):

                with open(os.path.join(reclass_dir, 'MCS.txt'), 'r') as f:
                    final_mcs_entry = f.readlines()
                    final_mcs_entry = int(final_mcs_entry[0])

                with open(os.path.join(reclass_dir, 'RECLASS_TUPLE.txt'), 'r') as f:
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
                    logger.info(f'Applying saved minimum cluster size: {final_mcs_entry}')

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=final_mcs_entry,
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
                        match_reference_implementation=False).fit(chunk[['emb1', 'emb2']])
                    chunk[f'cluster_{self.dimensionEmbeddingQC}d'] = clustering.labels_

                ###########################################################
                # add clean and noisy data (based on final reclass_tuple)
                # to reclass_storage_dict

                print()
                logger.info(
                    'Applying saved clean reclassification ' +
                    f'cutoff: {final_reclass_entry[0]}')
                logger.info(
                    'Applying saved noisy reclassification ' +
                    f'cutoff: {final_reclass_entry[1]}')

                chunk, clean, noisy = reclassify_chunk(
                    chunk=chunk, clean_cutoff=final_reclass_entry[0],
                    noisy_cutoff=final_reclass_entry[1]
                )

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
                logger.info(f"Reclassified clean tally: {len(reclass_storage_dict['clean'])}")
                logger.info(f"Reclassified noisy tally: {len(reclass_storage_dict['noisy'])}")
                print()

                ###########################################################
                # get clustermap

                # exit program if all cells are considered ambiguous by the
                # clustering algorithm (likely too few cells per chunk)
                if chunk[f'cluster_{self.dimensionEmbeddingQC}d'].eq(-1).all():
                    logger.warning(
                        f'Aborting; all cells in chunk {chunk_index + 1} ' +
                        'were deemed ambiguous by clustering algorithm ' +
                        '(i.e. assigned to cluster -1). '
                        + 'Try using larger batch size.')
                    sys.exit()

                clustermap_input = chunk[chunk[f'cluster_{self.dimensionEmbeddingQC}d'] != -1]

                cluster_heatmap_input = (
                    clustermap_input[
                        abx_channels +
                        [f'cluster_{self.dimensionEmbeddingQC}d']]
                    .groupby(f'cluster_{self.dimensionEmbeddingQC}d')
                    .mean()
                )

                sns.set(font_scale=0.8)
                for name, axis in zip(['within', 'across'], [0, 1]):
                    sns.clustermap(
                        cluster_heatmap_input, cmap='viridis', standard_scale=axis,
                        square=False, yticklabels=1, linewidth=0.1, cbar=True
                    )

                    plt.gcf().set_size_inches(8.0, 8.0)

                    plt.savefig(
                        os.path.join(chunk_dir, f'clustermap_{name}.pdf'),
                        bbox_inches='tight'
                    )
                    plt.close('all')

                ###########################################################
                # increment chunk_index
                chunk_index = chunk_index + 1

                with open(os.path.join(reclass_dir, 'chunk_index.txt'), 'w') as f:
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
                    pre_qc['CellID'].map(str) + '_' + pre_qc['Sample']
                )

        # create explicit global labels for
        # cleaned data (before reclassifiction)
        post_qc = module_dict[
            [i for i in module_dict.keys()][-1]][1].copy()
        post_qc['handle'] = (
            post_qc['CellID'].map(str) + '_' + post_qc['Sample']
        )

        # get raw values of cells in post_qc data
        cleaned_raw = pre_qc[pre_qc['handle'].isin(post_qc['handle'])]

        # convert clean data in predominantly noisy clusters to noisy
        # to yield final clean data
        drop = reclass_storage_dict['noisy'][
            reclass_storage_dict['noisy']['QC_status'] == 'clean'].copy()
        drop['handle'] = (
            drop['CellID'].map(str) + '_' + drop['Sample']
        )
        dropped = cleaned_raw[~cleaned_raw['handle'].isin(drop['handle'])]

        # convert noisy data in predominantly clean clusters to clean
        # to yield final replace data
        replace = reclass_storage_dict['clean'][
            reclass_storage_dict['clean']['QC_status'] == 'noisy'].copy()
        replace['handle'] = (
            replace['CellID'].map(str) + '_' + replace['Sample']
        )
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
            rescaled_data = scaler.transform(channel_data.values.reshape(-1, 1))

            rescaled_data = pd.DataFrame(
                data=rescaled_data, index=channel_data.index
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
            (qc_dict[qc_keys[i]] - qc_dict[qc_keys[i + 1]]) / qc_dict[qc_keys[0]]) * 100

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
        sizes, explode=explode, labels=labels, autopct='%1.2f%%', shadow=False, startangle=90
    )

    for i in range(len(texts)):
        texts[i].set_fontsize(4)
        autotexts[i].set_fontsize(4)

    ax1.axis('equal')

    plt.savefig(os.path.join(reclass_dir, 'censored_by_stage.pdf'), bbox_inches='tight')
    plt.close('all')

    # drop indexing handle from data before returning
    if self.metaQC:
        data.drop('handle', axis=1, inplace=True)

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
