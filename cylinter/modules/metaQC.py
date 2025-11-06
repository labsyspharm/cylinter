import os
import re
import sys
import glob
import yaml
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
from matplotlib.colors import TwoSlopeNorm
from matplotlib.colors import ListedColormap
from matplotlib.backends.qt_compat import QtWidgets

from qtpy.QtCore import QTimer

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

import hdbscan
from umap import UMAP
from sklearn.manifold import TSNE

from joblib import Memory

from magicgui import magicgui

import napari

from ..utils import (
    input_check, read_markers, matplotlib_warnings, categorical_cmap,
    SelectFromCollection, cluster_expression, marker_channel_number,
    single_channel_pyramid, get_filepath, reorganize_dfcolumns,
    sort_qc_report
)

logger = logging.getLogger(__name__)


def metaQC(data, self, args):

    print()

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # drop antibody channel exclusions for metaQC clustering
    abx_channels = [
        i for i in abx_channels
        if i not in self.channelExclusionsClusteringQC]

    # create metaQC directory if it hasn't already
    reclass_dir = os.path.join(
        self.outDir, 'metaQC')
    if not os.path.exists(reclass_dir):
        os.makedirs(reclass_dir)

    # read chunk index if it exists
    progress_path = os.path.join(reclass_dir, 'progress.yml')
    try:
        progress = yaml.safe_load(open(progress_path))
    except FileNotFoundError:
        pass

    # read QC report
    report_path = os.path.join(self.outDir, 'cylinter_report.yml')
    try:
        qc_report = yaml.safe_load(open(report_path))
        reload_report = False
        if qc_report is None:
            qc_report = {}
            reload_report = True
        if 'metaQC' not in qc_report or qc_report['metaQC'] is None:
            qc_report['metaQC'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(
                qc_report, module='metaQC', order=None
            )
            f = open(report_path, 'w')
            yaml.dump(
                qc_report_sorted, f, sort_keys=False, allow_unicode=False
            )
            qc_report = yaml.safe_load(open(report_path))
    except:
        logger.info(
            'Aborting; QC report missing from CyLinter output directory. '
            'Re-start pipeline from aggregateData module to '
            'initialize QC report.'
        )
        sys.exit()

    # specify the names of modules in the pipeline that perform
    # data redaction prior to the metaQC module
    modules = ['aggregateData', 'selectROIs', 'intensityFilter',
               'areaFilter', 'cycleCorrelation', 'pruneOutliers']

    # build a dictionary of clean data returned by each module
    module_dict = {}
    for module_idx, module in enumerate(modules):
        data = pd.read_parquet(
            os.path.join(
                self.outDir, f'checkpoints/{module}.parquet'))
        module_dict[module_idx] = [module, data]
    
    ###########################################################################

    if self.metaQC:
        
        # specify channels on which to perform metaQC clustering
        abx_channels_dna = [self.counterstainChannel] + abx_channels
        
        #######################################################################
        # Build QCData: QCData is a combination of clean (retained) 
        # and noisy (dropped) data
        
        # if QCData exists, read the parquet file
        if os.path.exists(os.path.join(reclass_dir, 'QCData.parquet')):
            
            QCData = pd.read_parquet(
                os.path.join(reclass_dir, 'QCData.parquet')
            )
        
        else:  # generate QCData
             
            # loop over modules to get data redacted by each module (noisy)
            for module_idx in [i for i in module_dict.keys()][:-1]:
                
                # isolate redacted data
                noise = module_dict[module_idx][1][
                    ~module_dict[module_idx][1].index.isin(
                        module_dict[module_idx + 1][1].index)].copy()

                # if noisy data exists, add a QC stage column
                if not noise.empty:
                    noise.loc[:, 'Filter'] = (
                        module_dict[module_idx + 1][0])

                # append as new value to module_dict key
                module_dict[module_idx + 1].append(noise)

            # combine noisy data from all modules into a single dataframe
            if self.delintMode:
                # if negative ROI selection, include gated cells as noisy data
                noisyData = pd.DataFrame()
                for module_idx, module in enumerate(modules):
                    if len(module_dict[module_idx]) == 3:
                        # number of dictionary value enteries should be 3:
                        # module name, clean data, noisy data, if noisy data
                        noisyData = pd.concat(
                            [noisyData, module_dict[module_idx][2]], axis=0
                        )
            else:
                # if postive ROI selection, don't count gated cells as noisy
                noisyData = pd.DataFrame()
                for module_idx, module in enumerate(modules):
                    if len(module_dict[module_idx]) == 3:
                        if not module == 'selectROIs':
                            noisyData = pd.concat(
                                [noisyData, module_dict[module_idx][2]], axis=0
                            )
        
            if noisyData.empty:
                logger.info(
                    'No data were filtered during prior QC stages. '
                    'Returning unfiltered data without reclassification.'
                )
                return data

            # create QC_status column for combined noisy data
            noisyData.loc[:, 'QC_status'] = 'noisy'

            # get raw version (untransformed, not rescaled) of data
            # from fully-redacted (clean) dataframe
            first_module_idx = modules.index(modules[0])
            last_module_idx = modules.index(modules[-1])

            cleanDataRescaled = module_dict[last_module_idx][1]

            cleanDataRaw = module_dict[first_module_idx][1][
                module_dict[first_module_idx][1].index.isin(
                    cleanDataRescaled.index)].copy()

            # create QC_status column for selected clean data
            cleanDataRaw.loc[:, 'QC_status'] = 'clean'

            # append noisyData to cleanDataRaw, row-wise
            QCData = pd.concat([cleanDataRaw, noisyData], axis=0)

            # shuffle QCData row order to randomize cells for clustering chunks
            QCData = QCData.sample(frac=1.0, random_state=5)
            QCData.reset_index(drop=True, inplace=True)

            # log-transform QCData
            QCData.loc[:, abx_channels_dna] += 0.001
            QCData.loc[:, abx_channels_dna] = np.log10(
                QCData[abx_channels_dna]
            )

            # rescale channel signal intensities (0-1) per channel per sample
            for channel in abx_channels_dna:
                for sample in natsorted(data['Sample'].unique()):
                    
                    sample_channel_data = QCData[
                        QCData['Sample'] == sample][channel]

                    scaler = (
                        MinMaxScaler(feature_range=(0, 1), copy=True)
                        .fit(sample_channel_data.values.reshape(-1, 1))
                    )
                    rescaled_data = scaler.transform(
                        sample_channel_data.values.reshape(-1, 1)
                    )

                    rescaled_data = pd.DataFrame(
                        data=rescaled_data,
                        index=sample_channel_data.index,
                    ).rename(columns={0: channel})

                    QCData.update(rescaled_data)

            # save QCData
            QCData.to_parquet(
                os.path.join(reclass_dir, 'QCData.parquet'), index=False
            )

        #######################################################################
        # handle re-starts

        def check_report_structure(report_path):
            
            expected_structure = {
                'MCS': 1,
                'reclassThresholds': [1.0, 1.0]
            }
            
            try:
                with open(report_path, 'r') as file:
                    qc_report = yaml.safe_load(file)
                    metaQC = qc_report.get('metaQC')
                
                # check if metaQC section is present 
                if metaQC is None:
                    return False
                
                # check that metaQC keys match the expected structure
                if set(metaQC.keys()) != set(expected_structure.keys()):
                    return False
                
                # check that the value types in metaQC are as expected
                for key, value in expected_structure.items():
                    if key == 'reclassThresholds':
                        if not (isinstance(metaQC.get(key), list) and 
                                len(metaQC.get(key)) == 2 and
                                all(isinstance(item, float) for 
                                item in metaQC.get(key))):
                            return False

                    elif not isinstance(metaQC.get(key), type(value)):
                        return False
                
                return True  # all checks pass
            
            except FileNotFoundError:
                return False

        if not check_report_structure(report_path):

            # initialize chunk index counter at 0 and store
            chunk_index = 0
            progress = {}
            progress['chunk_index'] = chunk_index
            f = open(progress_path, 'w')
            yaml.dump(progress, f, sort_keys=False, allow_unicode=False)
            
            parquet = glob.glob(os.path.join(reclass_dir, 'reclass_*.parquet'))
            for i in parquet:
                try:
                    os.remove(i)
                except FileNotFoundError:
                    pass

            try:
                del qc_report['metaQC']['MCS']
            except:
                pass
            try:
                del qc_report['metaQC']['reclassThresholds']
            except:
                pass
            
            # dump updated qc_report to YAML file
            qc_report_sorted = sort_qc_report(
                qc_report, module='metaQC', order=None
            )
            f = open(report_path, 'w')
            yaml.dump(
                qc_report_sorted, f, sort_keys=False, allow_unicode=False
            )
        
        else:
            try:
                chunk_index = progress['chunk_index']
                
                if chunk_index is None:
                    chunk_index = 0
                    progress['chunk_index'] = chunk_index
                    f = open(progress_path, 'w')
                    yaml.dump(
                        progress, f, sort_keys=False, allow_unicode=False
                    )
                    
                    parquet = glob.glob(
                        os.path.join(reclass_dir, 'reclass_*.parquet')
                    )
                    for i in parquet:
                        try:
                            os.remove(i)
                        except FileNotFoundError:
                            pass

                    try:
                        del qc_report['metaQC']['MCS']
                    except:
                        pass
                    try:
                        del qc_report['metaQC']['reclassThresholds']
                    except:
                        pass
            
                    # sort and dump updated qc_report to YAML file
                    qc_report_sorted = sort_qc_report(
                        qc_report, module='metaQC', order=None
                    )
                    f = open(report_path, 'w')
                    yaml.dump(
                        qc_report_sorted, f, sort_keys=False,
                        allow_unicode=False
                    )
            
            except (NameError, TypeError, KeyError):

                # initialize chunk index counter at 0 and store
                chunk_index = 0
                progress = {}
                progress['chunk_index'] = chunk_index
                f = open(progress_path, 'w')
                yaml.dump(progress, f, sort_keys=False, allow_unicode=False)
                
                parquet = glob.glob(
                    os.path.join(reclass_dir, 'reclass_*.parquet')
                )
                for i in parquet:
                    try:
                        os.remove(i)
                    except FileNotFoundError:
                        pass
                
                try:
                    del qc_report['metaQC']['MCS']
                except:
                    pass
                try:
                    del qc_report['metaQC']['reclassThresholds']
                except:
                    pass
                
                # dump updated qc_report to YAML file
                qc_report_sorted = sort_qc_report(
                    qc_report, module='metaQC', order=None
                )
                f = open(report_path, 'w')
                yaml.dump(
                    qc_report_sorted, f, sort_keys=False, allow_unicode=False
                )
        
        #######################################################################
        # initialize reclassification storage dataframe
        def check_for_reclass_dataframe():
            
            parquet = glob.glob(os.path.join(reclass_dir, 'reclass_*.parquet'))

            if len(parquet) == 0:
                
                # initialize reclassified data storage dataframe for 
                # reclassified clean/noisy data
                reclass = pd.DataFrame(
                    columns=['Reclass', 'CellID', 'Sample', 'QC_status', 
                             'Filter', 'ChunkIDX']
                )

                # initialize chunk index counter at 0 and store
                chunk_index = 0
                progress = {}
                progress['chunk_index'] = chunk_index
                f = open(progress_path, 'w')
                yaml.dump(progress, f, sort_keys=False, allow_unicode=False)
                
                try:
                    del qc_report['metaQC']['MCS']
                except:
                    pass
                try:
                    del qc_report['metaQC']['reclassThresholds']
                except:
                    pass

                return chunk_index, reclass, None

            if len(parquet) == 1:
                
                # reclassified data storage exists, open it
                reclass = pd.read_parquet(parquet[0])
                chunk_index = reclass['ChunkIDX'].max() + 1
                filename = os.path.basename(parquet[0])

                return chunk_index, reclass, filename
            
            elif len(parquet) > 1:
                logger.info(
                    'Aborting; more than one reclass parquet file exists. '
                    'To avoid ambiguity, only one should be in metaQC '
                    'output directory.'
                )
                sys.exit()

        chunk_index, reclass, filename = check_for_reclass_dataframe()

        if filename is not None:
            
            existing_percent_data = float(filename.split('_')[1])
            current_percent_data = self.percentDataPerChunk
            
            if not isinstance(current_percent_data, float):
                logger.info(
                    'Aborting; CyLinter configuration "percentDataPerChunk" '
                    'not in floating point format. Adjust and '
                    're-run metaQC module'
                )
                sys.exit()
            
            if existing_percent_data != current_percent_data:
                logger.info(
                    'Aborting; Chunk size of current reclass parquet is '
                    f'{existing_percent_data*100}%, but "percentDataPerChunk" '
                    'parameter in CyLinter config file is '
                    f'{current_percent_data*100}%. '
                    'Update CyLinter config to match or remove '
                    'reclass parquet. Then re-run metaQC module.'
                )
                sys.exit()

        #######################################################################
        # chunk QCData

        # specify the number of cells in each clustering batch (i.e. chunk)
        # (this limits clustering time and memory pressure)
        batch_size = len(QCData) * self.percentDataPerChunk

        # ensure minimun batch size equals the batch_size variable
        # otherwise cluster all data at once
        if len(QCData) < (batch_size) * 2:
            num_chunks = 1
            chunks = np.array_split(QCData, num_chunks)
        else:
            num_chunks = math.ceil(len(QCData) / batch_size)
            chunks = np.array_split(QCData, num_chunks)
        
        #######################################################################
        # loop over QCData chunks

        for chunk in chunks[chunk_index:]:
            if isinstance(chunk, pd.DataFrame):
                
                print()
                logger.info(
                    f'Working on data chunk {chunk_index + 1} of {len(chunks)}'
                )

                # make directory for current chunk if it hasn't already
                chunk_dir = os.path.join(
                    reclass_dir, f'chunk_{str(chunk_index + 1)}'
                )
                if not os.path.exists(chunk_dir):
                    os.makedirs(chunk_dir)

                # if embedding for current chunk has already been computed
                # apply emb1 and emb2 to chunk dataframe
                if os.path.exists(os.path.join(chunk_dir, 'embedding.npy')):

                    # recapitulate chunk index at the point of embedding
                    chunk = chunk[~chunk['Sample'].isin(
                        self.samplesToRemoveClusteringQC)]
                    chunk = chunk.sample(
                        frac=1.0, random_state=5)

                    logger.info(
                        f'Embedding for chunk {chunk_index + 1} of '
                        f'{len(chunks)} already exists; '
                        'appending embedding values to chunk.'
                    )

                    embedding = np.load(
                        os.path.join(chunk_dir, 'embedding.npy')
                    )
                    chunk['emb1'] = embedding[:, 0]
                    chunk['emb2'] = embedding[:, 1]

                else:
                    # compute embedding for chunk
                    startTime = datetime.now()

                    chunk = chunk[~chunk['Sample'].isin(
                        self.samplesToRemoveClusteringQC)]
                    chunk = chunk.sample(
                        frac=1.0, random_state=5)

                    print(chunk[abx_channels_dna])

                    if self.embeddingAlgorithmQC == 'TSNE':
                        logger.info('Computing TSNE embedding...')
                        embedding = TSNE(
                            n_components=2,
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
                            n_components=2,
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

                    logger.info(
                        'Embedding completed in '
                        f'{str(datetime.now() - startTime)}'
                    )

                    np.save(os.path.join(chunk_dir, 'embedding'), embedding)

                    chunk['emb1'] = embedding[:, 0]
                    chunk['emb2'] = embedding[:, 1]

                # define the point size for cells in the embedding
                point_size = 50000 / len(chunk)

                def reclassify_chunk(chunk, clean_cutoff, noisy_cutoff):

                    clean = pd.DataFrame()
                    noisy = pd.DataFrame()
                    for name, cluster in chunk.groupby('cluster_2d'):
                        if name != -1:
                            
                            # if a cluster contains >= n% clean data,
                            # reclassify all clustering cells as clean
                            if (
                               (len(cluster[cluster[
                                'QC_status'] == 'clean']) / len(cluster)) >=
                               clean_cutoff):
                                clean = pd.concat([clean, cluster], axis=0)
                            
                            # elif a cluster contains >= n% noisy data,
                            # reclassify all clustering cells as noisy
                            elif (
                                 (len(cluster[cluster[
                                  'QC_status'] == 'noisy']) / len(cluster)) >=
                                    noisy_cutoff):
                                noisy = pd.concat([noisy, cluster], axis=0)
                            
                            # else use original QC_status assignments
                            else:
                                clean = pd.concat(
                                    [clean, 
                                     cluster[cluster['QC_status'] == 'clean']],
                                    axis=0
                                )
                                noisy = pd.concat(
                                    [noisy, 
                                     cluster[cluster['QC_status'] == 'noisy']],
                                    axis=0
                                )

                    # consider -1 cells from clean data as noisy
                    clean_outliers = chunk[
                        (chunk['cluster_2d'] == -1) &
                        (chunk['QC_status'] == 'clean')].copy()
                    noisy = pd.concat([noisy, clean_outliers], axis=0)

                    # consider -1 cells from noisy data to be noisy
                    noisy_outliers = chunk[
                        (chunk['cluster_2d'] == -1) &
                        (chunk['QC_status'] == 'noisy')].copy()
                    noisy = pd.concat([noisy, noisy_outliers], axis=0)

                    clean['Reclass'] = 'clean'
                    noisy['Reclass'] = 'noisy'
                    
                    chunk.loc[
                        chunk.index.isin(clean.index), 'Reclass'] = 'clean'
                    chunk.loc[
                        chunk.index.isin(noisy.index), 'Reclass'] = 'noisy'

                    return chunk, clean, noisy

                # interact with plots to identify optimal min cluster size
                while not check_report_structure(report_path):
                        
                    if chunk_index == 0:
                    
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

                        #######################################################
                        @magicgui(
                            layout='horizontal',
                            call_button='Cluster and Plot',
                            MCS={'label': 'Min Cluster Size (MCS)', 'step': 1},
                        )
                        def cluster_and_plot(MCS: int = 200.0):

                            # placeholder for lasso selection
                            selector = None

                            sns.set_style('whitegrid')

                            fig = plt.figure(figsize=(3.5, 3.5))
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

                            chunk['cluster_2d'] = clustering.labels_

                            # scatter point selection tool assumes a
                            # sorted index, but the index of QCdata is
                            # shuffled to acheive a mix of clean and
                            # noisy data per chunk
                            chunk.sort_index(inplace=True)

                            print()
                            logger.info(
                                f'min_cluster_size = {MCS} '
                                f'{np.unique(clustering.labels_)}'
                            )
                           
                            # PLOT embedding
                            for color_by in ['cluster_2d', 'QC_status', 
                                             'Reclass', 'annotation']:

                                highlight = 'none'

                                if color_by == 'cluster_2d':

                                    # build cmap
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(
                                            chunk[color_by].unique()),
                                        numCatagories=10, cmap='tab10',
                                        continuous=False
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
                                        zip(natsorted(
                                            chunk[color_by].unique()),
                                            list(range(len(
                                                chunk[color_by].unique())))))

                                    c = [sample_dict[i] for i 
                                         in chunk[color_by]]

                                    ax_cluster.cla()

                                    cluster_paths = ax_cluster.scatter(
                                        chunk['emb1'], chunk['emb2'], c=c, 
                                        alpha=1.0, s=point_size,
                                        cmap=cmap, ec='k', linewidth=0.0
                                    )

                                    ax_cluster.set_title(
                                        'HDBSCAN (lasso)', fontsize=9)
                                    ax_cluster.axis('equal')
                                    ax_cluster.axes.xaxis.set_visible(False)
                                    ax_cluster.axes.yaxis.set_visible(False)
                                    ax_cluster.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(
                                         natsorted(chunk[color_by].unique())):

                                        hi_markers = cluster_expression(
                                            df=chunk, markers=abx_channels,
                                            cluster=i, num_proteins=3,
                                            clus_dim=2
                                        )

                                        legend_elements.append(
                                            Line2D([0], [0], marker='o',
                                                   color='none',
                                                   label=f'C{i}: {hi_markers}',
                                                   markerfacecolor=cmap.colors[
                                                    e],
                                                   markeredgecolor='none',
                                                   lw=0.001, markersize=6)
                                        )

                                    ax_cluster_lbs.legend(
                                        handles=legend_elements,
                                        prop={'size': 6}, loc='upper left',
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

                                    c = [sample_dict[i] for i in 
                                         chunk['QC_status']]

                                    ax_status.cla()

                                    ax_status.scatter(
                                        chunk['emb1'], chunk['emb2'], c=c, 
                                        cmap=cmap, alpha=1.0,
                                        s=point_size, ec='k', linewidth=0.0
                                    )

                                    ax_status.set_title(
                                        'QC Status', fontsize=9
                                    )
                                    ax_status.axis('equal')
                                    ax_status.axes.xaxis.set_visible(False)
                                    ax_status.axes.yaxis.set_visible(False)
                                    ax_status.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(natsorted(
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
                                                   markerfacecolor=cmap.colors[
                                                    e],
                                                   markeredgecolor=(
                                                    markeredgecolor),
                                                   lw=0.001,
                                                   markersize=8)
                                        )

                                    ax_status_lbs.legend(
                                        handles=legend_elements, 
                                        prop={'size': 8}, loc='upper left',
                                        frameon=False
                                    )

                                    ax_status_lbs.axis('off')

                                elif color_by == 'Reclass':

                                    # build cmap
                                    cmap = ListedColormap(
                                        np.array([[0.91, 0.29, 0.235],
                                                  [0.18, 0.16, 0.15]]))

                                    sample_dict = dict(
                                        zip(natsorted(
                                            chunk['QC_status'].unique()),
                                            list(range(len(chunk[
                                                'QC_status'].unique())))))

                                    c = [sample_dict[i] for i in
                                         chunk['QC_status']]

                                    ax_reclass.cla()

                                    ax_reclass.scatter(
                                        chunk['emb1'], chunk['emb2'], c=c,
                                        cmap=cmap, alpha=1.0,
                                        s=point_size, ec='k', linewidth=0.0
                                    )

                                    ax_reclass.set_title(
                                        'Reclassification', fontsize=9
                                    )
                                    ax_reclass.axis('equal')
                                    ax_reclass.axes.xaxis.set_visible(False)
                                    ax_reclass.axes.yaxis.set_visible(False)
                                    ax_reclass.grid(False)

                                    legend_elements = []
                                    for e, i in enumerate(natsorted(
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
                                                   markerfacecolor=cmap.colors[
                                                    e],
                                                   markeredgecolor=(
                                                    markeredgecolor),
                                                   lw=0.001,
                                                   markersize=8)
                                        )

                                    ax_reclass_lbs.legend(
                                        handles=legend_elements, 
                                        prop={'size': 8}, loc='upper left',
                                        frameon=False
                                    )

                                    ax_reclass_lbs.axis('off')

                                elif color_by == 'annotation':

                                    # build cmap
                                    cmap = categorical_cmap(
                                        numUniqueSamples=len(
                                            chunk[self.colormapAnnotationQC
                                                  ].unique()),
                                        numCatagories=10, cmap='tab10',
                                        continuous=False
                                    )

                                    sample_dict = dict(
                                        zip(natsorted(chunk[self.colormapAnnotationQC].unique()),
                                            list(range(len(chunk[self.colormapAnnotationQC].unique()))))
                                    )

                                    c = [sample_dict[i] for i in 
                                         chunk[self.colormapAnnotationQC]]

                                    ax_sample.cla()

                                    ax_sample.scatter(
                                        chunk['emb1'], chunk['emb2'], c=c,
                                        cmap=cmap, alpha=1.0,
                                        s=point_size, ec='k', linewidth=0.0
                                    )

                                    ax_sample.set_title(
                                        self.colormapAnnotationQC, fontsize=9
                                    )
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
                                                   markerfacecolor=cmap.colors[
                                                    e],
                                                   markeredgecolor=(
                                                    markeredgecolor),
                                                   lw=0.001, markersize=8)
                                        )

                                    ax_sample_lbs.legend(
                                        handles=legend_elements, 
                                        prop={'size': 8}, loc='upper left',
                                        frameon=False
                                    )

                                    ax_sample_lbs.axis('off')

                            count = cluster_layout.count()
                            # logger.info(
                            #     'widget children:', pruned_widget.children()
                            # )

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
                                NavigationToolbar(
                                    cluster_canvas, cluster_widget)
                            )
                            cluster_layout.addWidget(cluster_canvas)

                            # must call draw() before creating selector,
                            # or alpha setting doesn't work.
                            fig.canvas.draw()

                            if selector:
                                selector.disconnect()
                            selector = SelectFromCollection(
                                ax_cluster, cluster_paths)

                            ###################################################
                            @magicgui(
                                layout='horizontal',
                                call_button='View Lassoed Points',
                                sample_name={
                                 'choices': list(natsorted(
                                        data['Sample'].unique())), 
                                 'label': 'Sample'},
                            )
                            def sample_selector(sample_name: str):

                                return sample_name
                            ###################################################

                            ###################################################
                            @sample_selector.called.connect
                            def sample_selector_callback(value: str):

                                print()
                                logger.info(f'Sample selection: {value}')

                                # if cells lassoed
                                if not (selector.ind is not None and 
                                        len(selector.ind) >= 1):

                                    print()
                                    napari.utils.notifications.show_warning(
                                        'Cells must be lassoed in HDBSCAN '
                                        'plot before sample inspection.'
                                    )
                                    pass

                                else:
                                    # if valid sample name entered
                                    if value in chunk['Sample'].unique():

                                        # show highest expression channels
                                        lasso = chunk.loc[
                                            chunk.index.isin(
                                                selector.ind)].copy()

                                        # assign lassoed data a dummy
                                        # cluster variable and get highest
                                        # expressed markers across channels
                                        lasso['cluster_2d'] = 1000

                                        hi_markers = cluster_expression(
                                            df=lasso, markers=abx_channels,
                                            cluster=1000, num_proteins=3,
                                            clus_dim=2
                                        )

                                        print()
                                        logger.info(
                                            'Top three expressed markers '
                                            f'{hi_markers} '
                                            '(raw expression values, '
                                            'not z-scores'
                                        )

                                        # superimpose centroids of lassoed 
                                        # noisy cells colored by stage removed
                                        # over channel images

                                        # remove previous samples' image layers
                                        for i in reversed(
                                             range(0, len(viewer.layers))):
                                            viewer.layers.pop(i)

                                        if self.showAbChannels:

                                            for ch in reversed(abx_channels):
                                                channel_number = (
                                                    marker_channel_number(
                                                        self, markers, ch))

                                                # read antibody image
                                                file_path = get_filepath(
                                                    self, check, value, 'TIF'
                                                )
                                                img, min, max = single_channel_pyramid(
                                                    file_path,
                                                    channel=channel_number
                                                )
                                                viewer.add_image(
                                                    img, rgb=False, 
                                                    blending='additive',
                                                    colormap='green', 
                                                    visible=False, name=ch,
                                                    contrast_limits=(min, max)
                                                )
                                        
                                        # color noisy data points by
                                        # module used to redact them
                                        cmap = categorical_cmap(
                                            numUniqueSamples=len(modules[1:]),
                                            numCatagories=10, cmap='tab10',
                                            continuous=False
                                        )

                                        # reverse module order so they appear
                                        # in correct order in Napari
                                        QC_color_dict = dict(
                                            zip(modules[1:], cmap.colors)
                                        )

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
                                                visible=True, face_color=color,
                                                border_width=0.0, size=4.0)

                                        # read segmentation outlines,
                                        # add to Napari
                                        file_path = get_filepath(
                                            self, check, value, 'SEG'
                                        )
                                        seg, min, max = single_channel_pyramid(
                                            file_path, channel=0
                                        )
                                        viewer.add_image(
                                            seg, rgb=False,
                                            blending='additive',
                                            colormap='gray',
                                            visible=False, name='segmentation', 
                                            contrast_limits=(min, max)
                                        )
                                        
                                        # get ordered list of DNA cycles
                                        dna_moniker = (
                                            str(re.search(r'[^\W\d]+',
                                                self.counterstainChannel).group())
                                        )
                                        dna_cycles = natsorted(
                                            data.columns[
                                                data.columns.str.contains(
                                                    dna_moniker)]
                                        )
                                        
                                        # read last DNA, add to Napari
                                        last_dna = dna_cycles[-1]
                                        channel_number = marker_channel_number(
                                            self, markers, last_dna
                                        )
                                        file_path = get_filepath(
                                            self, check, value, 'TIF'
                                        )
                                        dna_last, min, max = single_channel_pyramid(
                                            file_path, channel=channel_number
                                        )
                                        viewer.add_image(
                                            dna_last, rgb=False,
                                            blending='additive',
                                            opacity=0.5, colormap='gray',
                                            visible=False, name=last_dna,
                                            contrast_limits=(min, max)
                                        )

                                        # read first DNA, add to Napari
                                        first_dna = dna_cycles[0]
                                        channel_number = marker_channel_number(
                                            self, markers, first_dna
                                        )
                                        file_path = get_filepath(
                                            self, check, value, 'TIF'
                                        )
                                        dna_first, min, max = single_channel_pyramid(
                                            file_path, channel=channel_number
                                        )
                                        viewer.add_image(
                                            dna_first, rgb=False,
                                            blending='additive',
                                            opacity=0.5, colormap='gray',
                                            visible=True, name=first_dna,
                                            contrast_limits=(min, max)
                                        )

                                        # apply previously defined contrast 
                                        # limits if they exist
                                        try:
                                            viewer.layers[
                                                f'{self.counterstainChannel}'].contrast_limits = (
                                                qc_report['setContrast'][self.counterstainChannel]
                                                [0],
                                                qc_report['setContrast'][self.counterstainChannel]
                                                [1]
                                            )

                                            for ch in reversed(abx_channels):
                                                viewer.layers[ch].contrast_limits = (
                                                    qc_report['setContrast'][ch][0],
                                                    qc_report['setContrast'][ch][1]
                                                ) 
                                            logger.info(
                                                'Existing contrast '
                                                'settings applied.'
                                            )
                                        except:
                                            pass

                                        print()
                                        napari.utils.notifications.show_info(
                                            f'Viewing sample {value}.'
                                        )

                                    else:
                                        print()
                                        napari.utils.notifications.show_warning(
                                            'Sample name not in filtered data.'
                                        )
                                        pass

                            ###################################################

                            sample_selector.native.setSizePolicy(
                                QtWidgets.QSizePolicy.Maximum,
                                QtWidgets.QSizePolicy.Maximum,
                            )

                            ###################################################
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
                            ###################################################

                            ###################################################
                            @reclass_selector.called.connect
                            def reclass_selector_callback(value: str):

                                print()
                                napari.utils.notifications.show_info(
                                    f'reclass cutoffs: clean={round(value[0], 2)}, ' 
                                    f'noisy={round(value[1], 2)}'
                                )

                                chunk, clean, noisy = reclassify_chunk(
                                    chunk=value[2], clean_cutoff=value[0],
                                    noisy_cutoff=value[1]
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
                                    chunk['emb1'], chunk['emb2'], c=c,
                                    cmap=cmap, alpha=1.0,
                                    s=point_size, ec='k', linewidth=0.0
                                )

                                ax_reclass.set_title(
                                    'Reclassification', fontsize=7
                                )
                                ax_reclass.axis('equal')
                                ax_reclass.axes.xaxis.set_visible(False)
                                ax_reclass.axes.yaxis.set_visible(False)
                                ax_reclass.grid(False)

                                fig.canvas.draw()

                            ###################################################

                            reclass_selector.native.setSizePolicy(
                                QtWidgets.QSizePolicy.Maximum,
                                QtWidgets.QSizePolicy.Maximum,
                            )

                            ###################################################

                            ###################################################
                            @magicgui(layout='horizontal', call_button='Save')
                            def save_selector(reclass):

                                current_MCS = cluster_and_plot[0].value
                                current_cleanReclass = reclass_selector[0].value
                                current_noisyReclass = reclass_selector[1].value
                                
                                # save reclass storage
                                reclass.to_parquet(
                                    os.path.join(reclass_dir, 
                                                 f'reclass_{self.percentDataPerChunk}'
                                                 f'_{current_MCS}'
                                                 f'_{current_cleanReclass}'
                                                 f'_{current_noisyReclass}.parquet')
                                )

                                # write current MCS to QC report 
                                logger.info(f'Saving current min cluster size: {current_MCS}')
                                qc_report['metaQC']['MCS'] = current_MCS

                                print()
                                # write current reclass cutoffs to QC report
                                logger.info(
                                    'Saving current clean reclassification ' +
                                    f'cutoff: {current_cleanReclass}'
                                )
                                logger.info(
                                    'Saving current noisy reclassification ' +
                                    f'cutoff: {current_noisyReclass}'
                                )
                                qc_report['metaQC']['reclassThresholds'] = [
                                    current_cleanReclass, current_noisyReclass
                                ]
                                
                                qc_report_sorted = sort_qc_report(
                                    qc_report, module='metaQC', order=None
                                )

                                # dump updated qc_report to YAML file
                                f = open(report_path, 'w')
                                yaml.dump(
                                    qc_report_sorted, f, sort_keys=False,
                                    allow_unicode=False
                                )

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
                            
                            # give save_selector access to reclass dataframe
                            save_selector.reclass.bind(reclass)
                            
                            ###################################################

                            save_selector.native.setSizePolicy(
                                QtWidgets.QSizePolicy.Maximum,
                                QtWidgets.QSizePolicy.Maximum,
                            )

                            ###################################################

                            cluster_layout.addWidget(sample_selector.native)

                            cluster_layout.addWidget(reclass_selector.native)

                            cluster_layout.addWidget(save_selector.native)

                            ###################################################

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
                        def sweep_MCS(lowerMCS: int = 200.0,
                                      upperMCS: int = 200.0,):

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

                                chunk['cluster_2d'] = clustering.labels_

                                logger.info(
                                    f'min_cluster_size = {i} '
                                    f'{np.unique(clustering.labels_)}'
                                )
                        #######################################################

                        sweep_MCS.native.setSizePolicy(
                            QtWidgets.QSizePolicy.Maximum,
                            QtWidgets.QSizePolicy.Maximum,
                        )

                        viewer.window.add_dock_widget(
                            cluster_and_plot, name='Plot Single MCS',
                            area='right')

                        viewer.window.add_dock_widget(
                            sweep_MCS, name='Sweep MCS Range',
                            area='right')

                        viewer.window.add_dock_widget(
                            cluster_widget, name='Clustering Result',
                            area='right'
                        )

                        viewer.scale_bar.visible = True
                        viewer.scale_bar.unit = 'um'

                        napari.run()
                    
                    else:
                        print()
                        logger.info(
                            'MCS and/or reclassThresholds not found in '
                            'QC report. Reinitializing reclass.parquet, '
                            'setting chunk index to 0, and opening Napari '
                            'window so these values can be assigned.'
                        )
                        
                        try:
                            os.remove(
                                os.path.join(reclass_dir, 'reclass.parquet')
                            )
                        except FileNotFoundError:
                            pass
                        
                        reclass = pd.DataFrame(
                            columns=['Reclass', 'CellID', 'Sample', 
                                     'QC_status', 'Filter', 'ChunkIDX']
                        )
                        
                        # initialize chunk index counter at 0 and store
                        chunk_index = 0
                        progress['chunk_index'] = chunk_index
                        f = open(progress_path, 'w')
                        yaml.dump(
                            progress, f, sort_keys=False, allow_unicode=False
                        )

                        # Napari window will open now that chunk_index = 0

                ###############################################################
                # once optimal MCS has been saved
                try:
                    if check_report_structure(report_path):
                        qc_report = yaml.safe_load(open(report_path))
                        final_mcs_entry = qc_report['metaQC']['MCS']
                        final_reclass_entry = qc_report[
                            'metaQC']['reclassThresholds']
                    else:
                        print()
                        logger.info(
                            'Aborting; metaQC keys in QC report not '
                            'formatted correctly. Reinitializing '
                            'reclass.parquet, setting chunk index to 0. '
                            'Please re-run metaQC module to make MCS and '
                            'reclass threshold selections.'
                        )

                        parquet = glob.glob(
                            os.path.join(reclass_dir, 'reclass_*.parquet')
                        )
                        for i in parquet:
                            try:
                                os.remove(i)
                            except FileNotFoundError:
                                pass
                        
                        # initialize chunk index counter at 0 and store
                        chunk_index = 0
                        progress['chunk_index'] = chunk_index
                        f = open(progress_path, 'w')
                        yaml.dump(
                            progress, f, sort_keys=False, allow_unicode=False
                        )
                        sys.exit()
                
                except FileNotFoundError:
                    print()
                    logger.info(
                        'Aborting; QC report not found. Ensure '
                        'cylinter_report.yml is stored at '
                        'top-level of CyLinter output file or '
                        're-start pipeline to start filtering data.'
                    )
                    sys.exit()

                expected_filename = (
                    f'reclass_{self.percentDataPerChunk}'
                    f'_{final_mcs_entry}'
                    f'_{final_reclass_entry[0]}'
                    f'_{final_reclass_entry[1]}.parquet'
                )
                
                if filename is None: 
                    # MCS and reclassThresholds just defined in Napari
                    filename = expected_filename

                if filename == expected_filename:
                    if chunk_index not in reclass['ChunkIDX'].unique():
                    
                        #######################################################
                        # cluster chunk using selected MCS

                        if 'cluster_2d' not in chunk.columns:
                            # (not applicable to first chunk,
                            # gets clustered during plotting above)
                            
                            print()
                            logger.info(
                                f'Applying saved minimum cluster size: '
                                f'{final_mcs_entry}'
                            )

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
                                match_reference_implementation=False).fit(
                                chunk[['emb1', 'emb2']]
                            )
                            chunk['cluster_2d'] = clustering.labels_
                        
                        # exit program if all cells are considered ambiguous 
                        # by the clustering algorithm (likely too few 
                        # cells per chunk)
                        if chunk['cluster_2d'].eq(-1).all():
                            logger.warning(
                                'Aborting; all cells in chunk '
                                f'{chunk_index + 1} '
                                'were deemed ambiguous by clustering '
                                'algorithm (i.e. assigned to cluster -1). '
                                'Try using larger batch size.')
                            sys.exit()
                        
                        #######################################################
                        # add clean and noisy data to reclass storage dataframe 
                        # according to reclassification thresholds
            
                        print()
                        logger.info(
                            'Applying saved clean reclassification cutoff ' 
                            f'{final_reclass_entry[0]} '
                            f'to chunk {chunk_index + 1} of {len(chunks)}'
                        )
                        logger.info(
                            'Applying saved clean reclassification cutoff '
                            f'{final_reclass_entry[1]} '
                            f'to chunk {chunk_index + 1} of {len(chunks)}'
                        )
                        
                        chunk, clean, noisy = reclassify_chunk(
                            chunk=chunk, clean_cutoff=final_reclass_entry[0],
                            noisy_cutoff=final_reclass_entry[1]
                        )
                        clean['ChunkIDX'] = chunk_index
                        noisy['ChunkIDX'] = chunk_index

                        #######################################################
                        # generate clustermap for chunk data

                        clustermap_input = chunk[
                            chunk['cluster_2d'] != -1].copy()

                        cluster_heatmap_input = (
                            clustermap_input[
                                abx_channels +
                                ['cluster_2d']]
                            .groupby('cluster_2d')
                            .mean()
                        )
                        
                        # compute per channel z-scores across clusters
                        cluster_heatmap_input = (
                            (cluster_heatmap_input - 
                             cluster_heatmap_input.mean()) / 
                            cluster_heatmap_input.std()
                        )
                        # assign NaNs (channels with no 
                        # variation in signal) to 0
                        cluster_heatmap_input[cluster_heatmap_input.isna()] = 0

                        # zero-center colorbar
                        norm = TwoSlopeNorm(
                            vcenter=0, vmin=cluster_heatmap_input.min().min(),
                            vmax=cluster_heatmap_input.max().max()
                        )

                        sns.set(font_scale=0.8)
                        g = sns.clustermap(
                            cluster_heatmap_input, cmap='coolwarm',
                            standard_scale=None,
                            square=False, yticklabels=1, linewidth=0.1,
                            cbar=True, norm=norm
                        )
                        
                        g.fig.suptitle(
                            'channel_z-scores.pdf}', y=0.995, fontsize=10
                        )
                        g.fig.set_size_inches(8.0, 8.0)

                        plt.savefig(
                            os.path.join(chunk_dir, 'channel_z-scores.pdf'),
                            bbox_inches='tight'
                        )
                        plt.close('all')

                        #######################################################
                        # update reclass storage and save
                        
                        reclass = pd.concat(
                            [reclass,
                             clean[['Reclass', 'CellID', 'Sample', 'QC_status',
                                    'Filter', 'ChunkIDX']]
                             ], axis=0
                        )
                        reclass = pd.concat(
                            [reclass,
                             noisy[['Reclass', 'CellID', 'Sample', 'QC_status',
                                    'Filter', 'ChunkIDX']]
                             ], axis=0
                        )

                        reclass.to_parquet(
                            os.path.join(reclass_dir, f'{filename}')
                        )

                        print()
                        
                        logger.info(
                            f"Clean data reclassified as noisy tally: "
                            f"{len(reclass[(reclass['Reclass'] == 'noisy') & (reclass['QC_status'] == 'clean')])}"
                        )
                        logger.info(
                            f"Noisy data reclassified as clean tally: "
                            f"{len(reclass[(reclass['Reclass'] == 'clean') & (reclass['QC_status'] == 'noisy')])}"
                        )
                    
                    else:
                        logger.info(
                            'Reclassified data for '
                            f'chunk {chunk_index + 1} already ' 
                            'added to reclass parquet. Moving to next chunk.'
                        )
                else:
                    print()
                    logger.info(
                        f'Aborting; reclass parquet was built using '
                        f'MCS = {filename.split("_")[2]}, '
                        f'clean reclass threshold = {filename.split("_")[3]}, and '
                        f'noisy reclass threshold = {filename.split("_")[4].split(".parquet")[0]}. '
                        f'However, current QC report values are MCS = {final_mcs_entry}, clean '
                        f'reclass threshold = {final_reclass_entry[0]}, and noisy reclass '
                        f'threshold = {final_reclass_entry[1]}. Update values in QC report to '
                        'match or remove QC metadata values. Then re-run metaQC module.'
                    )
                    sys.exit()

                ###############################################################
                # increment chunk_index
                
                chunk_index = chunk_index + 1
                progress['chunk_index'] = chunk_index
                f = open(progress_path, 'w')
                yaml.dump(progress, f, sort_keys=False, allow_unicode=False)

                ###############################################################
            
        print()

        #######################################################################
        # ensure all chunks have been reclassified before 
        # reclassifying final dataframe
        if reclass['ChunkIDX'].max() == len(chunks) - 1:

            # create unique cell labels for raw data
            for module_idx, module in enumerate(modules):
                if module_dict[module_idx][0] == 'aggregateData':
                    pre_qc = module_dict[module_idx][1].copy()
                    pre_qc['handle'] = (
                        pre_qc['CellID'].map(str) + '_' + pre_qc['Sample']
                    )

            # create unique cell labels for cleaned
            # data (before reclassifiction)
            post_qc = module_dict[
                [i for i in module_dict.keys()][-1]][1].copy()
            post_qc['handle'] = (
                post_qc['CellID'].map(str) + '_' + post_qc['Sample']
            )

            # get raw values of cells in post_qc data
            cleaned_raw = pre_qc[pre_qc['handle'].isin(post_qc['handle'])]

            # isolate clean cells reclassified as noisy
            drop = reclass[
                (reclass['Reclass'] == 'noisy') & 
                (reclass['QC_status'] == 'clean')
            ].copy()

            if not drop.empty:    
                drop['handle'] = drop['CellID'].map(str) + '_' + drop['Sample']
            else:
                drop['handle'] = pd.Series()
            
            dropped = cleaned_raw[~cleaned_raw['handle'].isin(drop['handle'])]

            # isolate noisy cells reclassified as clean
            replace = reclass[
                (reclass['Reclass'] == 'clean') & 
                (reclass['QC_status'] == 'noisy')
            ].copy()
            
            if not replace.empty:
                replace['handle'] = (
                    replace['CellID'].map(str) + '_' + 
                    replace['Sample']
                )
            else:
                replace['handle'] = pd.Series()
            
            replaced = pre_qc[pre_qc['handle'].isin(replace['handle'])]

            data = pd.concat([dropped, replaced], axis=0).reset_index(drop=True)

            ###################################################################
            # transform data

            data.loc[:, abx_channels] += 0.001
            data.loc[:, abx_channels] = np.log10(data[abx_channels])

            ###################################################################
            # rescale channel signal intensities (0-1) per channel per sample

            for channel in abx_channels:
                for sample in natsorted(data['Sample'].unique()):
                    
                    sample_channel_data = data[
                        data['Sample'] == sample][channel]

                    # rescale channel signal intensities
                    scaler = (
                        MinMaxScaler(feature_range=(0, 1), copy=True)
                        .fit(sample_channel_data.values.reshape(-1, 1))
                    )
                    rescaled_data = scaler.transform(
                        sample_channel_data.values.reshape(-1, 1)
                    )

                    rescaled_data = pd.DataFrame(
                        data=rescaled_data,
                        index=sample_channel_data.index
                    ).rename(columns={0: channel})

                    # data to return from metaQC module
                    data.update(rescaled_data)

        else:
            logger.info(
                f'Aborting; only {reclass["ChunkIDX"].max() + 1} of '
                f'{len(chunks)} chunks '
                'have been reclassified. Updating chunk_index in QC report to ' 
                f'{reclass["ChunkIDX"].max() + 2}. Please re-run metaQC ' 
                'module to reclassify remaining chunks.'
            )
            progress['chunk_index'] = reclass['ChunkIDX'].max() + 1
            f = open(progress_path, 'w')
            yaml.dump(progress, f, sort_keys=False, allow_unicode=False)
            sys.exit()
    
    ###########################################################################
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
            (qc_dict[qc_keys[i]] - 
             qc_dict[qc_keys[i + 1]]) / qc_dict[qc_keys[0]]) * 100

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
        os.path.join(reclass_dir, 'censored_by_stage.pdf'), 
        bbox_inches='tight'
    )
    plt.close('all')

    # drop indexing handle from data before returning
    if self.metaQC:
        data.drop('handle', axis=1, inplace=True)

    data = reorganize_dfcolumns(data, markers, 2)

    print()
    print()
    return data
