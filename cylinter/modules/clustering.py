import os
import sys
import yaml
import logging

import numpy as np
import pandas as pd

from natsort import natsorted
from itertools import product
import math

from datetime import datetime

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from qtpy.QtCore import QTimer
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

from umap import UMAP
from sklearn.manifold import TSNE
import hdbscan
from joblib import Memory

from sklearn.metrics import silhouette_samples, silhouette_score

from magicgui import magicgui
import napari

from ..utils import (
    input_check, read_markers, matplotlib_warnings, categorical_cmap, SelectFromCollection, 
    cluster_expression, marker_channel_number, single_channel_pyramid, log_banner,
    log_multiline, get_filepath, reorganize_dfcolumns, sort_qc_report
)

logger = logging.getLogger(__name__)


def clustering(data, self, args):

    print()

    check, markers_filepath = input_check(self)
    
    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # drop antibody channel exclusions for clustering
    abx_channels = [i for i in abx_channels if i not in self.channelExclusionsClustering]

    # create N-dimensional clustering subdirectory if it hasn't already
    dim_dir = os.path.join(self.outDir, 'clustering', f'{self.dimensionEmbedding}d')
    if not os.path.exists(dim_dir):
        os.makedirs(dim_dir)

    # read QC report
    report_path = os.path.join(self.outDir, 'cylinter_report.yml')
    try:
        qc_report = yaml.safe_load(open(report_path))
        reload_report = False
        if qc_report is None:
            qc_report = {}
            reload_report = True
        if 'clustering' not in qc_report or qc_report['clustering'] is None:
            qc_report['clustering'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(qc_report, module='clustering', order=None)
            f = open(report_path, 'w')
            yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
            qc_report = yaml.safe_load(open(report_path))
    except:
        logger.info(
            'Aborting; QC report missing from CyLinter output directory. Re-start pipeline '
            'from aggregateData module to initialize QC report.'
        )
        sys.exit()
    
    def check_report_structure(report_path):
        
        expected = {
            'MCS_2d': 1,
            'MCS_3d': 1
        }
        
        try:
            with open(report_path, 'r') as file:
                qc_report = yaml.safe_load(file)
                clustering = qc_report.get('clustering')
            
            # check if top-level key is empty 
            if clustering is None:
                return False
            
            # check that subkeys match the expected structure
            if (f'MCS_{self.dimensionEmbedding}d' not in 
                    set(expected.keys())):
                return False
            
            # check that value types are formatted correctly
            for key, value in expected.items():
                if key == f'MCS_{self.dimensionEmbedding}d':
                    if not isinstance(clustering.get(key), type(value)):
                        return False
            
            return True  # all checks pass
        
        except FileNotFoundError:
            return False
    
    ##############################################################################################
    
    # recapitulate df index at the point of embedding
    data = data[~data['Sample'].isin(self.samplesToRemoveClustering)]

    # pick a random seed for reproducibility of
    # data sampling via self.fracForEmbedding
    random_state = 5

    if self.normalizeTissueCounts:
        logger.info('Performing weighted sampling of cells from tissues...')
        logger.info('Check that resulting cell counts below are similar across samples.')
        logger.info(
            'If not, try embedding a smaller fraction of data '
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

        log_banner(logger.info, 'Cell counts:')
        log_multiline(logger.info, data.groupby(['Sample']).size().to_string(index=True))
        print()

    else:

        data = data.sample(frac=self.fracForEmbedding, random_state=random_state)

    data.reset_index(drop=True, inplace=True)
    print(data[abx_channels])
    print()

    # if embedding has already been computed
    if os.path.exists(os.path.join(dim_dir, 'embedding.npy')):

        embedding = np.load(os.path.join(dim_dir, 'embedding.npy'))

        if embedding.shape[1] == 2:
            try:
                data['emb1'] = embedding[:, 0]
                data['emb2'] = embedding[:, 1]
                clustering_input = data[['emb1', 'emb2']]
            except ValueError:
                logger.info(
                    'Aborting; the number of rows in the dataframe differs from embedding.npy. '
                    'Ensure that the dataframe contains the same cells referred to in the ' 
                    'embedding or re-embed the current dataframe after removing the existing '
                    'embedding.'
                )
                sys.exit()

        elif embedding.shape[1] == 3:
            try:
                data['emb1'] = embedding[:, 0]
                data['emb2'] = embedding[:, 1]
                data['emb3'] = embedding[:, 2]
                clustering_input = data[['emb1', 'emb2', 'emb3']]
            except ValueError:
                logger.info(
                    'Aborting; the number of rows in the dataframe differs from the embedding. '
                    'Ensure that the dataframe contains the same cells referred to in the ' 
                    'embedding or re-embed the current dataframe after removing the existing '
                    'embedding.'
                )
                sys.exit()

    else:
        # exit program if dimensionEmbedding configuration is not 2 or 3
        if self.dimensionEmbedding not in [2, 3]:
            print()
            print('Embedding dimension must be set to 2 or 3.')
            sys.exit()

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

        logger.info('Embedding completed in ' + str(datetime.now() - startTime))

        np.save(os.path.join(dim_dir, 'embedding.npy'), embedding)

        if embedding.shape[1] == 2:
            data['emb1'] = embedding[:, 0]
            data['emb2'] = embedding[:, 1]
            clustering_input = data[['emb1', 'emb2']]

        elif embedding.shape[1] == 3:
            data['emb1'] = embedding[:, 0]
            data['emb2'] = embedding[:, 1]
            data['emb3'] = embedding[:, 2]
            clustering_input = data[['emb1', 'emb2', 'emb3']]
        print()

    # define the point size for cells in the embedding
    if len(data) >= 20000:
        point_size = 5000 / len(data)
    else:
        point_size = 4
    
    ##############################################################################################
    # show abx intensity for each marker on UMAP embedding
    
    if not os.path.exists(os.path.join(dim_dir, 'emb_channels.png')):

        logger.info('Coloring embedding by channel')

        ncols = 10
        nrows = math.ceil(len(abx_channels) / ncols)

        fig = plt.figure(figsize=(ncols, nrows))

        # grid specifications
        gs = plt.GridSpec(nrows=nrows, ncols=ncols, figure=fig)

        subsample_size = 50000
        if len(data) < subsample_size:
            plot_sample = data.sample(n=len(data))
        else:
            plot_sample = data.sample(n=subsample_size)
        
        for e, ax in enumerate(product(range(nrows), range(ncols))):

            if e < len(abx_channels):

                ch = abx_channels[e]
                
                if embedding.shape[1] == 3:
                    ax = fig.add_subplot(gs[ax[0], ax[1]], projection='3d')
                    ax.view_init(azim=-130, elev=10)  # all clusters view

                    plot = ax.scatter(
                        plot_sample['emb1'], plot_sample['emb2'], plot_sample['emb3'],
                        c=plot_sample[ch], s=point_size, lw=0.0, alpha=1.0,
                        cmap='magma'
                    )
                    ax.set_aspect('auto')
                    ax.xaxis.set_ticklabels([])
                    ax.yaxis.set_ticklabels([])
                    ax.zaxis.set_ticklabels([])
                
                else:
                    ax = fig.add_subplot(gs[ax[0], ax[1]])
                    plot = ax.scatter(
                        plot_sample['emb1'], plot_sample['emb2'], c=plot_sample[ch],
                        s=point_size, lw=0.0, alpha=1.0, cmap='magma'
                    )
                    ax.set_aspect('equal')
                    ax.axis('off')

                ax.set_title(ch, fontdict={'fontsize': 6}, pad=-0.1)

            if e == len(abx_channels) - 1:  # last iteration in loop
                divider = make_axes_locatable(ax)

                if embedding.shape[1] == 2:
                    cax = divider.append_axes('right', size='5%', pad=0.0)
                    cbar = fig.colorbar(plot, cax=cax)
                    cax.tick_params(labelsize=4)
                    cbar.outline.set_edgecolor('k')
                    cbar.outline.set_linewidth(0.2)
                    cbar.ax.tick_params(which='both', labelsize=3, width=0.2, length=2, pad=0)

                # TO IMPLEMENT CBAR FOR 3D PLOTS
                # cax = divider.append_axes('bottom', size='5%', pad=0.0)
                # cbar = fig.colorbar(plot, orientation='horizontal', cax=cax)
        
        # plt.subplots_adjust(
        #     left=0.01, right=0.99, bottom=0.01, top=0.99, hspace=-0.55, wspace=0.1
        # )
        plt.tight_layout()
        
        plt.savefig(os.path.join(dim_dir, 'emb_channels.png'), dpi=800)
        plt.close('all')

        print()
    
    ##############################################################################################
    # interact with plots to identify optimal minimum cluster size
    
    while not check_report_structure(report_path):

        # initialize Napari viewer without images
        viewer = napari.Viewer(title='CyLinter')

        # generate Qt widget
        cluster_widget = QtWidgets.QWidget()

        # generate vertical widget layout
        cluster_layout = QtWidgets.QVBoxLayout(cluster_widget)

        cluster_widget.setSizePolicy(
            QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Maximum,
        )

        ###################################################################
        @magicgui(
            layout='horizontal',
            call_button='Cluster and Plot',
            MCS={'label': 'Min Cluster Size (MCS)', 'step': 1},
        )
        def cluster_and_plot(MCS: int = 200):

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

            data[f'cluster_{self.dimensionEmbedding}d'] = clustering.labels_

            logger.info(f'min_cluster_size = {MCS} {np.unique(clustering.labels_)}')

            ###################################################
            # plot silhouette scores

            silho_input = clustering_input.copy()
            silho_input[f'cluster_{self.dimensionEmbedding}d'] = clustering.labels_
            silho_input = silho_input[
                silho_input[f'cluster_{self.dimensionEmbedding}d'] != -1
            ]

            all_minus_one = (silho_input[f'cluster_{self.dimensionEmbedding}d'] == -1).all()
            
            if not all_minus_one:
                
                # subsample clustered data for silhouette analysis
                silho_subset = silho_input.head(50000)

                cmap = categorical_cmap(
                    numUniqueSamples=len(
                        silho_input[f'cluster_{self.dimensionEmbedding}d'].unique()),
                    numCatagories=10, cmap='tab10', continuous=False
                )

                cluster_centers = pd.DataFrame(
                    index=sorted(silho_subset[f'cluster_{self.dimensionEmbedding}d'].unique())
                )
                
                embed_cols = [i for i in silho_subset.columns if 'emb' in i]
                for clus in sorted(silho_subset[f'cluster_{self.dimensionEmbedding}d'].unique()):
                    group = silho_subset[
                        silho_subset[f'cluster_{self.dimensionEmbedding}d'] == clus]
                    for emb_dim in embed_cols:
                        emb_mean = group[emb_dim].mean()
                        cluster_centers.loc[clus, emb_dim] = emb_mean

                n_clusters = len(cluster_centers.index.unique())

                silhouette_spacer = 1000

                if embedding.shape[1] == 2:
                    fig_silho, (ax1_silho, ax2_silho) = plt.subplots(1, 2)
                    fig_silho.set_size_inches(18, 7)
                elif embedding.shape[1] == 3:
                    fig_silho = plt.figure(figsize=(18, 7))
                    gs = plt.GridSpec(1, 2, figure=fig_silho)
                    ax1_silho = fig_silho.add_subplot(gs[0, 0], projection=None)
                    ax2_silho = fig_silho.add_subplot(gs[0, 1], projection='3d')

                ax1_silho.set_xlim([-1, 1])

                ax1_silho.set_ylim([0, len(silho_subset) + (n_clusters + 1) * silhouette_spacer])

                sample_silhouette_values = silhouette_samples(
                    silho_subset[embed_cols], silho_subset[f'cluster_{self.dimensionEmbedding}d']
                )

                y_lower = silhouette_spacer
                for i in cluster_centers.index.unique():

                    ith_cluster_silhouette_values = (
                        sample_silhouette_values[silho_subset[
                            f'cluster_{self.dimensionEmbedding}d'] == i]
                    )
                    ith_cluster_silhouette_values.sort()

                    size_cluster_i = ith_cluster_silhouette_values.shape[0]

                    y_upper = y_lower + size_cluster_i

                    color = cmap.colors[i]

                    ax1_silho.fill_betweenx(
                        np.arange(y_lower, y_upper), 0, ith_cluster_silhouette_values,
                        facecolor=color, edgecolor=color, alpha=0.7
                    )

                    ax1_silho.text(
                        0.0, y_lower + 0.5 * size_cluster_i, str(i),
                        fontdict={'size': 20 / np.log(n_clusters)}
                    )

                    # compute the new y_lower for next plot
                    y_lower = y_upper + silhouette_spacer

                silhouette_avg = silhouette_score(
                    silho_subset[embed_cols], silho_subset[f'cluster_{self.dimensionEmbedding}d']
                )

                ax1_silho.set_title('Silhouette Plot')
                ax1_silho.set_xlabel('Silhouette Coefficients')
                ax1_silho.set_ylabel('Cluster label')

                ax1_silho.axvline(x=silhouette_avg, color='r', lw=0.75, linestyle='--')

                ax1_silho.set_yticks([])

                if embedding.shape[1] == 2:
                    ax2_silho.scatter(
                        silho_subset['emb1'], silho_subset['emb2'], marker='.',
                        s=30, edgecolor='k', lw=0, alpha=1.0,
                        c=[cmap.colors[i] for i in
                           silho_subset[f'cluster_{self.dimensionEmbedding}d']]
                    )

                    ax2_silho.scatter(
                        cluster_centers['emb1'], cluster_centers['emb2'],
                        marker='o', c='white', alpha=1, s=125, edgecolor='k'
                    )

                    for i in cluster_centers.iterrows():
                        ax2_silho.scatter(
                            i[1]['emb1'], i[1]['emb2'], marker='$%d$' % i[0],
                            alpha=1, s=40, edgecolor='k'
                        )

                elif embedding.shape[1] == 3:
                    ax2_silho.scatter(
                        silho_subset['emb1'], silho_subset['emb2'],
                        silho_subset['emb3'], marker='.', s=30, edgecolor='k',
                        lw=0, alpha=1.0,
                        c=[cmap.colors[i] for i in
                           silho_subset[f'cluster_{self.dimensionEmbedding}d']]
                    )

                    ax2_silho.scatter(
                        cluster_centers['emb1'], cluster_centers['emb2'],
                        cluster_centers['emb3'], marker='o', c='white',
                        alpha=1, s=125, edgecolor='k'
                    )

                    for i in cluster_centers.iterrows():
                        ax2_silho.scatter(
                            i[1]['emb1'], i[1]['emb2'], i[1]['emb3'],
                            marker='$%d$' % i[0], zorder=100, alpha=1, s=40,
                            edgecolor='k'
                        )

                ax2_silho.set_title('Clustering')
                ax2_silho.set_xlabel('Feature Space 1')
                ax2_silho.set_ylabel('Feature Space 2')

                total_clusters = len(silho_input[f'cluster_{self.dimensionEmbedding}d'].unique())

                fig_silho.suptitle(
                    ('MCS=%d, average silhouette score for %d/%d clusters '
                     'in silhouette data subset is %f'
                     % (MCS, n_clusters, total_clusters, silhouette_avg)),
                    fontsize=14, fontweight='bold'
                )

                # show silhouette plot after cluster widget
                # is added to Napari window (below)
            
            else:
                logger.info(
                    'All data points were determined by HDBSCAN to be ambiguous '
                    ', skipping silhouette plot.'
                )
            
            ###################################################

            if embedding.shape[1] == 2:

                # placeholder for lasso selection
                selector = None

                sns.set_style('whitegrid')

                fig = plt.figure(figsize=(3, 8))
                matplotlib_warnings(fig)

                gs = plt.GridSpec(2, 3, figure=fig)

                # define axes
                ax_cluster = fig.add_subplot(gs[0, 0])
                ax_gate = fig.add_subplot(gs[0, 1])
                ax_sample = fig.add_subplot(gs[0, 2])

                ax_cluster_lbs = fig.add_subplot(gs[1, 0])
                ax_gate_lbs = fig.add_subplot(gs[1, 1])
                ax_sample_lbs = fig.add_subplot(gs[1, 2])

                plt.subplots_adjust(
                    left=0.01, right=0.99, bottom=0.0,
                    top=0.95, wspace=0.0, hspace=0.0)

                # PLOT embedding
                for color_by in [f'cluster_{self.dimensionEmbedding}d', 'class', 'annotation']:

                    highlight = 'none'

                    if color_by == f'cluster_{self.dimensionEmbedding}d':

                        # build cmap
                        cmap = categorical_cmap(
                            numUniqueSamples=len(data[color_by].unique()),
                            numCatagories=10, cmap='tab10', continuous=False
                        )

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

                        c = [sample_dict[i] for i in data[color_by]]

                        ax_cluster.cla()

                        cluster_paths = ax_cluster.scatter(
                            data['emb1'], data['emb2'], c=c, alpha=1.0,
                            s=point_size, cmap=cmap, ec='k', linewidth=0.0)

                        ax_cluster.set_title('HDBSCAN (lasso)', fontsize=7)
                        ax_cluster.set_aspect('equal')
                        ax_cluster.axes.xaxis.set_visible(False)
                        ax_cluster.axes.yaxis.set_visible(False)
                        ax_cluster.grid(False)

                        legend_elements = []
                        for e, i in enumerate(natsorted(data[color_by].unique())):

                            hi_markers = cluster_expression(
                                df=data, markers=abx_channels, cluster=i, num_proteins=3,
                                clus_dim=self.dimensionEmbedding
                            )

                            legend_elements.append(
                                Line2D([0], [0], marker='o',
                                       color='none', label=f'Cluster {i}: {hi_markers}',
                                       markerfacecolor=cmap.colors[e], markeredgecolor='none',
                                       lw=0.001, markersize=3)
                            )

                        ax_cluster_lbs.legend(
                            handles=legend_elements, prop={'size': 3},
                            loc='upper left', frameon=False)

                        ax_cluster_lbs.axis('off')

                    elif color_by == 'class':
                        title = 'Gating'
                        
                        if 'class' not in data.columns:
                            c = 'gainsboro'
                            legend_elements = []
                        
                        else:
                            
                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=data['class'].nunique(),
                                numCatagories=10, cmap='tab10', continuous=False
                            )

                            if 'unclassified' in data['class'].unique():
                                cmap = ListedColormap(
                                    np.insert(arr=cmap.colors, obj=0,
                                              values=[0.0, 0.0, 0.0], axis=0)
                                )
                                cmap = ListedColormap(cmap.colors[:-1])

                                order = (
                                    ['unclassified'] +
                                    list(
                                        data['class']
                                        .value_counts()
                                        .drop('unclassified')
                                        .index)
                                )
                                class_dict = dict(zip(order, list(range(len(order)))))
                            else:
                                class_dict = dict(
                                    zip(data['class'].value_counts().index,
                                        list(range(len(data['class'].nunique()))))
                                )

                            c = [class_dict[i] for i in data['class']]

                            legend_elements = []
                            for vector_name, color_idx in class_dict.items():
                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none', label=vector_name,
                                           markerfacecolor=cmap.colors[color_idx],
                                           markeredgecolor='none',
                                           lw=0.001, markersize=3)
                                )
                        
                        ax_gate.cla()

                        ax_gate.scatter(
                            data['emb1'], data['emb2'], c=c, cmap=cmap,
                            alpha=1.0, s=point_size, ec='k', linewidth=0.0)

                        ax_gate.set_title(title, fontsize=7)
                        ax_gate.set_aspect('equal')
                        ax_gate.axes.xaxis.set_visible(False)
                        ax_gate.axes.yaxis.set_visible(False)
                        ax_gate.grid(False)

                        ax_gate_lbs.legend(
                            handles=legend_elements, prop={'size': 3},
                            loc='upper left', frameon=False
                        )

                        ax_gate_lbs.axis('off')

                    elif color_by == 'annotation':

                        # build cmap
                        cmap = categorical_cmap(
                            numUniqueSamples=len(
                                data[self.colormapAnnotationClustering].unique()),
                            numCatagories=10, cmap='tab10', continuous=False
                        )

                        sample_dict = dict(
                            zip(natsorted(
                                data[self.colormapAnnotationClustering].unique()),
                                list(range(len(
                                    data[self.colormapAnnotationClustering]
                                     .unique())))))

                        c = [sample_dict[i] for i in data[self.colormapAnnotationClustering]]

                        ax_sample.cla()

                        ax_sample.scatter(
                            data['emb1'], data['emb2'], c=c, cmap=cmap,
                            alpha=1.0, s=point_size, ec='k', linewidth=0.0
                        )

                        ax_sample.set_title(self.colormapAnnotationClustering, fontsize=7)
                        ax_sample.set_aspect('equal')
                        ax_sample.axes.xaxis.set_visible(False)
                        ax_sample.axes.yaxis.set_visible(False)
                        ax_sample.grid(False)

                        legend_elements = []
                        for e, i in enumerate(natsorted(
                            data[self.colormapAnnotationClustering].unique()
                        )):

                            if i == highlight:
                                markeredgecolor = 'k'
                            else:
                                markeredgecolor = 'none'

                            legend_elements.append(
                                Line2D([0], [0], marker='o', color='none', label=i,
                                       markerfacecolor=cmap.colors[e],
                                       markeredgecolor=markeredgecolor,
                                       lw=0.001, markersize=3)
                            )

                        ax_sample_lbs.legend(
                            handles=legend_elements, prop={'size': 3},
                            loc='upper left', frameon=False)

                        ax_sample_lbs.axis('off')

            elif embedding.shape[1] == 3:

                # initialize figure and add to FigureCanvas
                # before rendering plot in 3D
                sns.set_style('whitegrid')
                fig = plt.figure(figsize=(4, 8))
                matplotlib_warnings(fig)

            count = cluster_layout.count()
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
                ax_gate = fig.add_subplot(gs[0, 1], projection='3d')
                ax_sample = fig.add_subplot(gs[0, 2], projection='3d')

                ax_cluster_lbs = fig.add_subplot(gs[1, 0])
                ax_gate_lbs = fig.add_subplot(gs[1, 1])
                ax_sample_lbs = fig.add_subplot(gs[1, 2])

                plt.subplots_adjust(
                    left=0.0, right=0.99, bottom=0.0, top=0.9, wspace=0.0, hspace=0.0
                )

                for color_by in [
                    f'cluster_{self.dimensionEmbedding}d', 'class', 'annotation'
                ]:

                    highlight = 'none'

                    if color_by == f'cluster_{self.dimensionEmbedding}d':

                        # build cmap
                        cmap = categorical_cmap(
                            numUniqueSamples=len(data[color_by].unique()),
                            numCatagories=10, cmap='tab10', continuous=False
                        )

                        # make black the first color to specify
                        # unclustered cells (cluster -1)
                        if -1 in data[color_by].unique():
                            cmap = ListedColormap(
                                np.insert(arr=cmap.colors, obj=0,
                                          values=[0.0, 0.0, 0.0], axis=0)
                            )

                            # trim qualitative cmap to number
                            # of unique samples
                            trim = len(cmap.colors) - len(data[color_by].unique())
                            cmap = ListedColormap(cmap.colors[:-trim])

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

                        ax_cluster.set_title('HDBSCAN', fontsize=7)
                        ax_cluster.set_xticklabels([])
                        ax_cluster.set_yticklabels([])
                        ax_cluster.set_zticklabels([])
                        ax_cluster.xaxis._axinfo['grid'].update({'linewidth': 0.5})
                        ax_cluster.yaxis._axinfo['grid'].update({'linewidth': 0.5})
                        ax_cluster.zaxis._axinfo['grid'].update({'linewidth': 0.5})

                        legend_elements = []
                        for e, i in enumerate(natsorted(data[color_by].unique())):

                            hi_markers = cluster_expression(
                                df=data, markers=abx_channels, cluster=i, num_proteins=3,
                                clus_dim=self.dimensionEmbedding
                            )

                            legend_elements.append(
                                Line2D([0], [0], marker='o',
                                       color='none', label=f'Cluster {i}: {hi_markers}',
                                       markerfacecolor=cmap.colors[e], markeredgecolor='none',
                                       lw=0.001, markersize=3)
                            )

                        ax_cluster_lbs.legend(
                            handles=legend_elements, prop={'size': 3}, 
                            loc='upper left', frameon=False
                        )
                        
                        ax_cluster_lbs.axis('off')

                    elif color_by == 'class':

                        title = 'Gating'
                        
                        if 'class' not in data.columns:
                            c = 'gainsboro'
                            legend_elements = []
                        
                        else:
                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=data['class'].nunique(),
                                numCatagories=10, cmap='tab10', continuous=False
                            )

                            if 'unclassified' in data['class'].unique():
                                cmap = ListedColormap(
                                    np.insert(arr=cmap.colors, obj=0,
                                              values=[0.0, 0.0, 0.0], axis=0)
                                )
                                cmap = ListedColormap(cmap.colors[:-1])

                                order = (
                                    ['unclassified'] +
                                    list(
                                        data['class']
                                        .value_counts()
                                        .drop('unclassified')
                                        .index)
                                )
                                class_dict = dict(zip(order, list(range(len(order)))))
                            else:
                                class_dict = dict(
                                    zip(data['class'].value_counts().index,
                                        list(range(len(data['class'].nunique()))))
                                )

                            c = [class_dict[i] for i in data['class']]

                            legend_elements = []
                            for vector_name, color_idx in class_dict.items():
                                legend_elements.append(
                                    Line2D([0], [0], marker='o',
                                           color='none', label=vector_name,
                                           markerfacecolor=cmap.colors[color_idx],
                                           markeredgecolor='none', lw=0.001, markersize=3)
                                )
                        
                        ax_gate.cla()
                        
                        ax_gate.scatter(
                            data['emb1'], data['emb2'], data['emb3'],
                            c=c, cmap=cmap, alpha=1.0, s=point_size,
                            ec='k', linewidth=0.0)

                        ax_gate.set_title(title, fontsize=7)
                        ax_gate.set_xticklabels([])
                        ax_gate.set_yticklabels([])
                        ax_gate.set_zticklabels([])

                        ax_gate.xaxis._axinfo['grid'].update(
                            {'linewidth': 0.5})
                        ax_gate.yaxis._axinfo['grid'].update(
                            {'linewidth': 0.5})
                        ax_gate.zaxis._axinfo['grid'].update(
                            {'linewidth': 0.5})

                        ax_gate_lbs.legend(
                            handles=legend_elements, prop={'size': 3},
                            loc='upper left', frameon=False
                        )

                        ax_gate_lbs.axis('off')

                    elif color_by == 'annotation':

                        # build cmap
                        cmap = categorical_cmap(
                            numUniqueSamples=len(
                                data[self.colormapAnnotationClustering].unique()),
                            numCatagories=10, cmap='tab10', continuous=False
                        )

                        sample_dict = dict(
                            zip(natsorted(
                                data[
                                    self.colormapAnnotationClustering].unique()),
                                list(range(
                                    len(data[self.colormapAnnotationClustering]
                                     .unique())))))

                        c = [sample_dict[i] for i in data[self.colormapAnnotationClustering]]

                        ax_sample.scatter(
                            data['emb1'], data['emb2'], data['emb3'], c=c, cmap=cmap,
                            alpha=1.0, s=point_size, ec='k', linewidth=0.0
                        )

                        ax_sample.set_title(self.colormapAnnotationClustering, fontsize=7)
                        ax_sample.set_xticklabels([])
                        ax_sample.set_yticklabels([])
                        ax_sample.set_zticklabels([])
                        ax_sample.xaxis._axinfo['grid'].update({'linewidth': 0.5})
                        ax_sample.yaxis._axinfo['grid'].update({'linewidth': 0.5})
                        ax_sample.zaxis._axinfo['grid'].update({'linewidth': 0.5})

                        legend_elements = []
                        for e, i in enumerate(natsorted(
                            data[self.colormapAnnotationClustering].unique())
                        ):

                            if i == highlight:
                                markeredgecolor = 'k'
                            else:
                                markeredgecolor = 'none'

                            legend_elements.append(
                                Line2D([0], [0], marker='o', color='none', label=i,
                                       markerfacecolor=cmap.colors[e], 
                                       markeredgecolor=markeredgecolor, lw=0.001, markersize=3)
                            )

                        ax_sample_lbs.legend(
                            handles=legend_elements, prop={'size': 3},
                            loc='upper left', frameon=False
                        )

                        ax_sample_lbs.axis('off')

            cluster_layout.addWidget(NavigationToolbar(cluster_canvas, cluster_widget))
            cluster_layout.addWidget(cluster_canvas)

            if not all_minus_one:
                fig_silho.show()

            if embedding.shape[1] == 2:

                # must call draw() before creating selector,
                # or alpha setting doesn't work.
                fig.canvas.draw()

                if selector:
                    selector.disconnect()
                selector = SelectFromCollection(
                    ax_cluster, cluster_paths)

            ###############################################################
            @magicgui(
                layout='horizontal',
                call_button='View Lassoed Points',
                sample_name={'label': 'Sample Name'},
            )
            def sample_selector(sample_name: str):

                return sample_name

            ###############################################################
            @sample_selector.called.connect
            def sample_selector_callback(value: str):

                print()
                napari.utils.notifications.show_info(f'Sample selection: {value}')

                # if cells lassoed
                if selector.ind is None or len(selector.ind) == 0:

                    print()
                    napari.utils.notifications.show_warning(
                        'Cells must be lassoed in HDBSCAN plot before sample inspection.'
                    )
                    pass

                else:

                    # superimpose centroids of lassoed noisy cells
                    # colored by stage removed over channel images

                    # if valid sample name entered
                    if value in data['Sample'].unique():

                        # show highest expression channels
                        lasso = data.loc[data.index.isin(selector.ind)].copy()

                        # assign lassoed data a dummy
                        # cluster variable and get highest
                        # expressed markers across channels
                        lasso[f'cluster_{self.dimensionEmbedding}d'] = 1000

                        hi_markers = cluster_expression(
                            df=lasso, markers=abx_channels, cluster=1000, num_proteins=3,
                            clus_dim=self.dimensionEmbedding
                        )

                        # remove previous samples' image layers
                        for i in reversed(range(0, len(viewer.layers))):
                            viewer.layers.pop(i)

                        if self.showAbChannels:

                            for ch in reversed(abx_channels):
                                channel_number = marker_channel_number(self, markers, ch)

                                # read antibody image
                                file_path = get_filepath(self, check, value, 'TIF')
                                img, min, max = single_channel_pyramid(
                                    file_path, channel=channel_number
                                )
                                viewer.add_image(
                                    img, rgb=False, blending='additive', colormap='green',
                                    visible=False, name=ch, contrast_limits=(min, max)
                                )

                        centroids = data[['Y_centroid', 'X_centroid']][
                            (data.index.isin(selector.ind)) & (data['Sample'] == value)
                        ]

                        viewer.add_points(
                            centroids, name='lassoed cells', visible=True,
                            face_color='yellow', edge_width=0.0, size=4.0
                        )

                        # read segmentation outlines, add to Napari
                        file_path = get_filepath(self, check, value, 'SEG')
                        seg, min, max = single_channel_pyramid(file_path, channel=0)
                        viewer.add_image(
                            seg, rgb=False, blending='additive', colormap='gray', visible=False,
                            name='segmentation', contrast_limits=(min, max)
                        )

                        # read first DNA, add to Napari
                        file_path = get_filepath(self, check, value, 'TIF')
                        channel_number = marker_channel_number(
                            self, markers, self.counterstainChannel
                        )
                        dna_first, min, max = single_channel_pyramid(
                            file_path, channel=channel_number
                        )
                        viewer.add_image(
                            dna_first, rgb=False, blending='additive',
                            opacity=0.5, colormap='gray', visible=True,
                            name=self.counterstainChannel, contrast_limits=(min, max)
                        )

                        # apply previously defined contrast limits if they exist
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
                                'Existing contrast settings applied.'
                            )
                        except:
                            pass
                        
                        print()
                        logger.info(
                            f'Top three expressed markers {hi_markers} '
                            '(raw expression values, not z-scores)'
                        )
                        napari.utils.notifications.show_info(
                            f'Viewing lassoed cells in sample {value}'
                        )

                    else:
                        print()
                        napari.utils.notifications.show_warning(
                            'Sample name not in filtered data.'
                        )
                        pass

            ###############################################################

            sample_selector.native.setSizePolicy(
                QtWidgets.QSizePolicy.Maximum,
                QtWidgets.QSizePolicy.Maximum,
            )

            ###############################################################
            
            @magicgui(layout='horizontal', call_button='Save')
            def save_selector():

                current_MCS = cluster_and_plot[0].value

                print()
                logger.info(f'Saving current min cluster size: {current_MCS}')
                qc_report['clustering'][f'MCS_{self.dimensionEmbedding}d'] = current_MCS

                # sort and dump updated qc_report to YAML file
                qc_report_sorted = sort_qc_report(
                    qc_report, module='clustering', order=None
                )
                f = open(report_path, 'w')
                yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
                
                # ensure silhouette plot is closed 
                plt.close('all')  
                
                if embedding.shape[1] == 2:
                    print()
                    logger.info('Saving 2D cluster plot')

                    if selector:
                        selector.disconnect()

                    fig.savefig(
                        os.path.join(dim_dir, f'{self.embeddingAlgorithm}_{current_MCS}.png'),
                        bbox_inches='tight', dpi=1000
                    )

                    fig_silho.savefig(
                        os.path.join(dim_dir, 
                                     f'{self.embeddingAlgorithm}_{current_MCS}_silho.png'),
                        bbox_inches='tight', dpi=1000
                    )
                    fig_silho.savefig(
                        os.path.join(dim_dir, 
                                     f'{self.embeddingAlgorithm}_{current_MCS}_silho.pdf'),
                        bbox_inches='tight'
                    )

                elif embedding.shape[1] == 3:
                    print()
                    logger.info('Saving 3D cluster plot')

                    def animate(i):
                        ax_cluster.view_init(elev=10., azim=i)
                        ax_gate.view_init(elev=10., azim=i)
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
            QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Maximum,
        )

        ###################################################################
        @magicgui(
            layout='horizontal',
            call_button='Sweep Range',
            lowerMCS={'label': 'Lower MCS', 'step': 1},
            upperMCS={'label': 'Upper MCS', 'step': 1},
        )
        def sweep_MCS(lowerMCS: int = 200, upperMCS: int = 200,):

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

                data[f'cluster_{self.dimensionEmbedding}d'] = clustering.labels_

                logger.info(f'min_cluster_size = {i} {np.unique(clustering.labels_)}')

                silho_input = clustering_input.copy()
                silho_input[f'cluster_{self.dimensionEmbedding}d'] = clustering.labels_
                silho_input = silho_input[
                    silho_input[f'cluster_{self.dimensionEmbedding}d'] != -1
                ]

                all_minus_one = (silho_input[f'cluster_{self.dimensionEmbedding}d'] == -1).all()
            
                if not all_minus_one:
                    
                    # subsample clustered data for silhouette analysis
                    silho_subset = silho_input.head(50000)

                    embed_cols = [i for i in silho_subset.columns if 'emb' in i]

                    silhouette_avg = silhouette_score(
                        silho_subset[embed_cols],
                        silho_subset[f'cluster_{self.dimensionEmbedding}d']
                    )

                    total_clusters = len(
                        silho_input[f'cluster_{self.dimensionEmbedding}d'].unique()
                    )

                    n_clusters = len(silho_subset[f'cluster_{self.dimensionEmbedding}d'].unique())

                    logger.info(
                        ('Average silhouette score for %d/%d clusters in '
                         'silhouette data subset: %f'
                         % (n_clusters, total_clusters, silhouette_avg))
                    )
                    print()

                else:
                    logger.info(
                        'HDBSCAN determined that all data points are ambiguous (-1), '
                        'skipping silhouette scoring.'
                    )

        ##########################################################################################

        sweep_MCS.native.setSizePolicy(
            QtWidgets.QSizePolicy.Maximum,
            QtWidgets.QSizePolicy.Maximum,
        )

        viewer.window.add_dock_widget(cluster_and_plot, name='Plot Single MCS', area='right')

        viewer.window.add_dock_widget(sweep_MCS, name='Sweep MCS Range', area='right')

        viewer.window.add_dock_widget(cluster_widget, name='Clustering Result', area='right')

        viewer.scale_bar.visible = True
        viewer.scale_bar.unit = 'um'

        napari.run()

        print()

    ##########################################################################################
    # apply final MCS and return data from clustering module
    try:
        if check_report_structure(report_path):
            qc_report = yaml.safe_load(open(report_path))
            final_mcs_entry = qc_report['clustering'][f'MCS_{self.dimensionEmbedding}d']
        else:
            print()
            logger.info(
                'Aborting; clustering keys in QC report not formatted correctly.'
                'Reinitializing reclass.parquet, setting chunk index to 0. Please '
                're-run clustering module to make MCS and reclass threshold selections.'
            )
    except FileNotFoundError:
        print()
        logger.info(
            'Aborting; QC report not found. Ensure cylinter_report.yml is stored at '
            'top-level of CyLinter output file or re-start pipeline '
            'to start filtering data.'
        )
        sys.exit()

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
        match_reference_implementation=False).fit(clustering_input)
    data[f'cluster_{self.dimensionEmbedding}d'] = clustering.labels_

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    # save dataframe in standard CSV format for analysis outside CyLinter
    data.to_csv(os.path.join(self.outDir, 'checkpoints', 'clustering.csv'), index=False)
    
    print()
    print()
    return data
