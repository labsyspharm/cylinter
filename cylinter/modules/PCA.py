import os
import logging

import math
import numpy as np
import pandas as pd

from natsort import natsorted

from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects

from sklearn.decomposition import PCA as PCA_MODULE

from ..utils import input_check, read_markers, categorical_cmap, reorganize_dfcolumns

logger = logging.getLogger(__name__)


def PCA(data, self, args):

    print()
    
    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
    )

    # drop antibody channel exclusions for PCA
    abx_channels = [
        i for i in abx_channels if i not in self.channelExclusionsPCA]

    # create PCA directory if it doesn't already exist
    pca_dir = os.path.join(self.outDir, 'PCA')
    if not os.path.exists(pca_dir):
        os.makedirs(pca_dir)

    ######################################################################

    sns.set_style("whitegrid", {'axes.grid': False})
    gs = plt.GridSpec(len(abx_channels), 1)

    for plot in ['ridgeplots', 'ridgeplots_persample']:
        fig = plt.figure(figsize=(2, 7))
        ax_objs = []
        for i, channel in enumerate(abx_channels):

            # creating new axes object
            ax_objs.append(fig.add_subplot(gs[i:i + 1, 0:]))

            if plot == 'ridgeplots':
                # plotting the distribution
                n, bins, patches = ax_objs[-1].hist(
                    data[channel], bins=100, density=True,
                    histtype='stepfilled', linewidth=2.0,
                    ec='k', alpha=1.0, color='k'
                )

                # setting uniform x and y lims
                ax_objs[-1].set_xlim(0, 1)
                ax_objs[-1].set_ylim(0, math.ceil(n.max()) + 1)

            elif plot == 'ridgeplots_persample':
                y_vals = []
                for sample in sorted(data['Sample'].unique()):

                    # plotting the distribution
                    n, bins, patches = ax_objs[-1].hist(
                        data[channel][data['Sample'] == sample], bins=100,
                        density=True, histtype='stepfilled', linewidth=0.0,
                        ec='k', alpha=1.0
                    )
                    y_vals.extend(n.tolist())

                # setting uniform x and y lims
                ax_objs[-1].set_xlim(0, 1)
                ax_objs[-1].set_ylim(0, math.ceil(max(y_vals)) + 1)

            # make background transparent
            rect = ax_objs[-1].patch
            rect.set_alpha(0)

            # remove borders, axis ticks, and labels
            ax_objs[-1].set_yticklabels([])

            if i == len(abx_channels) - 1:
                ax_objs[-1].set_xlabel(
                    'Intensity', fontsize=10, fontweight='normal', labelpad=5
                )
            else:
                ax_objs[-1].set_xticks([])
                ax_objs[-1].set_xticklabels([])

            ax_objs[-1].set_yticks([])

            spines = ['top', 'right', 'left']
            for s in spines:
                ax_objs[-1].spines[s].set_visible(False)

            ax_objs[-1].tick_params(axis='x', width=2)

            ax_objs[-1].text(
                -0.02, 0, channel, fontweight='normal', fontsize=8, ha='right'
            )

        gs.update(hspace=0.3)
        plt.subplots_adjust(left=0.3, bottom=0.1, right=0.9, top=0.95)
        fig.savefig(os.path.join(pca_dir, f'{plot}.png'), dpi=800)
        plt.close('all')

    #######################################################################

    # Horn's parallel analysis to determine the number of non-random PCs
    unshuffled = data[abx_channels]

    unshuffled = unshuffled[~unshuffled.index.isin(self.samplesToRemovePCA)]

    evr = pd.DataFrame()

    n_components = min(len(data['Sample'].unique()), len(abx_channels))

    if n_components > 1:
        logging.info(
            "Performing Horn's parallel analysis to determine number of non-random PCs..."
        )

        for iteration in range(1, n_components + 1):
            logging.info(f'iteration {iteration}/{n_components}')
            shuffled = unshuffled.copy()
            for e, col in enumerate(shuffled.columns):
                shuffled[col] = shuffled[col].sample(
                    frac=1.0, random_state=e + iteration).values

            # specify PCA parameters
            pca = PCA_MODULE(n_components=n_components, random_state=1)

            # mean-center data (axis=0)
            for c in shuffled.columns:
                shuffled[c] = shuffled[c] - shuffled[c].mean()

            # fit PCA parameters to data
            pca.fit_transform(shuffled)

            evr[iteration] = pca.explained_variance_ratio_

        # used 1-based indexing for PCs
        evr.index = [i + 1 for i in evr.index]

        # specify PCA parameters for unshuffled data
        pca = PCA_MODULE(n_components=n_components, random_state=1)

        # mean-center unshuffled data (axis=0)
        for c in unshuffled.columns:
            unshuffled[c] = unshuffled[c] - unshuffled[c].mean()

        # apply PCA parameters to unshuffled data
        projected = pca.fit_transform(unshuffled)

        sns.set_style('whitegrid')
        fig1, ax1 = plt.subplots()

        x = np.arange(1, n_components + 1, step=1)
        
        ax1.plot(evr.index, evr.mean(axis=1), color='tab:orange', marker='o')
        ax1.plot(x, pca.explained_variance_ratio_, color='tab:blue', marker='o')

        ax1.set_xticks(x)

        ax1.set_xlabel('PC', fontsize=10, labelpad=7.0)
        ax1.set_ylabel('Explained Variance Ratio', fontsize=10, labelpad=4.0)

        legend_handles = []
        legend_handles.append(
            Line2D([0], [0], marker=None, color='tab:blue',
                   label='unshuffled', markeredgewidth=0.7,
                   markersize=5.0, linewidth=5)
        )
        legend_handles.append(
            Line2D([0], [0], marker=None, color='tab:orange',
                   label='average shuffled', markeredgewidth=0.7,
                   markersize=5.0, linewidth=5)
        )
        ax1.legend(handles=legend_handles, prop={'size': 10.0}, bbox_to_anchor=[0.95, 1.0])
        fig1.savefig(os.path.join(pca_dir, 'variance.pdf'), bbox_inches='tight')
        plt.close(fig1)

        ###################################################################
        # plot cell distribution across PCs 1 and 2

        # generate dataframe for plot input
        scatter_input = pd.DataFrame(data=projected, index=unshuffled.index)
        
        scatter_input.rename(columns={0: 'PC1', 1: 'PC2'}, inplace=True)

        fig2, ax2 = plt.subplots()
        sns.scatterplot(
            data=scatter_input, x='PC1', y='PC2', color='k', linewidth=0.0,
            s=30000 / len(data), alpha=1.0, legend=False
        )

        ax2.set_xlabel(
            f'PC1 ({round((pca.explained_variance_ratio_[0] * 100), 2)}'
            '% of variance)', fontsize=10, labelpad=7.0
        )
        ax2.set_ylabel(
            f'PC2 ({round((pca.explained_variance_ratio_[1] * 100), 2)}'
            '% of variance)', fontsize=10, labelpad=4.0
        )
        ax2.tick_params(axis='both', which='major', labelsize=7.0)

        fig2.savefig(
            os.path.join(pca_dir, 'pcaScoresPlotCells.png'), dpi=600, bbox_inches='tight'
        )
        plt.close(fig2)

        ###################################################################
        # compute median channel intensites for samples and
        # plot their distribution across PCs 1 and 2

        if len(data['Sample'].unique()) > 1:

            # compute median antibody expression per sample
            # samples (rows) x features (columns)
            medians = data.groupby(['Sample']).median(numeric_only=True)[abx_channels]

            # drop sample exclusions for PCA
            medians = medians[~medians.index.isin(self.samplesToRemovePCA)]

            # mean-center data
            medians = medians - medians.mean(axis=0)

            # sort medians index (sample names) naturally
            medians = medians.reindex(natsorted(medians.index))

            # initialize PCA
            pca = PCA_MODULE(self.dimensionPCA, random_state=1)

            # fit PCA to data
            projected = pca.fit_transform(medians)

            # generate dataframe of the projection coordinates
            # to be used as plot input
            scatter_input = pd.DataFrame(data=projected, index=medians.index)
            scatter_input.rename(columns={0: 'PC1', 1: 'PC2'}, inplace=True)

            # assign row index (sample names) as a column
            scatter_input.reset_index(drop=False, inplace=True)

            # get sample file names (sampleMetadata keys) from config.yml
            # based on "Sample" column (element 1 of sampleMetadata vals)
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

            # create a naturally-sorted list of samples NOT be grayed
            ordered_conditions = natsorted(
                set(scatter_input.index.unique()).difference(
                    set(self.conditionsToSilhouette))
            )

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
            for e, (condition_name, sample_scores) in enumerate(scatter_input.iterrows()):

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
                    alpha=1.0, legend=False
                )

            g.grid(color='gray', linewidth=0.05, linestyle='-', alpha=1.0)
            plt.setp(g.spines.values(), color='k', lw=0.5)

            # assign row index (condition abbreviations) as column
            scatter_input = scatter_input.reset_index().rename(columns={'index': 'abbreviation'})

            # annotate data points
            if self.labelPoints is True:

                # generate squareform distance matrix
                sq = squareform(pdist(scatter_input[['PC1', 'PC2']], metric='euclidean'))

                # add numerical row and column indices
                df = pd.DataFrame(sq, index=scatter_input.index, columns=scatter_input.index)

                # isolate values from upper triangle
                df1 = df.where(np.triu(np.ones(df.shape)).astype('bool'))

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
                    list(df5['sample_id1']) + list(df5['sample_id2'])
                )

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
                                    list(df6['sample_id1']) + list(df6['sample_id2'])
                                )

                                # add neighboring indices to overall
                                # neighbors_set
                                neighbors_set = neighbors_set.union(neighbors)

                                # slice scatter_input to get samples
                                # proximal to data point e

                                neighbors_df = scatter_input.loc[list(neighbors)]

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
                                [path_effects.Stroke(linewidth=0.75, foreground='k'),
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
            g.legend(
                handles=legend_handles, prop={'size': 5.0},
                loc='upper left', bbox_to_anchor=[1.01, 1.0]
            )

            # update x and y axis labels
            plt.xlabel(
                'PC1 '
                f'({round((pca.explained_variance_ratio_[0] * 100), 2)}'
                '% of variance)', fontsize=10, labelpad=7.0)
            plt.ylabel(
                'PC2 '
                f'({round((pca.explained_variance_ratio_[1] * 100), 2)}'
                '% of variance)', fontsize=10, labelpad=4.0)

            # modify x and y axis ticks
            plt.tick_params(axis='both', which='major', labelsize=7.0)

            # save figure
            plt.savefig(
                os.path.join(pca_dir, 'pcaScoresPlotSamples.pdf'),
                bbox_inches='tight')
            plt.close('all')

            data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)
    else:
        logging.info("n_components = 1, skipping PCA and Horn's parallel analysis.")

    print()
    print()
    return data
