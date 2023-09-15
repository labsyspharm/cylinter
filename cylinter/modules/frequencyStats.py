import os
import logging

import numpy as np
import pandas as pd

import math
import natsort
from natsort import natsorted
from itertools import product

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from scipy.stats import ttest_ind

from ..utils import (
    input_check, read_markers, categorical_cmap, fdrcorrection, reorganize_dfcolumns
)

logger = logging.getLogger(__name__)


def frequencyStats(data, self, args):

    print()
    
    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
    )

    for type in ['class', f'cluster_{self.dimensionEmbedding}d']:
        if type in data.columns:

            stats_input = data[['Sample', 'Replicate', type]]

            # loop over comma-delimited binary declarations
            for i in range(len(list(self.sampleStatuses.values())[0].split(', '))):

                # get unique declaration categories (should be 2 per test)
                comparison = set(
                    [j.split(', ')[i] for j in self.sampleStatuses.values()
                     if '-UNK' not in j.split(', ')[i]])

                if len(comparison) > 1:

                    # assign test and control groups
                    test = [
                        i for i in comparison if i not in
                        self.controlGroups][0]
                    control = [
                        i for i in comparison if i in self.controlGroups][0]

                    # create frequency stats directory if it hasn't already
                    frequency_dir = os.path.join(
                        self.outDir, 'clustering',
                        f'{self.dimensionEmbedding}d',
                        'frequency_stats', type, f'{test}_v_{control}'
                    )
                    if not os.path.exists(frequency_dir):
                        os.makedirs(frequency_dir)

                    # create single-column dataFrame with all sample names
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

                    # loop over populations
                    for cluster, group in natsorted(stats_input.groupby(type)):
                            
                        if cluster not in [-1, 'unclassified']:
                            
                            logger.info(
                                f'Calculating log2({test}/{control}) of mean cell '
                                f'density for cluster {str(cluster)}.')

                            group = (
                                group.groupby(['Sample', 'Replicate', type])
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
                                0 if np.isnan(i) else int(i) for
                                i in group['count']
                            ]

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
                            group['Replicate'] = [self.sampleReplicates[i] for i in file_names]

                            group[type] = cluster

                            # drop samples for which a declaration cannot be made
                            group = group[~group['status'].str.contains('-UNK')]

                            group.reset_index(drop=True, inplace=True)

                            # get denominator cell count for each sample
                            if self.denominatorCluster is None:
                                group['tissue_count'] = [
                                    len(stats_input[stats_input['Sample'] == i])
                                    for i in group['Sample']]
                            else:
                                group['tissue_count'] = [
                                    len(stats_input[(stats_input['Sample'] == i) &
                                        (stats_input[type] == self.denominatorCluster)])
                                    for i in group['Sample']]

                            # compute density of cells per sample
                            group['density'] = group['count'] / group['tissue_count']

                            # append group data to catplot_input
                            catplot_input = pd.concat([catplot_input, group], axis=0)

                            # isolate test and control group values
                            cnd1_values = group['density'][group['status'] == test]
                            cnd2_values = group['density'][group['status'] == control]

                            # perform Welch's t-test (equal_var=False)
                            stat, pval = ttest_ind(
                                cnd1_values, cnd2_values, axis=0, equal_var=False,
                                nan_policy='propagate'
                            )

                            # round resulting values
                            stat = round(stat, 3)
                            pval = round(pval, 3)

                            # compute mean of test and control group values
                            cnd1_mean = np.mean(cnd1_values)
                            cnd2_mean = np.mean(cnd2_values)

                            # compute mean ratio
                            ratio = np.log2(
                                (cnd1_mean + 0.00000000001) / (cnd2_mean + 0.00000000001)
                            )

                            # compute mean difference
                            dif = cnd1_mean - cnd2_mean

                            cluster_list.append(cluster)
                            ratio_list.append(ratio)
                            dif_list.append(dif)
                            pval_list.append(pval)

                    # create stats dataframe
                    statistics = pd.DataFrame(
                        list(zip(cluster_list, ratio_list, dif_list, pval_list)),
                        columns=[type, 'ratio', 'dif', 'pval']
                    ).sort_values(by=type)

                    # compute FDR p-val corrections
                    # (uses statsmodels.stats.multitest implementation)
                    rejected, p_adjust = fdrcorrection(
                        statistics['pval'].tolist(), alpha=0.05, method='indep', is_sorted=False
                    )

                    statistics['qval'] = p_adjust

                    # save total stats table
                    statistics.to_csv(
                        os.path.join(frequency_dir, 'stats_total.csv'), index=False
                    )

                    if self.FDRCorrection:
                        stat = 'qval'
                    else:
                        stat = 'pval'

                    # isolate statistically significant stat values
                    significant = statistics[statistics[stat] <= 0.05].sort_values(by=stat)

                    # save significant stats table
                    significant.to_csv(
                        os.path.join(frequency_dir, 'stats_sig.csv'), index=False
                    )

                    # plot
                    sns.set_style('whitegrid')
                    fig, ax = plt.subplots()
                    plt.scatter(abs(significant['dif']), significant['ratio'], s=9.0, c='tab:red')

                    for label, qval, x, y in zip(
                        significant[type], significant[stat],
                            abs(significant['dif']), significant['ratio']):

                        plt.annotate(
                            (label, f'{stat[0]}=' + str(qval)), size=3,
                            xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom',
                            bbox=dict(boxstyle='round,pad=0.1', fc='yellow',
                                      alpha=0.0)
                        )

                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(8)
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(8)

                    plt.title(f'{test} vs. {control} ({stat[0]}<0.05)', fontsize=9)
                    plt.xlabel(f'abs({test} - {control})', fontsize=8)
                    plt.ylabel(f'log2({test} / {control})', fontsize=8)
                    plt.savefig(os.path.join(frequency_dir, 'plot.pdf'))
                    plt.close()

                    catplot_input.reset_index(drop=True, inplace=True)

                    catplot_input[stat] = [
                        'ns' if i not in
                        significant[type].unique() else
                        significant[stat][
                            significant[type] == i].values[0]
                        for i in catplot_input[type]]

                    # filter catplot_input to plot only significant differences
                    catplot_input = catplot_input[catplot_input[stat] != 'ns']

                    if not catplot_input.empty:
                        # build cmap
                        cmap = categorical_cmap(
                            numUniqueSamples=len(catplot_input['Sample'].unique()),
                            numCatagories=10, cmap='tab10', continuous=False
                        )

                        sample_color_dict = dict(
                            zip(natsorted(catplot_input['Sample'].unique()),
                                cmap.colors))

                        catplot_input[type] = (
                            catplot_input[type].astype(str) +
                            f'; {stat} = ' + catplot_input[stat].astype(str)
                        )

                        catplot_input.sort_values(
                            by=[stat, 'status', 'density'], key=lambda x:
                            natsort.natsort_keygen(
                                alg=natsort.ns.LOCALE |
                                natsort.ns.IGNORECASE)(x), inplace=True
                        )

                        sns.set(font_scale=0.3)
                        sns.set_style('whitegrid')
                        ncols = 5
                        nrows = math.ceil(len(catplot_input[type].unique()) / ncols)

                        fig = plt.figure(figsize=(ncols + 2, nrows))

                        # grid specifications
                        gs = plt.GridSpec(nrows=nrows, ncols=ncols, figure=fig)

                        for (name, group), ax in zip(
                            catplot_input.groupby(type, sort=False),
                                product(range(nrows), range(ncols))):

                            ax = fig.add_subplot(gs[ax[0], ax[1]])
                            
                            g = sns.barplot(
                                data=group, x='status', y='density', hue='Sample', 
                                palette=sample_color_dict, width=0.8, lw=0.0, ax=ax
                            )
                            ax.grid(lw=0.5)
                            [x.set_linewidth(0.5) for x in ax.spines.values()]
                            
                            g.set_xticklabels(
                                [i.get_text().split('-')[1] for i in g.get_xticklabels()]
                            )
                            plt.tick_params(axis='x', pad=-3)
                            ax.set(xlabel=None)
                            plt.tick_params(axis='y', pad=-3)
                            ax.yaxis.labelpad = 2
                            ax.set_title(name, size=2, pad=2)
                            ax.legend_.remove()

                        plt.tight_layout()
                         
                        file_names = [
                            get_key(i) for i in natsorted(catplot_input['Sample'].unique())
                        ]

                        sample_conds = [self.sampleConditions[i] for i in file_names]

                        sample_abbrs = [self.sampleConditionAbbrs[i] for i in file_names]

                        cond_abbr = [f'{i}-{j}' for i, j in zip(sample_conds, sample_abbrs)]

                        handles_dict = dict(zip(
                            natsorted(catplot_input['Sample'].unique()), cond_abbr)
                        )

                        legend_handles = []
                        for k, v in handles_dict.items():
                            legend_handles.append(
                                Line2D([0], [0], marker='o', color='none',
                                       label=v, markerfacecolor=sample_color_dict[k],
                                       markeredgecolor='k', markeredgewidth=0.2,
                                       markersize=5.0)
                            )

                        fig.legend(
                            handles=legend_handles, prop={'size': 5.0}, loc='upper left',
                            bbox_to_anchor=[1.0, 1.0]
                        )

                        plt.savefig(
                            os.path.join(frequency_dir, 'catplot.pdf'), bbox_inches='tight'
                        )
                        plt.close('all')

                        print()

                else:
                    logger.info(
                        'Only one binary declaration ' +
                        f'class represented for {list(comparison)[0]}. ' +
                        'Statistics will not be computed.')
                    print()
        print()

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    return data
