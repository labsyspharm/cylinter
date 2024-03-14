import os
import logging

import math

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from ..utils import input_check, read_markers, reorganize_dfcolumns

logger = logging.getLogger(__name__)


def clustermap(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # create clustering dimension directory if it hasn't already
    dim_dir = os.path.join(self.outDir, 'clustering', f'{self.dimensionEmbedding}d')
    if not os.path.exists(dim_dir):
        os.makedirs(dim_dir)

    # drop antibody channel exclusions for clustering
    abx_channels = [i for i in abx_channels if i not in self.channelExclusionsClustering]

    ######################################################################

    sns.set_style("whitegrid", {'axes.grid': False})
    gs = plt.GridSpec(len(abx_channels), 1)
    fig = plt.figure(figsize=(2, 7))

    ax_objs = []
    for i, channel in enumerate(abx_channels):

        # creating new axes object
        ax_objs.append(fig.add_subplot(gs[i:i + 1, 0:]))
  
        # plotting the distribution
        n, bins, patches = ax_objs[-1].hist(
            data[channel], bins=50, density=True, histtype='stepfilled',
            linewidth=2.0, ec='k', alpha=1.0, color='k'
        )

        # setting uniform x and y lims
        ax_objs[-1].set_xlim(0, 1)
        ax_objs[-1].set_ylim(0, math.ceil(n.max()) + 1)

        # make background transparent
        rect = ax_objs[-1].patch
        rect.set_alpha(0)

        # remove borders, axis ticks, and labels
        ax_objs[-1].set_yticklabels([])

        if i == len(abx_channels) - 1:
            ax_objs[-1].set_xlabel(
                'Intensity', fontsize=11, fontweight='normal', labelpad=10
            )
        else:
            ax_objs[-1].set_xticks([])
            ax_objs[-1].set_xticklabels([])

        ax_objs[-1].set_yticks([])

        spines = ['top', 'right', 'left']
        for s in spines:
            ax_objs[-1].spines[s].set_visible(False)

        ax_objs[-1].tick_params(axis='x', width=2)

        ax_objs[-1].text(-0.02, 0, channel, fontweight='normal', fontsize=8, ha='right')
  
    gs.update(hspace=0.3)
    plt.subplots_adjust(left=0.3, bottom=0.1, right=0.9, top=0.95)
    plt.savefig(os.path.join(dim_dir, 'ridgeplots.pdf'))
    plt.close('all')

    ##############################################################################################

    for type in [f'cluster_{self.dimensionEmbedding}d', 'class']:
        if type in data.columns:

            if type == f'cluster_{self.dimensionEmbedding}d':

                clustermap_input = data[data[type] != -1]

                # compute mean antibody signals for clusters
                clustermap_input = clustermap_input[abx_channels + [type]].groupby(type).mean()

            elif type == 'class':
                
                clustermap_input = data[data[type] != 'unclassified']
                
                # compute mean antibody signals for clusters
                clustermap_input = clustermap_input[abx_channels + [type]].groupby(type).mean()

            if len(clustermap_input) > 1:
                
                sns.set(font_scale=0.7)

                # Compute per channel z-scores across clusters
                clustermap_input = (
                    (clustermap_input - clustermap_input.mean()) / clustermap_input.std()
                )
                # assign NaNs (channels with no variation in signal) to 0
                clustermap_input[clustermap_input.isna()] = 0
                
                # Zero-center colorbar
                norm = TwoSlopeNorm(
                    vcenter=0, vmin=clustermap_input.min().min(),
                    vmax=clustermap_input.max().max()
                )

                g = sns.clustermap(
                    clustermap_input, cmap='coolwarm', standard_scale=None, square=False,
                    xticklabels=1, yticklabels=1, linewidth=0.0, cbar=True, norm=norm 
                )
                
                # g = sns.clustermap(
                #     clustermap_input, cmap='viridis', standard_scale=1, square=False,
                #     xticklabels=1, yticklabels=1, linewidth=0.0, cbar=True 
                # )

                g.fig.suptitle('channel_z-scores.pdf', y=0.995, fontsize=10)
                g.fig.set_size_inches(6.0, 6.0)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)

                plt.savefig(
                    os.path.join(dim_dir, f'{type}_channel_z-scores.pdf'), bbox_inches='tight'
                )
            else:
                logger.info(
                    f' {type} clustermap cannot be generated with only one cell population.'
                )

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
