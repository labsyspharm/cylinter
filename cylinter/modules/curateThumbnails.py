import os
import sys
import ast
import yaml
import logging

import numpy as np
import pandas as pd

from natsort import natsorted

from subprocess import run

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from skimage.color import gray2rgb
from skimage.util.dtype import img_as_float

import zarr

from ..utils import (
    input_check, read_markers, cluster_expression, gate_expression,
    marker_channel_number, get_filepath, reorganize_dfcolumns
)

logger = logging.getLogger(__name__)


def curateThumbnails(data, self, args):

    print()

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # read QC report
    report_path = os.path.join(self.outDir, 'cylinter_report.yml')
    try:
        qc_report = yaml.safe_load(open(report_path))
    except FileNotFoundError:
        print()
        logger.info(
            'Aborting; QC report not found. Ensure cylinter_report.yml is '
            'stored at top-level of CyLinter output file or re-start pipeline '
            'to start filtering data.'
        )
        sys.exit()

    for type in ['class', f'cluster_{self.dimensionEmbedding}d']:

        df = data.copy()
        
        if type in df.columns:

            if type == f'cluster_{self.dimensionEmbedding}d':
                
                # drop antibody channel exclusions for clustering
                abx_channels = [
                    i for i in abx_channels if i 
                    not in self.channelExclusionsClustering
                ]

                # drop unclustered cells from data
                df = df[df[type] != -1]

            elif type == 'class':
                
                # drop antibody channel exclusions for gating
                abx_channels = [
                    i for i in abx_channels if i 
                    not in self.channelExclusionsGating
                ]

                # drop unclassified cells from data
                df = df[df[type] != 'unclassified']

            # create thumbnails directory
            thumbnails_dir = os.path.join(
                self.outDir, 'clustering', f'{self.dimensionEmbedding}d', 
                'thumbnails', type
            )
            if not os.path.exists(thumbnails_dir):
                os.makedirs(thumbnails_dir)

            # create zarr directory
            zarr_dir = os.path.join(thumbnails_dir, 'zarrs')
            if not os.path.exists(zarr_dir):
                os.makedirs(zarr_dir)

            ###################################################################

            # read the indices of clusters that have already been run
            if os.path.exists(os.path.join(thumbnails_dir, 'completed.txt')):

                with open(
                  os.path.join(thumbnails_dir, 'completed.txt'), 'r') as f:
                    completed = f.readlines()
                    completed = ast.literal_eval(completed[0])

                if type == f'cluster_{self.dimensionEmbedding}d':

                    completed = set([int(i) for i in completed])

                elif type == 'class':

                    completed = set([i for i in completed])

                total = set(df[type].unique())

                to_run = natsorted(total.difference(completed))

                # convert completed from a set to a list to append
                # to while looping over populations
                completed = list(completed)

                logger.info(
                    f'{len(to_run)} {type} image gallery(s) to generate.'
                )

            else:
                # create a list of populations to run
                completed = []
                to_run = natsorted(df[type].unique())

                logger.info(
                    f'{len(to_run)} {type} image gallery(s) to create.'
                )

            ###################################################################
            for pop in to_run:

                # create dataframe to collect thumbnail images and metadata
                long_table = pd.DataFrame()

                if type == f'cluster_{self.dimensionEmbedding}d':

                    # identify top expressed markers
                    hi_markers = cluster_expression(
                        df=df, markers=abx_channels, cluster=pop,
                        num_proteins=3,
                        clus_dim=self.dimensionEmbedding
                    )
                elif type == 'class':

                    # select non-negated markers
                    gate_dir = os.path.join(self.outDir, 'gating')
                    hi_markers = gate_expression(pop=pop, gate_dir=gate_dir)
                
                # combine DNA1 with top expressed markers
                markers_to_show = [self.counterstainChannel] + hi_markers

                # get marker channel numbers from markers.csv
                channel_nums = []
                for marker in markers_to_show:
                    channel_number = marker_channel_number(
                        self, markers, marker
                    )
                    # cellcutter 1-based indexing
                    channel_nums.append(str(channel_number + 1))   

                # create marker LUT
                color_dict = {}
                for i, j in zip(
                    markers_to_show,
                    [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0),
                     (0.0, 0.0, 1.0), (1.0, 1.0, 0), (1.0, 0.0, 1.0),
                     (0.0, 1.0, 1.0), (1.0, 0.5, 0.0), (0.5, 0.0, 1.0),
                     (0.0, 0.5, 1.0), (0.5, 1.0, 0.0), (1.0, 0.0, 0.5),
                     (0.0, 1.0, 0.5), (0.5, 0.5, 0.0), (0.0, 0.5, 0.5),
                     (0.5, 0.0, 0.0)]):

                    color_dict[i] = j

                for sample in [i for i in df['Sample'].unique()]:

                    # isolate data for current cluster and sample
                    cellcutter_input = df[
                        (df[type] == pop) & (df['Sample'] == sample)]

                    if len(cellcutter_input) == 0:
                        continue

                    # randomly select example cells
                    if len(cellcutter_input) > self.numThumbnails:

                        cellcutter_input = cellcutter_input.sample(
                            n=self.numThumbnails, random_state=1
                        )

                    # generate cellcutter zarr file if it doesn't already exist
                    # if not os.path.exists(
                    #    os.path.join(
                    #         zarr_dir, f'{type}_{pop}_sample_{sample}'
                    #         f'_win{self.windowSize}.zarr', '0.0.0.0')):

                    # write cellcutter_input to disk
                    cellcutter_input.to_csv(
                        os.path.join(thumbnails_dir, 'csv_data.csv'),
                        index=False
                    )

                    print()
                    
                    if type == f'cluster_{self.dimensionEmbedding}d':
                        logger.info(
                            f'Cutting {len(cellcutter_input)} '
                            f'{type} {pop} cell(s) '
                            f'showing {markers_to_show} channels '
                            f'from sample {sample}...'
                        )
                    else:
                        logger.info(
                            f'Cutting {len(cellcutter_input)} '
                            f'{type} {pop} cell(s) '
                            f'showing {markers_to_show} channels '
                            f'from sample {sample}...'
                        )

                    # run cellcutter on sample image
                    tif_file_path = get_filepath(self, check, sample, 'TIF')
                    mask_file_path = get_filepath(self, check, sample, 'MASK')
                    run(
                        ["cut_cells", "--force", "--window-size",
                         f"{self.windowSize}", "--cells-per-chunk",
                         "200", "--cache-size", "57711",
                         f"{tif_file_path}",
                         f"{mask_file_path}",
                         f"{thumbnails_dir}/csv_data.csv",
                         f"{zarr_dir}/{type}_{pop}_sample_{sample}"
                         f"_win{self.windowSize}.zarr",
                         "--channels"] + channel_nums,
                        check=True
                    )

                    # read multi-channel zarr file created by cellcutter
                    z_path_img = os.path.join(
                        zarr_dir, f'{type}_{pop}_sample_{sample}'
                        f'_win{self.windowSize}.zarr'
                    )

                    z_img = zarr.open(z_path_img, mode='r')

                    if self.segOutlines:

                        # if not os.path.exists(
                        #    os.path.join(
                        #         zarr_dir, f'{type}_{pop}_sample_{sample}'
                        #         f'_win{self.windowSize}_seg.zarr', '0.0.0.0')):

                        # write cellcutter_input to disk
                        cellcutter_input.to_csv(
                            os.path.join(thumbnails_dir, 'csv_data.csv'), 
                            index=False
                        )

                        print()
                        
                        if type == f'cluster_{self.dimensionEmbedding}d':
                            logger.info(
                                'Cutting segmentation outlines for '
                                f'{len(cellcutter_input)} '
                                f'{type} {pop} cell(s) '
                                f'showing {markers_to_show} channels '
                                f'from sample {sample}...'
                            )
                        else:
                            logger.info(
                                'Cutting segmentation outlines for '
                                f'{len(cellcutter_input)} '
                                f'{type} {pop} cell(s) '
                                f'showing {markers_to_show} channels '
                                f'from sample {sample}...'
                            )

                        # run cellcutter on segmentation outlines image
                        seg_file_path = get_filepath(
                            self, check, sample, 'SEG'
                        )
                        mask_file_path = get_filepath(
                            self, check, sample, 'MASK'
                        )
                        run(
                            ["cut_cells", "--force", "--window-size",
                             f"{self.windowSize}",
                             "--cells-per-chunk", "200",
                             "--cache-size", "57711",
                             f"{seg_file_path}",
                             f"{mask_file_path}",
                             f"{thumbnails_dir}/csv_data.csv",
                             f"{zarr_dir}/{type}_{pop}_sample_{sample}"
                             f"_win{self.windowSize}_seg.zarr",
                             "--channels", "1"],
                            check=True
                        )

                        # read segmentation outlines zarr file
                        # created by cellcutter
                        z_path_seg = os.path.join(
                            zarr_dir,
                            f'{type}_{pop}_sample_{sample}_win{self.windowSize}_seg.zarr'
                        )
                        z_seg = zarr.open(z_path_seg, mode='r')

                    if os.path.exists(
                      os.path.join(thumbnails_dir, 'csv_data.csv')):
                        # remove cellcutter_input file after cells are cut
                        os.remove(os.path.join(thumbnails_dir, 'csv_data.csv'))

                    # create composite thumbnail images
                    for cell in range(z_img.shape[1]):

                        # create blank image with same x/y dims as thumbnail
                        blank_img = np.zeros((z_img.shape[2], z_img.shape[3]))

                        # add centroid point at the center of the image
                        blank_img[
                            int(z_img.shape[2] / 2):int(z_img.shape[2] / 2) + 1,
                            int(z_img.shape[3] / 2):int(z_img.shape[3] / 2) + 1
                        ] = 1

                        # convert blank image to rgb and colorize
                        blank_img = gray2rgb(blank_img)

                        # loop over markers to show
                        for ch, marker in enumerate(markers_to_show):

                            # slice marker channel of current thumbnail
                            slice = img_as_float(z_img[ch, cell, :, :])

                            # apply image contrast settings
                            if z_img.dtype == 'uint16':
                                divisor = 65535
                            elif z_img.dtype == 'uint8':
                                divisor = 255
                            
                            # apply contrast settings if in QC report
                            try:  
                                slice -= (qc_report['setContrast'][marker][0] /
                                          divisor)
                                slice /= (
                                    (qc_report[
                                        'setContrast'][marker][1] / divisor) -
                                    (qc_report[
                                        'setContrast'][marker][0] / divisor)
                                )
                            except:
                                pass
                            
                            slice = np.clip(slice, 0, 1)

                            # convert channel slice to RGB, colorize,
                            # and add to blank image
                            slice = gray2rgb(slice)
                            slice = (slice * color_dict[marker])
                            blank_img += slice

                        if self.segOutlines:

                            # get segmentation thumbnails
                            seg = img_as_float(z_seg[0, cell, :, :])

                            # ensure segmentation outlines are normalized 0-1
                            seg = (seg - np.min(seg)) / np.ptp(seg)

                            # convert segmentation thumbnail to RGB
                            # and add to blank image
                            seg = gray2rgb(seg) * 0.25  # decrease alpha

                            blank_img += seg

                        # append merged thumbnail image to long_table
                        append_df = pd.DataFrame.from_dict(
                            {'sample': sample,
                             'example': int(cell + 1),
                             'image': blank_img}, orient='index').T
                        
                        long_table = pd.concat(
                            [long_table, append_df], axis=0, ignore_index=True
                        )

                # plot facet grid of thumbnails for current cluster
                if not long_table.empty:
                    fig, ax = plt.subplots()

                    g = sns.FacetGrid(
                        long_table, row='sample', col='example',
                        sharex=False, sharey=False,
                        gridspec_kws={'hspace': 0.1, 'wspace': 0.1}
                    )

                    g.map(
                        lambda image, **kwargs: (
                            plt.imshow(np.clip(image.values[0], 0, 1)), 
                            plt.grid(False)), 'image'
                    )  # image clipping prevents matplotlib warning

                    for ax in g.axes.flatten():
                        ax.get_xaxis().set_ticks([])
                        ax.set_xlabel('')
                        ax.get_yaxis().set_ticks([])
                        ax.set_ylabel('')

                    g.set_titles(
                        col_template="Ex. {col_name}", row_template="Smpl. {row_name}",
                        fontweight='bold', size=6
                    )

                    custom_lines = []
                    for k, v in color_dict.items():

                        custom_lines.append(Line2D([0], [0], color=v, lw=6))

                    ax.legend(
                        custom_lines,
                        list(color_dict.keys()), prop={'size': 12},
                        bbox_to_anchor=(
                            1.05, len(long_table['sample'].unique()) + 0.3),
                        loc='upper left'
                    )

                    plt.savefig(
                        os.path.join(thumbnails_dir, str(pop) + '.pdf'), 
                        bbox_inches='tight'
                    )

                    plt.close('all')

                # update completed clusters list
                completed.append(pop)

                # overwrite completed_clusters.txt file
                with open(
                  os.path.join(thumbnails_dir, 'completed.txt'), 'w') as f:
                    f.write(str(completed))

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    return data
