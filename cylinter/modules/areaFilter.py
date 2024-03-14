import os
import sys
import yaml
import logging

import pandas as pd
import numpy as np

from natsort import natsorted

import napari
from magicgui import magicgui

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.widgets import Slider, Button
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

from qtpy.QtCore import QTimer

from ..utils import (
    input_check, read_markers, marker_channel_number, sort_qc_report, 
    single_channel_pyramid, get_filepath, reorganize_dfcolumns, compute_gmm 
)

logger = logging.getLogger(__name__)

arbitrary_selection_toggle = False
sample_index = 1


def callback(self, viewer, sample, samples_to_run, data, initial_callback, selection_widget, selection_layout, hist_widget, hist_layout, area_dir, qc_report, report_path):

    if sample in data['Sample'].unique():
        
        print()
        
        # FOR ALI
        # data.drop(columns='Area', inplace=True)
        # data.rename(columns={'N_distance': 'Area'}, inplace=True)
        # data['Area'] = np.log(data['Area'])
        
        check, markers_filepath = input_check(self)

        # read marker metadata
        markers, abx_channels = read_markers( 
            markers_filepath=markers_filepath,
            counterstain_channel=self.counterstainChannel,
            markers_to_exclude=self.markersToExclude, data=None
        )

        # clear existing channels from Napari window if they exist
        viewer.layers.clear()

        # read segmentation outlines and add to Napari viewer
        file_path = get_filepath(self, check, sample, 'SEG')
        seg, min, max = single_channel_pyramid(file_path, channel=0)
        viewer.add_image(
            seg, rgb=False, visible=False, colormap='gray', opacity=0.5,
            name='segmentation', contrast_limits=(min, max)
        )

        # read DNA1 and add to Napari viewer
        file_path = get_filepath(self, check, sample, 'TIF')
        channel_number = marker_channel_number(self, markers, self.counterstainChannel)
        dna, min, max = single_channel_pyramid(file_path, channel=channel_number)
        viewer.add_image(
            dna, rgb=False, blending='additive',
            name=self.counterstainChannel, contrast_limits=(min, max)
        )
        
        # remove hist_widget and layout attributes from Napari viewer if they exist
        if not initial_callback:
            viewer.window.remove_dock_widget(hist_widget)

            count = hist_layout.count()
            for i in range(count - 1, -1, -1):
                item = hist_layout.itemAt(i)
                widget = item.widget()
                if widget:
                    widget.setParent(None)

        # generate a blank figure canvas
        canvas = FigureCanvas(Figure(figsize=(5.5, 4.5)))

        # add navigation tool bar and figure canvas to hist_layout
        hist_layout.addWidget(NavigationToolbar(canvas, hist_widget))
        hist_layout.addWidget(canvas)

        ###########################################################################
        # plot histogram
        
        sns.set_style('whitegrid')
        
        # get figure object from canvas
        fig = canvas.figure

        fig.subplots_adjust(left=0.25, bottom=0.25)

        fig.suptitle(f'Sample = {sample}', size=10)

        # get axis object from canvas
        ax = canvas.figure.subplots()
        ax2 = ax.twinx()
        
        # grab sample data
        group = data[data['Sample'] == sample].copy()

        # avoiding log(0) errors
        group['Area'] = group['Area'] + 0.001  
        
        n, bins, patches = ax.hist(
            np.log(group['Area']), bins=self.numBinsArea,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype='stepfilled',
            range=None, label='before'
        )

        ax.set_ylabel('Count')
        ax2.set_ylabel('GMM density')

        # add sliders to plot
        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = fig.add_axes(
            [0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = fig.add_axes(
            [0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        # specify data range
        rnge = [bins.min(), bins.max()]

        # convert histogram data into 2D numpy array with 1 column to pass to GMM
        gmm_data = np.log(group['Area']).values.reshape(-1, 1)
        
        # compute GMM
        comp_lower_percentile, comp_upper_percentile = compute_gmm(
            data=gmm_data, x_min=rnge[0], x_max=rnge[1], ax=ax2
        )
        
        try:
            lowerCutoff, upperCutoff = qc_report['areaFilter'][sample]
            if lowerCutoff is None or upperCutoff is None:
                lowerCutoff, upperCutoff = (comp_lower_percentile, comp_upper_percentile)

        except (TypeError, KeyError, ValueError):
            # use GMM to assign default lower and upper thresholds
            lowerCutoff, upperCutoff = (comp_lower_percentile, comp_upper_percentile)

        # add slider functionality
        sLower = Slider(
            axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
            valinit=lowerCutoff, valstep=(rnge[1] / 100000))
        sLower.label.set_fontsize(11)
        sLower.label.set_color('b')
        
        sUpper = Slider(
            axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
            valinit=upperCutoff, valstep=(rnge[1] / 100000))
        sUpper.label.set_fontsize(11)
        sUpper.label.set_color('r')

        # function for updating sliders
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
        
        # add vbars to plot
        update(val=None)
        
        # add button to show selected centroids in Napari viewer
        button_ax = fig.add_axes([0.65, 0.025, 0.25, 0.06])
        button = Button(button_ax, 'Plot Points', color=axcolor, hovercolor='0.975')
        button.label.set_fontsize(11)

        def apply_cutoffs(event):

            # get current cutoffs
            lowerCutoff, upperCutoff = sLower.val, sUpper.val

            # apply lower and upper cutoffs
            group_filtered = group[
                (np.log(group['Area']) > lowerCutoff) & (np.log(group['Area']) < upperCutoff)
            ]

            # isolate x, y coordinates of selected centroids
            centroids = group_filtered[['Y_centroid', 'X_centroid']]

            # isolate segmentation area values and assign
            # as quantitative point properties
            cell_area = np.log(group_filtered['Area']).values
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

        # maintain reference to button after exiting callback()
        button_ax._button = button
        
        # dock (or re-dock) hist_widget to Napari window 
        viewer.window.add_dock_widget(
            hist_widget, name='Cell Segmentation Area Histogram', area='right'
        )

        # remove and re-dock selection_widget if it exists 
        # so hist_widget appears first in Napari window
        if not initial_callback:
            viewer.window.remove_dock_widget(selection_widget)
            viewer.window.add_dock_widget(
                selection_widget, name='Arbitrary Sample Selector', area='right'
            )
        
        ###########################################################################
        
        @magicgui(
            layout='horizontal',
            call_button='Apply Gates and Move to Next Sample -->'
        )
        def next_sample(sample):

            global arbitrary_selection_toggle
            global sample_index
            
            # get current cutoffs
            lowerCutoff, upperCutoff = update(val=None)

            if lowerCutoff <= upperCutoff:
           
                # store cutoffs in QC report
                qc_report['areaFilter'][sample] = [float(lowerCutoff), float(upperCutoff)]

                # sort and dump updated qc_report to YAML file
                qc_report_sorted = sort_qc_report(qc_report, module='areaFilter', order=None)
                f = open(report_path, 'w')
                yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
                
                napari.utils.notifications.show_info(
                    f'Sliders updated to ({lowerCutoff:.3f}, {upperCutoff:.3f})'
                )
                
                # go to next sample
                try:
                    if arbitrary_selection_toggle:
                        sample_index -= 1 

                    sample = samples_to_run[sample_index]
                    
                    initial_callback = False
                    callback(
                        self, viewer, sample, samples_to_run, data, initial_callback,
                        selection_widget, selection_layout, hist_widget, hist_layout,
                        area_dir, qc_report, report_path
                    )

                    sample_index += 1
                    arbitrary_selection_toggle = False
                
                except IndexError:

                    print()
                    logger.info('Thresholding complete!')
                    QTimer().singleShot(0, viewer.close)
            
            else:
                napari.utils.notifications.show_warning(
                    'LowerCutoff (blue) must be lower than upperCutoff (red).'
                )
                pass
        
        next_sample.native.setSizePolicy(
            QtWidgets.QSizePolicy.Maximum,
            QtWidgets.QSizePolicy.Maximum,
        )
        
        # give next_sample access to sample variable passed to callback
        next_sample.sample.bind(sample)

        hist_layout.addWidget(next_sample.native)
        
        ###########################################################################

        @magicgui(
            layout='vertical', call_button='Enter',
            sample={'label': 'Sample Name'}
        )
        def sample_selector(sample: str):

            return sample

        sample_selector.native.setSizePolicy(
            QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Fixed
        )
        
        if initial_callback:  
            selection_layout.addWidget(sample_selector.native)

        # call connect
        @sample_selector.called.connect
        def sample_callback(value: str):

            global arbitrary_selection_toggle

            sample = value

            initial_callback = False
            callback(
                self, viewer, sample, samples_to_run, data, initial_callback,
                selection_widget, selection_layout, hist_widget, hist_layout,
                area_dir, qc_report, report_path
            )

            arbitrary_selection_toggle = True
        
        ###########################################################################
        napari.utils.notifications.show_info(f'Viewing Sample {sample}')  
    
    else:
        print()
        napari.utils.notifications.show_warning(
            'Sample name not in filtered data.'
        )


# main
def areaFilter(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # create area directory if it doesn't already exist
    area_dir = os.path.join(self.outDir, 'area')
    if not os.path.exists(area_dir):
        os.makedirs(area_dir)

    # read QC report
    report_path = os.path.join(self.outDir, 'cylinter_report.yml')
    try:
        qc_report = yaml.safe_load(open(report_path))
        reload_report = False
        if qc_report is None:
            qc_report = {}
            reload_report = True
        if 'areaFilter' not in qc_report or qc_report['areaFilter'] is None:
            qc_report['areaFilter'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(qc_report, module='areaFilter', order=None)
            f = open(report_path, 'w')
            yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
            qc_report = yaml.safe_load(open(report_path))
    except:
        logger.info(
            'Aborting; QC report missing from CyLinter output directory. Re-start pipeline '
            'from aggregateData module to initialize QC report.'
        )
        sys.exit()

    # initialize Napari viewer
    viewer = napari.Viewer(title='CyLinter')

    # generate arbitrary sample selection Qt widget
    selection_widget = QtWidgets.QWidget()
    selection_layout = QtWidgets.QVBoxLayout(selection_widget)
    selection_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Minimum,
        QtWidgets.QSizePolicy.Fixed,
    )
    
    # generate histogram Qt widget
    hist_widget = QtWidgets.QWidget()
    hist_layout = QtWidgets.QVBoxLayout(hist_widget)
    hist_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Minimum,
        QtWidgets.QSizePolicy.Maximum
    )
    
    # make a list of all samples in batch, select the first one
    # for which cutoffs have not been previously assigned, and pass it to the callback
    samples = natsorted(data['Sample'].unique())
    
    if len(samples) == 0:
        print()
        logger.info(
            'Aborting; Single-cell data have not been selected for analysis. '
            'Re-run selectROIs module to select data points.'
        )
        sys.exit()
    
    else:
        samples_to_run = []
        for sample in samples:
            
            try:  
                lowerCutoff, upperCutoff = qc_report['areaFilter'][sample]
                if not isinstance(lowerCutoff, float) or not isinstance(upperCutoff, float):
                    samples_to_run.append(sample)
                
            except: 
                samples_to_run.append(sample)
    
    try:
        sample = samples_to_run[0]
        
        initial_callback = True
        callback(
            self, viewer, sample, samples_to_run, data, initial_callback,
            selection_widget, selection_layout, hist_widget, hist_layout,
            area_dir, qc_report, report_path
        )
        
        viewer.window.add_dock_widget(
            selection_widget, name='Arbitrary Sample Selector', area='right'
        )
        
        viewer.scale_bar.visible = True
        viewer.scale_bar.unit = 'um'

        napari.run()

    except IndexError:
        logger.info('Thresholding complete!')
        QTimer().singleShot(0, viewer.close)

    print()

    ###########################################################################
    # save histogram plots with selected data highlighted

    plot_dir = os.path.join(area_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    
    idxs_to_drop = {}
    for sample in samples_to_run:
        group = data[data['Sample'] == sample]
        
        lowerCutoff, upperCutoff = qc_report['areaFilter'][sample]

        # select all data points if sliders were not adjusted
        if lowerCutoff == upperCutoff:
            logger.info(f'All data points selected for sample {sample}.')
            lowerCutoff = np.log(group[self.counterstainChannel]).min()
            upperCutoff = np.log(group[self.counterstainChannel]).max()
        else:
            logger.info(
                f'Applying cutoffs ({lowerCutoff:.3f}, '
                f'{upperCutoff:.3f}) to sample {sample}'
            )

        # plot cell segmentation area histogram BEFORE filtering
        fig, ax = plt.subplots()

        n, bins, patches = plt.hist(
            np.log(group['Area']), bins=self.numBinsArea,
            density=False, color='b', ec='none',
            alpha=0.5, histtype='stepfilled',
            range=None, label='before'
        )

        # apply lower and upper cutoffs
        group_filtered = group.copy()[
            (np.log(group['Area']) > lowerCutoff) & (np.log(group['Area']) < upperCutoff)]

        # plot cell segmentation area histogram AFTER filtering
        plt.hist(
            np.log(group_filtered['Area']), bins=bins,
            density=False, color='r', ec='none', 
            alpha=0.5, histtype='stepfilled',
            range=None, label='after'
        )
        plt.xlabel('Cell Segementation Area')
        plt.ylabel('Count')
        plt.title(f'Sample = {sample}', size=10)

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
        plt.savefig(os.path.join(plot_dir, f'{sample}.pdf'))
        plt.close('all')

        # isolate sample data to drop
        data_to_drop = group.copy()[
            (np.log(group['Area']) < lowerCutoff) | (np.log(group['Area']) > upperCutoff)]

        if not data_to_drop.empty:
            # create a column of unique IDs for cells to drop from current sample
            data_to_drop['handle'] = (
                data_to_drop['CellID'].map(str) + '_' + data_to_drop['Sample']
            )

            # add IDs to idxs_to_drop dictionary
            idxs_to_drop[sample] = data_to_drop['handle']
        else:
            idxs_to_drop[sample] = pd.Series()
    
    # create a column of unique IDs for cells in the full dataframe
    data['handle'] = data['CellID'].map(str) + '_' + data['Sample']

    # create an overall list of indices to drop from the dataframe
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
