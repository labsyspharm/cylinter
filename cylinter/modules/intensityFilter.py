import os
import sys
import logging
import pickle

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
    input_check, read_markers, napari_notification, 
    single_channel_pyramid, get_filepath, reorganize_dfcolumns, 
)

logger = logging.getLogger(__name__)

arbitrary_selection_toggle = False
sample_index = 1


def callback(self, viewer, sample, samples, data, initial_callback, selection_widget, selection_layout, hist_widget, hist_layout, intensity_dir): 

    if sample in data['Sample'].unique():
        
        print()
        
        check, markers_filepath = input_check(self)

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
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
        dna, min, max = single_channel_pyramid(file_path, channel=0)
        viewer.add_image(
            dna, rgb=False, blending='additive',
            name=f'{dna1}: {sample}', contrast_limits=(min, max)
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

        fig.suptitle(f'Sample={sample} mean DNA intensity', size=10)

        # get axis object from canvas
        ax = canvas.figure.subplots()

        # grab sample data
        group = data[data['Sample'] == sample].copy()
        group[dna1] = group[dna1] + 0.00000000001  # avoiding log(0) errors
        
        n, bins, patches = ax.hist(
            np.log(group[dna1]), bins=self.numBinsIntensity,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype='stepfilled',
            range=None, label='before'
        )

        ax.set_ylabel('Count')

        # add sliders to plot
        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = fig.add_axes(
            [0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = fig.add_axes(
            [0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        # specify data range
        rnge = [bins.min(), bins.max()]

        # load cutoffs dictionary if it exists
        if os.path.exists(os.path.join(intensity_dir, 'cutoffs.pkl')):
            f = open(os.path.join(intensity_dir, 'cutoffs.pkl'), 'rb')
            cutoffs_dict = pickle.load(f)
        else:
            # create cutoffs dictionary
            cutoffs_dict = {}
            
        try:
            lowerCutoff, upperCutoff = cutoffs_dict[sample]
            vbars = True  # toggle to show vertical red/blue bars on plot
            if lowerCutoff == upperCutoff:
                vbars = False  # cutoffs were negated
        except KeyError:
            lowerCutoff, upperCutoff = (0.0, 0.0)
            vbars = False

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

            napari_notification(f'Sliders updated to ({lowerCutoff:.3f}, {upperCutoff:.3f})')
            
            return lowerCutoff, upperCutoff

        # update sliders when moved
        sLower.on_changed(update)
        sUpper.on_changed(update)
        
        # add vbars to plot
        if vbars:
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
                (np.log(group[dna1]) > lowerCutoff) & (np.log(group[dna1]) < upperCutoff)
            ]

            # isolate x, y coordinates of selected centroids
            centroids = group_filtered[['Y_centroid', 'X_centroid']]

            # isolate cycle 1 DNA intensity values and assign
            # as quantitative point properties
            dna_intensity = np.log(group_filtered[dna1]).values
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

        # maintain reference to button after exiting callback()
        button_ax._button = button
        
        # dock (or re-dock) hist_widget to Napari window 
        viewer.window.add_dock_widget(
            hist_widget, name='DNA Intensity Histogram', area='right'
        )

        # remove and re-dock selection_widget if it exists 
        # so hist_widget appears first in Napari window
        if not initial_callback:
            viewer.window.remove_dock_widget(selection_widget)
            viewer.window.add_dock_widget(
                selection_widget, name='Arbitrary Sample Selection', area='right'
            )
        
        ###########################################################################
        
        @magicgui(
            layout='horizontal',
            call_button=' Apply Gates and Move to Next Sample -->'
        )
        def next_sample(sample):

            global arbitrary_selection_toggle
            global sample_index
            
            # get current cutoffs
            lowerCutoff, upperCutoff = update(val=None)

            # add cutoffs to dictionary and store
            cutoffs_dict[sample] = (lowerCutoff, upperCutoff)
            f = open(os.path.join(intensity_dir, 'cutoffs.pkl'), 'wb')
            pickle.dump(cutoffs_dict, f)
            f.close()

            # go to next sample
            try:
                if arbitrary_selection_toggle:
                    sample_index -= 1 

                sample = samples[sample_index]
                
                initial_callback = False
                callback(
                    self, viewer, sample, samples, data, initial_callback,
                    selection_widget, selection_layout, hist_widget, hist_layout,
                    intensity_dir 
                )

                sample_index += 1
                arbitrary_selection_toggle = False
            
            except IndexError:

                print()
                napari_notification('Gating complete!')
                QTimer().singleShot(0, viewer.close)

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
            QtWidgets.QSizePolicy.Minimum,
            QtWidgets.QSizePolicy.Maximum
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
                self, viewer, sample, samples, data, initial_callback,
                selection_widget, selection_layout, hist_widget, hist_layout,
                intensity_dir
            )

            arbitrary_selection_toggle = True
        
        ###########################################################################
        napari_notification(f'Working on Sample {sample}')
        
    else:
        print()
        napari_notification('Invalid entry.')
        pass


# main
def intensityFilter(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
    )

    # create intensity directory if it doesn't already exist
    intensity_dir = os.path.join(self.outDir, 'intensity')
    if not os.path.exists(intensity_dir):
        os.makedirs(intensity_dir)

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
    
    # make a list of samples, select the first one, and pass it to the callback
    samples = natsorted(data['Sample'].unique())
    sample = samples[0] 
    
    initial_callback = True
    callback(
        self, viewer, sample, samples, data, initial_callback,
        selection_widget, selection_layout, hist_widget, hist_layout,
        intensity_dir
    )
    
    viewer.window.add_dock_widget(
        selection_widget, name='Arbitrary Sample Selection', area='right'
    )
    
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = 'um'

    napari.run()
    
    print()
    
    ###########################################################################
    # load cutoffs dictionary if it exists
    
    if os.path.exists(os.path.join(intensity_dir, 'cutoffs.pkl')):
        f = open(os.path.join(intensity_dir, 'cutoffs.pkl'), 'rb')
        cutoffs_dict = pickle.load(f)

    else:
        print()
        logger.info(
            'Aborting; DNA intensity cutoffs dictionary does not exist. '
            'Please re-run intensityFilter module to select cutoffs.'
        )
        sys.exit()
    
    # save plots of selected data points
    plot_dir = os.path.join(intensity_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    
    idxs_to_drop = {}
    for sample in samples:

        group = data[data['Sample'] == sample]
        
        try:
            lowerCutoff, upperCutoff = cutoffs_dict[sample]
            if lowerCutoff == upperCutoff:
                logger.info(f'All data points selected for sample {sample}.')
            else:
                logger.info(
                    f'Applying cutoffs ({lowerCutoff:.3f}, '
                    f'{upperCutoff:.3f}) to sample {sample}'
                )
        except KeyError:
            print()
            logger.info(
                f'Aborting; Cutoffs have not been ' 
                f'selected for sample {sample}. '
                'Please re-run intensityFilter module to select '
                'cutoffs for this sample.'
            )
            sys.exit()

        # plot DNA intensity histogram BEFORE filtering
        fig, ax = plt.subplots()

        n, bins, patches = plt.hist(
            np.log(group[dna1]), bins=self.numBinsIntensity,
            density=False, color='b', ec='none',
            alpha=0.5, histtype='stepfilled',
            range=None, label='before'
        )
        
        # select all data points if sliders were not adjusted
        if lowerCutoff == upperCutoff:
            lowerCutoff = np.log(group[dna1]).min()
            upperCutoff = np.log(group[dna1]).max()
        
        # apply lower and upper cutoffs
        group_filtered = group.copy()[
            (np.log(group[dna1]) > lowerCutoff) & (np.log(group[dna1]) < upperCutoff)]

        # plot DNA intensity histogram AFTER filtering
        plt.hist(
            np.log(group_filtered[dna1]), bins=bins,
            density=False, color='r', ec='none', 
            alpha=0.5, histtype='stepfilled',
            range=None, label='after'
        )
        plt.xlabel('Mean DNA intensity')
        plt.ylabel('Count')
        plt.title(f'Sample={sample} mean DNA intensity', size=10)

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
            (np.log(group[dna1]) < lowerCutoff) | (np.log(group[dna1]) > upperCutoff)]

        # create a column of unique IDs for cells to drop from current sample
        data_to_drop['handle'] = data_to_drop['CellID'].map(str) + '_' + data_to_drop['Sample']

        # add IDs to idxs_to_drop dictionary
        # sidxs_to_drop[sample] = [i for i in data_to_drop['handle']]
        idxs_to_drop[sample] = data_to_drop['handle']

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
