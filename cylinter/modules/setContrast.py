import os
import sys
import yaml
import logging

from matplotlib.backends.qt_compat import QtWidgets
from qtpy.QtCore import QTimer

import napari
from magicgui import magicgui

from ..utils import (
    input_check, read_markers, napari_notification, marker_channel_number,
    single_channel_pyramid, get_filepath, reorganize_dfcolumns
)

logger = logging.getLogger(__name__)

channels_to_samples = {}
arbitrary_selection_toggle = False
sample_index = 1


def callback(self, viewer, channel, sample, data, initial_callback, next_widget, next_layout, arbitrary_widget, arbitrary_layout, contrast_dir):

    print()

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
    )

    # clear existing channels from Napari window if they exist
    viewer.layers.clear()
    
    # remove next_widget and arbitrary_widget docks and layout attributes from Napari viewer
    if not initial_callback:
        viewer.window.remove_dock_widget(next_widget)
        count = next_layout.count()
        for i in range(count - 1, -1, -1):
            item = next_layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.setParent(None)

        viewer.window.remove_dock_widget(arbitrary_widget)
        count = arbitrary_layout.count()
        for i in range(count - 1, -1, -1):
            item = arbitrary_layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.setParent(None)

    # loop over antibody channels (except target channel,
    # which will be added to the Napari viewer last) and add them to Napari viewer
    for ch in reversed(abx_channels):
        if ch != channel:
            channel_number = marker_channel_number(markers, ch)

            # read antibody image
            file_path = get_filepath(self, check, sample, 'TIF')
            img, min, max = single_channel_pyramid(file_path, channel=channel_number)
            viewer.add_image(
                img, rgb=False, blending='additive', colormap='green',
                visible=False, name=ch, contrast_limits=(min, max)
            )

    # read DNA1 channel
    file_path = get_filepath(self, check, sample, 'TIF')
    dna, min, max = single_channel_pyramid(file_path, channel=0)
    viewer.add_image(
        dna, rgb=False, blending='additive', colormap='gray',
        name=dna1, contrast_limits=(min, max)
    )

    # read target antibody image
    channel_number = marker_channel_number(markers, channel)
    file_path = get_filepath(self, check, sample, 'TIF')
    img, min, max = single_channel_pyramid(file_path, channel=channel_number)
    viewer.add_image(
        img, rgb=False, blending='additive', colormap='green',
        visible=True, name=channel, contrast_limits=(min, max)
    )
    
    # apply previously defined contrast limits if they exist
    if os.path.exists(os.path.join(contrast_dir, 'contrast_limits.yml')):

        contrast_limits = yaml.safe_load(open(f'{contrast_dir}/contrast_limits.yml'))   

        for ch in list(contrast_limits.keys()):
            viewer.layers[ch].contrast_limits = (
                contrast_limits[ch][0], contrast_limits[ch][1])

    # dock (or re-dock) next_widget and arbitrary_widget to Napari window
    viewer.window.add_dock_widget(next_widget, name=f'Sample: {sample}', area='right')
    viewer.window.add_dock_widget(
        arbitrary_widget, name='Arbitrary Sample Selection', area='right'
    )

    #######################################################################
    
    @magicgui(
        layout='horizontal',
        call_button='Apply Limits and Move to Next Channel -->'
    )
    def next_sample(channel):

        global channels_to_samples
        global arbitrary_selection_toggle
        global sample_index

        # update channel contrast yaml with selected constrast limits
        if os.path.exists(os.path.join(contrast_dir, 'contrast_limits.yml')):

            contrast_limits = yaml.safe_load(open(f'{contrast_dir}/contrast_limits.yml'))
        else:
            contrast_limits = {}
            
        contrast_limits[dna1] = [int(i) for i in viewer.layers[dna1].contrast_limits]
        contrast_limits[channel] = [int(i) for i in viewer.layers[channel].contrast_limits]

        with open(f'{contrast_dir}/contrast_limits.yml', 'w') as file:
            yaml.dump(contrast_limits, file)
        
        # go to next sample
        try:
            if arbitrary_selection_toggle:
                sample_index -= 1 

            channel = list(channels_to_samples.keys())[sample_index]
            sample = channels_to_samples[list(channels_to_samples.keys())[sample_index]] 
            
            initial_callback = False
            callback(
                self, viewer, channel, sample, data, initial_callback, 
                next_widget, next_layout, arbitrary_widget, arbitrary_layout,
                contrast_dir
            )

            sample_index += 1
            arbitrary_selection_toggle = False
        
        except IndexError:

            print()
            napari_notification('Contrast Adjustments Complete!')
            QTimer().singleShot(0, viewer.close)
    
    next_sample.native.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Maximum,
    )

    # give next_sample access to channel passed to callback
    next_sample.channel.bind(channel)
    
    next_layout.addWidget(next_sample.native)
    
    #######################################################################

    @magicgui(layout='vertical', call_button='Enter', sample={'label': 'Sample Name'})
    def sample_selector(sample: str):

        return sample

    sample_selector.native.setSizePolicy(
        QtWidgets.QSizePolicy.Minimum,
        QtWidgets.QSizePolicy.Maximum
    )

    arbitrary_layout.addWidget(sample_selector.native)

    # call connect
    @sample_selector.called.connect
    def sample_callback(value: str):

        global arbitrary_selection_toggle
        
        sample = value

        # update channel contrast yaml with selected constrast limits
        if os.path.exists(os.path.join(contrast_dir, 'contrast_limits.yml')):

            contrast_limits = yaml.safe_load(
                open(f'{contrast_dir}/contrast_limits.yml'))
        else:
            contrast_limits = {}
            
        contrast_limits[dna1] = [int(i) for i in viewer.layers[dna1].contrast_limits]
        contrast_limits[channel] = [int(i) for i in viewer.layers[channel].contrast_limits]

        with open(f'{contrast_dir}/contrast_limits.yml', 'w') as file:
            yaml.dump(contrast_limits, file)

        initial_callback = False
        callback(
            self, viewer, channel, sample, data, initial_callback, 
            next_widget, next_layout, arbitrary_widget, arbitrary_layout,
            contrast_dir
        )

        arbitrary_selection_toggle = True
    
    #######################################################################
    napari_notification(f'Viewing channel {channel} in sample {sample}')


# main
def setContrast(data, self, args):

    global channels_to_samples

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
    )

    # create contrast directory if it hasn't already
    contrast_dir = os.path.join(self.outDir, 'contrast')
    if not os.path.exists(contrast_dir):
        os.makedirs(contrast_dir)

    viewer = napari.Viewer(title='CyLinter')

    # generate next sample selection Qt widget
    next_widget = QtWidgets.QWidget()
    next_layout = QtWidgets.QVBoxLayout(next_widget)
    next_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Fixed,
    )

    # generate arbitrary sample selection Qt widget
    arbitrary_widget = QtWidgets.QWidget()
    arbitrary_layout = QtWidgets.QVBoxLayout(arbitrary_widget)
    arbitrary_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Fixed,
    )

    # identify samples with 85th percentile of median cell signal intensity  
    # (try to avoid outliers associated with max values)
    for ch in abx_channels:
        medians = data[['Sample', ch]].groupby('Sample').median()
        percentile_value = medians.quantile(0.85).item()
        sample = (medians[ch] == percentile_value).idxmax()
        channels_to_samples[ch] = sample

    # pass first channel and sample in channels_to_samples to callback
    channel = list(channels_to_samples.keys())[0]
    sample = channels_to_samples[list(channels_to_samples.keys())[0]] 
    
    initial_callback = True
    callback(
        self, viewer, channel, sample, data, initial_callback, 
        next_widget, next_layout, arbitrary_widget, arbitrary_layout,
        contrast_dir
    )
    
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = 'um'

    napari.run()

    print()

    ##############################################################################################
    # print current channel contrast limits dictionary and exit

    if os.path.exists(os.path.join(contrast_dir, 'contrast_limits.yml')):
        
        contrast_limits = yaml.safe_load(open(f'{contrast_dir}/contrast_limits.yml'))
        
        if set(list(contrast_limits.keys())) == set(abx_channels + [dna1]):
            logger.info('Current channel contrast settings are as follows:')
            for k, v in contrast_limits.items():
                logger.info(f'{k}: {v}')
        else:
            logger.info(
                'Aborting; Contrast limits dictionary does not contain values for all channels. '
                'Please ensure limits are selected for all channels.'
            )
            sys.exit()
        
        data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

        print()
        print()
        return data
    
    else:
        logger.info(
            'Aborting; Contrast limits dictionary does not exist. '
            'Please select contrast limits.'
        )
        sys.exit()

    
