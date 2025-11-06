import os
import sys
import yaml
import logging

from natsort import natsorted

from matplotlib.backends.qt_compat import QtWidgets
from qtpy.QtCore import QTimer

import napari
from magicgui import magicgui

from ..utils import (
    input_check, read_markers, marker_channel_number, single_channel_pyramid,
    get_filepath, reorganize_dfcolumns, sort_qc_report
)

logger = logging.getLogger(__name__)

channels_to_samples = {}
arbitrary_selection_toggle = False
sample_index = 1


def callback(self, viewer, channel, sample, data, initial_callback, next_widget, next_layout, arbitrary_widget, arbitrary_layout, qc_report, report_path):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # clear existing channels from Napari window if they exist
    viewer.layers.clear()
    
    # remove next_widget and arbitrary_widget docks and layout
    # attributes from Napari viewer
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

    # read segmentation outlines, add to Napari
    file_path = get_filepath(self, check, sample, 'SEG')
    seg, min, max = single_channel_pyramid(file_path, channel=0)
    viewer.add_image(
        seg, rgb=False, blending='additive', colormap='gray',
        visible=False, name='segmentation', contrast_limits=(min, max)
    )
    
    # read DNA1 channel
    file_path = get_filepath(self, check, sample, 'TIF')
    channel_number = marker_channel_number(
        self, markers, self.counterstainChannel
    )
    dna, min, max = single_channel_pyramid(file_path, channel=channel_number)
    viewer.add_image(
        dna, rgb=False, blending='additive', colormap='gray',
        name=self.counterstainChannel, contrast_limits=(min, max)
    )

    # read target antibody image
    if channel != self.counterstainChannel:
        channel_number = marker_channel_number(self, markers, channel)
        file_path = get_filepath(self, check, sample, 'TIF')
        img, min, max = single_channel_pyramid(
            file_path, channel=channel_number
        )
        viewer.add_image(
            img, rgb=False, blending='additive', colormap='green',
            visible=True, name=channel, contrast_limits=(min, max)
        )

    # apply previously defined contrast limits if they exist 
    try:
        viewer.layers[self.counterstainChannel].contrast_limits = (
            qc_report['setContrast'][
                self.counterstainChannel][0], qc_report['setContrast'][
                self.counterstainChannel][1]
        )
    except KeyError:
        pass

    try:
        viewer.layers[channel].contrast_limits = (
            qc_report['setContrast'][channel][0], qc_report['setContrast'][
                channel][1])
    except KeyError:
        pass

    # dock (or re-dock) next_widget and arbitrary_widget to Napari window
    viewer.window.add_dock_widget(
        next_widget, name=f'Channel: {channel},  Sample: {sample}',
        area='right'
    )
    viewer.window.add_dock_widget(
        arbitrary_widget, name='Sample Selector', area='right'
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
        qc_report['setContrast'][self.counterstainChannel] = (
            [int(i) for i in 
             viewer.layers[self.counterstainChannel].contrast_limits]
        )
        qc_report['setContrast'][channel] = [
            int(i) for i in viewer.layers[channel].contrast_limits
        ]

        # sort and dump updated qc_report to YAML file
        qc_report_sorted = sort_qc_report(
            qc_report, module='setContrast', 
            order=[self.counterstainChannel] + abx_channels
        )
        f = open(report_path, 'w')
        yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
        
        # go to next sample
        try:
            if arbitrary_selection_toggle:
                sample_index -= 1 

            channel = list(channels_to_samples.keys())[sample_index]
            sample = channels_to_samples[
                list(channels_to_samples.keys())[sample_index]] 
            
            initial_callback = False
            callback(
                self, viewer, channel, sample, data, initial_callback, 
                next_widget, next_layout, arbitrary_widget, arbitrary_layout,
                qc_report, report_path
            )

            sample_index += 1
            arbitrary_selection_toggle = False
        
        except IndexError:

            print()
            logger.info('Contrast Adjustments Complete!')
            QTimer().singleShot(0, viewer.close)
    
    next_sample.native.setSizePolicy(
        QtWidgets.QSizePolicy.Fixed,
        QtWidgets.QSizePolicy.Fixed,
    )

    # give next_sample access to channel passed to callback
    next_sample.channel.bind(channel)

    next_layout.addWidget(next_sample.native)
    
    #######################################################################

    @magicgui(layout='vertical', call_button='Enter', 
              sample={'choices': list(natsorted(data['Sample'].unique())),
                      'label': 'Go to Sample'})
    def sample_selector(sample: str):

        return sample

    sample_selector.native.setSizePolicy(
        QtWidgets.QSizePolicy.Fixed,
        QtWidgets.QSizePolicy.Fixed
    )
    
    # Set the current sample value immediately after creating the widget
    sample_selector.sample.value = sample

    arbitrary_layout.addWidget(sample_selector.native)

    # call connect
    @sample_selector.called.connect
    def sample_callback(value: str):

        global arbitrary_selection_toggle
        
        sample = value

        print()
        if sample not in data['Sample'].unique():
            napari.utils.notifications.show_warning(
                'Sample name not in filtered data.'
            )
            pass
        else:
            # update channel contrast yaml with selected constrast limits
                
            qc_report['setContrast'][self.counterstainChannel] = (
                [int(i) for i in 
                 viewer.layers[self.counterstainChannel].contrast_limits]
            )
            qc_report['setContrast'][channel] = [
                int(i) for i in viewer.layers[channel].contrast_limits
            ]

            # dump updated qc_report to YAML file
            qc_report_sorted = sort_qc_report(
                qc_report, module='setContrast', 
                order=[self.counterstainChannel] + abx_channels
            )
            f = open(report_path, 'w')
            yaml.dump(
                qc_report_sorted, f, sort_keys=False, allow_unicode=False
            )

            initial_callback = False
            callback(
                self, viewer, channel, sample, data, initial_callback, 
                next_widget, next_layout, arbitrary_widget, arbitrary_layout,
                qc_report, report_path
            )

            arbitrary_selection_toggle = True
    
    #######################################################################
    # napari.utils.notifications.show_info(
    #     f'Viewing marker {channel} in sample {sample}'
    # )


# main
def setContrast(data, self, args):

    global channels_to_samples

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
        reload_report = False
        if qc_report is None:
            qc_report = {}
            reload_report = True
        if 'setContrast' not in qc_report or qc_report['setContrast'] is None:
            qc_report['setContrast'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(
                qc_report, module='setContrast', order=None
            )
            f = open(report_path, 'w')
            yaml.dump(
                qc_report_sorted, f, sort_keys=False, allow_unicode=False
            )
            qc_report = yaml.safe_load(open(report_path))
    except:
        logger.info(
            'Aborting; QC report missing from CyLinter output directory. '
            'Re-start pipeline from aggregateData module to initialize '
            'QC report.'
        )
        sys.exit()

    viewer = napari.Viewer(title='CyLinter')

    # generate next sample selection Qt widget
    next_widget = QtWidgets.QWidget()
    next_layout = QtWidgets.QVBoxLayout(next_widget)
    next_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Fixed,
        QtWidgets.QSizePolicy.Fixed,
    )

    # generate arbitrary sample selection Qt widget
    arbitrary_widget = QtWidgets.QWidget()
    arbitrary_layout = QtWidgets.QVBoxLayout(arbitrary_widget)
    arbitrary_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Fixed,
        QtWidgets.QSizePolicy.Fixed,
    )

    # identify samples with 85th percentile of median cell signal intensity  
    # (trying to avoid outliers associated with max values)
    for ch in [self.counterstainChannel] + abx_channels:
        medians = data[['Sample', ch]].groupby('Sample').median()
        percentile_value = medians.quantile(0.85).item()
        differences = abs(medians - percentile_value)
        # select sample whose median channel value is closest to quantile
        selected_sample = differences.idxmin().item()  
        channels_to_samples[ch] = selected_sample

    # pass first channel and sample in channels_to_samples to callback
    channel = list(channels_to_samples.keys())[0]
    sample = channels_to_samples[channel] 
    
    initial_callback = True
    callback(
        self, viewer, channel, sample, data, initial_callback, 
        next_widget, next_layout, arbitrary_widget, arbitrary_layout,
        qc_report, report_path
    )
    
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = 'um'

    napari.run()

    print()

    ###########################################################################
    # print current channel contrast limits and exit
     
    if (
        set(list(qc_report['setContrast'].keys())) == 
        set(abx_channels + [self.counterstainChannel])
         ):
        logger.info('Current channel contrast settings are as follows:')
        for k, v in qc_report['setContrast'].items():
            logger.info(f'{k}: {v}')
    else:
        logger.info(
            'Aborting; QC report does not contain contrast settings '
            'for all channels. Please ensure limits are selected for '
            'all channels.'
        )
        sys.exit()
    
    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
