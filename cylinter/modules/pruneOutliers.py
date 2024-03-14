import os
import sys
import yaml
import logging

import numpy as np
import pandas as pd
import math

from natsort import natsorted

import napari
from magicgui import magicgui

from sklearn.preprocessing import MinMaxScaler

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

from qtpy.QtCore import QTimer

from ..utils import (
    input_check, read_markers, marker_channel_number, sort_qc_report, 
    single_channel_pyramid, get_filepath, reorganize_dfcolumns
)

logger = logging.getLogger(__name__)

arbitrary_selection_toggle = False
arbitrary_marker_idx = None
marker_index = 1
lowerCutoff = 0.0
upperCutoff = 100.0
dfTest = None


def keys_before_key(dictionary, abx_channels, target_channel):
    
    keys = sorted(dictionary.keys(), key=lambda x: abx_channels.index(x))
    index = keys.index(target_channel)
    subset_dict = {key: dictionary[key] for key in keys[:index] if key in dictionary}
    
    return subset_dict


def callback(self, viewer, channel, dfTrim, data, initial_callback, percentiles_widget,percentiles_layout, arbitrary_widget, arbitrary_layout, plot_widget, plot_layout, pruning_dir, plot_dir, qc_report, report_path):
    
    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    ) 

    napari.utils.notifications.show_info(f'Viewing data for marker {channel}')
    
    # clear existing channels from Napari window if they exist
    viewer.layers.clear()

    ##########################################################################################
    # plot raw signal distributions

    hist_facet = (
        dfTrim[['Sample', 'Condition', 'Area'] + [channel]]
        .sample(frac=1.0)
        .melt(id_vars=['Sample', 'Condition', 'Area'],
              var_name='channel', value_name='signal')
    )

    # naturally sort hist_facet by sample
    hist_facet['Sample'] = pd.Categorical(
        hist_facet['Sample'], ordered=True,
        categories=natsorted(hist_facet['Sample'].unique()))
    hist_facet.sort_values('Sample', inplace=True)
    hist_facet['Sample'] = hist_facet['Sample'].astype('str')

    # create column for facet labels
    hist_facet['for_plot'] = hist_facet['Sample'] + ', ' + hist_facet['Condition']

    # plot raw facets
    col_wrap = 5
    sns.set_style('white')

    g_raw = sns.FacetGrid(
        hist_facet, col='for_plot', col_wrap=col_wrap, height=1.27,
        aspect=(1.27 / 1.27), sharex=True, sharey=False
    )

    # use hexbins for plotting
    if self.hexbins:
        g_raw.map(
            plt.hexbin, 'signal', 'Area', gridsize=self.hexbinGridSize,
            linewidths=0.02, color='dimgrey'
        )

    # use scatter points for plotting
    else:
        g_raw.map(plt.scatter, 'signal', 'Area', s=0.05, linewidths=0.0, color='k')

    g_raw.set_titles(
        col_template="{col_name}", fontweight='bold',
        size=np.log(650 / len(g_raw.axes.flatten())), pad=0.0
    )

    for ax in g_raw.axes.flatten():
        ax.tick_params(axis='both', which='major', labelsize=5.0, pad=-2)

        ax.set_xticks([])
        ax.set_yticks([])

        ax.xaxis.label.set_size(np.log(750 / len(g_raw.axes.flatten())))
        ax.yaxis.label.set_size(np.log(750 / len(g_raw.axes.flatten())))

        ax.xaxis.labelpad = 1.0
        ax.yaxis.labelpad = 1.0

        if self.hexbins:
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
        else:
            ax.spines['left'].set_linewidth(0.1)
            ax.spines['bottom'].set_linewidth(0.1)

    # linear interpol. (between 1 and 30 rows at col_wrap = 5 plots/row)
    # for top/bottom subplot adjustment
    pad = (
        1 + (math.ceil(len(g_raw.axes.flatten()) / col_wrap) - 1) * (10 - 1) / (30 - 1)
    )

    plt.tight_layout(pad=pad)
    plt.subplots_adjust(left=0.04, right=0.98, hspace=0.8, wspace=0.3)
    
    # save raw figure dimensions to apply to trimmed figure in "next_channel" function below
    fig_size_inches = g_raw.fig.get_size_inches()  
    
    plt.savefig(os.path.join(plot_dir, f'{channel}_raw.png'), dpi=600, bbox_inches='tight')

    # remove plot_widget and layout attributes from Napari viewer if they exist
    if not initial_callback:
        viewer.window.remove_dock_widget(plot_widget)

        count = plot_layout.count()
        for i in range(count - 1, -1, -1):
            item = plot_layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.setParent(None)
    
    # add raw FacetGrid object to a blank canvas
    raw_canvas = FigureCanvas(g_raw.fig)

    # add navigation tool bar and figure canvas to widget
    plot_layout.addWidget(NavigationToolbar(raw_canvas, plot_widget))
    plot_layout.addWidget(raw_canvas)

    ##########################################################################################
    
    @magicgui(
        layout='vertical', call_button='Apply Cutoffs and Move to Next Marker -->'
    )
    def next_channel(channel, fig_size_inches):

        global lowerCutoff
        global upperCutoff
        
        global arbitrary_selection_toggle
        global marker_index

        global dfTest

        # store cutoffs in QC report
        qc_report['pruneOutliers'][channel] = [float(lowerCutoff), float(upperCutoff)]
        
        # sort and dump updated qc_report to YAML file
        qc_report_sorted = sort_qc_report(
            qc_report, module='pruneOutliers', order=abx_channels
        )
        f = open(report_path, 'w')
        yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)

        napari.utils.notifications.show_info(
            f'Sliders updated to ({lowerCutoff:.3f}, {upperCutoff:.3f})'
        )
        print()
        
        # restore default slider values
        lowerCutoff = 0.0  
        upperCutoff = 100.0

        # raw plot has already been closed in percentile_selector, so the current figure 
        # is the trimmed version.

        # ensure raw and trimmed figures have the same dimensions 
        fig_width, fig_height = fig_size_inches
        plt.gcf().set_size_inches(fig_width, fig_height)
        
        plt.savefig(
            os.path.join(plot_dir, f'{channel}_trimmed.png'), dpi=600, bbox_inches='tight'
        )
        
        # go to next sample
        try:
            if arbitrary_selection_toggle:
                marker_index = arbitrary_marker_idx + 1

            channel = abx_channels[marker_index]

            # update dfTrim to dfTest
            if dfTest is not None:
                dfTest.to_parquet(os.path.join(pruning_dir, 'dfTrim.parquet'))
            
            # read trimmed dataframe
            dfTrim = pd.read_parquet(os.path.join(pruning_dir, 'dfTrim.parquet'))

            # restore dfTest back to None before running next marker 
            dfTest = None
            
            initial_callback = False
            callback(
                self, viewer, channel, dfTrim, data, initial_callback, percentiles_widget,
                percentiles_layout, arbitrary_widget, arbitrary_layout, plot_widget,
                plot_layout, pruning_dir, plot_dir, qc_report, report_path 
            )

            marker_index += 1
            arbitrary_selection_toggle = False
        
        except IndexError:
            logger.info('Thresholding complete!')
            QTimer().singleShot(0, viewer.close)
            print()

    # give next_channel access to channel variable passed to callback
    next_channel.channel.bind(channel)

    # give next_channel access to fig_size_inches variable passed to callback
    next_channel.fig_size_inches.bind(fig_size_inches)
    
    next_channel.native.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Fixed,
    )

    plot_layout.addWidget(next_channel.native)

    ##########################################################################################
    
    @magicgui(
        layout='horizontal', call_button='View Outliers', sample={'label': 'Sample Name'}
    )
    def sample_selector(sample: str, channel):

        return sample, channel
    
    # give sample_selector access to channel variable passed to callback
    sample_selector.channel.bind(channel)

    sample_selector.native.setSizePolicy(
        QtWidgets.QSizePolicy.Maximum,
        QtWidgets.QSizePolicy.Maximum,
    )

    ##########################################################################################
    @magicgui(
        layout='vertical',
        call_button=f'Check {channel} Cutoffs',
        lower_cutoff={'widget_type': 'FloatSlider', 'max': 100.0,
                      'label': 'Lower Percentile'},
        upper_cutoff={'widget_type': 'FloatSlider', 'max': 100.0,
                      'label': 'Upper Percentile'},
    )
    def percentile_selector(lower_cutoff: float = 0.0, upper_cutoff: float = 100.0):
        
        global total_low_idxs
        global total_high_idxs
        
        global dfTest
        
        global lowerCutoff
        global upperCutoff

        lowerCutoff = lower_cutoff
        upperCutoff = upper_cutoff

        if not lowerCutoff > upperCutoff:

            # avoid RuntimeWarning: More than 20 figures have been opened.
            # while keeping event loop running.
            plt.close()
            
            ##################################################################################
            # plot trimmed signal distributions

            # create a copy of dfTrim for testing cutoffs
            dfTest = dfTrim.copy()

            # create lists of indices removed according
            # to lower and upper cutoffs (used for data viz in Napari)
            total_low_idxs = []
            total_high_idxs = []

            # apply current percentile cutoffs to individual samples
            for sample in natsorted(dfTest['Sample'].unique()):

                # drop cells < lower cutoff and > than upper cutoff
                indices_to_drop = []

                sample_channel_data = dfTest[dfTest['Sample'] == sample][channel]

                low_drop_idxs = sample_channel_data.index[
                    sample_channel_data < np.percentile(
                        sample_channel_data, lower_cutoff)]
                indices_to_drop.extend(low_drop_idxs)

                high_drop_idxs = sample_channel_data.index[
                    sample_channel_data > np.percentile(
                        sample_channel_data, upper_cutoff)]
                indices_to_drop.extend(high_drop_idxs)

                dfTest.drop(labels=set(indices_to_drop), axis=0, inplace=True, errors='raise')

                # rescale residual signal intensities
                trimmed_data = dfTest[dfTest['Sample'] == sample][channel]

                scaler = (
                    MinMaxScaler(feature_range=(0, 1), copy=True)
                    .fit(trimmed_data.values.reshape(-1, 1)))
                rescaled_data = scaler.transform(trimmed_data.values.reshape(-1, 1))
                rescaled_data = pd.DataFrame(
                    data=rescaled_data, index=trimmed_data.index
                ).rename(columns={0: channel})

                # update residual signal intensities
                dfTest.update(rescaled_data)

                # update lists of total indices
                total_low_idxs.extend(low_drop_idxs)
                total_high_idxs.extend(high_drop_idxs)

            # melt trimmed and rescaled dfTest
            dfTest_channel = dfTest[['Sample', 'Condition', 'Area'] + [channel]].copy()
            hist_facet = (
                dfTest_channel
                .sample(frac=1.0)
                .melt(id_vars=['Sample', 'Condition', 'Area'],
                      var_name='channel', value_name='signal')
            )

            # naturally sort hist_facet by sample
            hist_facet['Sample'] = pd.Categorical(
                hist_facet['Sample'], ordered=True,
                categories=natsorted(hist_facet['Sample'].unique()))
            hist_facet.sort_values('Sample', inplace=True)
            hist_facet['Sample'] = hist_facet['Sample'].astype('str')

            # create column for facet labels
            hist_facet['for_plot'] = hist_facet['Sample'] + ', ' + hist_facet['Condition']
            
            sns.set_style('white')

            # plot trimmed facets
            g_trimmed = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=col_wrap,
                height=1.27, aspect=(1.27 / 1.27), sharex=True, sharey=False
            )

            # use hexbins for plotting
            if self.hexbins:
                g_trimmed.map(
                    plt.hexbin, 'signal', 'Area', gridsize=self.hexbinGridSize,
                    linewidths=0.02, color='dimgrey'
                )

            # use scatter plots for plotting
            else:
                g_trimmed.map(plt.scatter, 'signal', 'Area', s=0.05, linewidths=0.0, color='k')

            g_trimmed.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=np.log(650 / len(g_trimmed.axes.flatten())), pad=0.0)

            for ax in g_trimmed.axes.flatten():
                ax.tick_params(axis='both', which='major', labelsize=5.0, pad=-2)

                ax.set_xticks([])
                ax.set_yticks([])

                ax.xaxis.label.set_size(np.log(750 / len(g_trimmed.axes.flatten())))
                ax.yaxis.label.set_size(np.log(750 / len(g_trimmed.axes.flatten())))

                ax.xaxis.labelpad = 1.0
                ax.yaxis.labelpad = 1.0

                if self.hexbins:
                    ax.spines['left'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                else:
                    ax.spines['left'].set_linewidth(0.1)
                    ax.spines['bottom'].set_linewidth(0.1)

                # linear interpol. (between 1 and 30 rows at 5 plots/row)
                # for top/bottom subplot adjustment
                pad = (
                    1 + (math.ceil(len(g_trimmed.axes.flatten()) / col_wrap) - 1) *
                    (10 - 1) / (30 - 1)
                )

                plt.tight_layout(pad=pad)
                plt.subplots_adjust(left=0.04, right=0.98, hspace=0.8, wspace=0.3)
            
            # remove old widgets from plot_layout
            count = plot_layout.count()
            # print('widget children:', plot_widget.children())
            
            for i in range(count - 1, -1, -1):
                item = plot_layout.itemAt(i)
                widget = item.widget()
                if widget:
                    widget.setParent(None)

            # add updated plot widgets to plot_layout
            trimmed_canvas = FigureCanvas(g_trimmed.fig)

            # add navigation tool bar and figure canvas to widget
            plot_layout.addWidget(NavigationToolbar(trimmed_canvas, plot_widget))
            plot_layout.addWidget(trimmed_canvas)
            
            # add sample_selector and next_channel button to plot_layout 
            plot_layout.addWidget(sample_selector.native)
            plot_layout.addWidget(next_channel.native)

        else:
            print()
            napari.utils.notifications.show_warning(
                'lowerCutoff must be less than upperCutoff (red).'
            )
            pass
    ##########################################################################################
    
    @sample_selector.called.connect
    def view_points(value: str):

        sample = value[0]
        ch = value[1]
        
        channel_number = marker_channel_number(self, markers, ch)

        if sample in dfTrim['Sample'].unique():

            # remove existing layers
            for i in reversed(range(len(viewer.layers))):
                viewer.layers.pop(i)

            # add DNA1 image to Napari viewer
            file_path = get_filepath(self, check, sample, 'TIF')
            channel_number = marker_channel_number(self, markers, self.counterstainChannel)
            dna, min, max = single_channel_pyramid(file_path, channel=channel_number)
            viewer.add_image(
                dna, rgb=False, opacity=1.0, name=self.counterstainChannel,
                contrast_limits=(min, max)
            )

            # read antibody channel
            file_path = get_filepath(self, check, sample, 'TIF')
            channel_number = marker_channel_number(self, markers, ch)
            channel, min, max = single_channel_pyramid(file_path, channel=channel_number)
            viewer.add_image(
                channel, rgb=False, blending='additive', colormap='green',
                visible=False, name=ch, contrast_limits=(min, max)
            )

            # read segmentation outlines
            file_path = get_filepath(self, check, sample, 'SEG')
            seg, min, max = single_channel_pyramid(file_path, channel=0)
            viewer.add_image(
                seg, rgb=False, blending='additive', opacity=1.0, colormap='gray',
                visible=False, name='segmentation', contrast_limits=(min, max)
            )

            # grab centroids of low signal intensity outliers
            low_centroids = dfTrim[
                ['Y_centroid', 'X_centroid']][
                (dfTrim.index.isin(total_low_idxs)) &
                (dfTrim['Sample'] == sample)]

            # grab centroids of high signal intensity outliers
            high_centroids = dfTrim[
                ['Y_centroid', 'X_centroid']][
                (dfTrim.index.isin(total_high_idxs)) &
                (dfTrim['Sample'] == sample)]

            viewer.add_points(
                low_centroids, name='low centroids', properties=None,
                face_color='magenta', edge_color='k', edge_width=0.0, size=8.0)

            viewer.add_points(
                high_centroids, name='high centroids', properties=None,
                face_color='cyan', edge_color='k', edge_width=0.0, size=8.0)

            napari.utils.notifications.show_info(f'Viewing outliers in sample {sample}')
        else:
            print()
            napari.utils.notifications.show_warning('Sample name not in filtered data.')
            pass

    percentile_selector.native.setSizePolicy(
        QtWidgets.QSizePolicy.Minimum,
        QtWidgets.QSizePolicy.Fixed,
    )

    if initial_callback:  
        percentiles_layout.addWidget(percentile_selector.native)
    
    ##########################################################################################
   
    @magicgui(
        layout='vertical', call_button='Enter', channel={'label': 'Re-start from Marker'}
    )
    def channel_selector(channel: str):

        return channel

    channel_selector.native.setSizePolicy(
        QtWidgets.QSizePolicy.Fixed,
        QtWidgets.QSizePolicy.Fixed
    )
    
    if initial_callback:  
        arbitrary_layout.addWidget(channel_selector.native)

    @channel_selector.called.connect
    def channel_callback(value: str):

        print()

        global arbitrary_selection_toggle
        global arbitrary_marker_idx
        global marker_index
        
        channel = value

        if channel in abx_channels:
        
            dfTrim = data.copy()

            if channel in qc_report['pruneOutliers'].keys():

                # apply subset of cutoffs to dfTrim
                subset_dict = keys_before_key(
                    qc_report['pruneOutliers'], abx_channels, channel
                )
                
                if subset_dict:
                    for ch, (lower_cutoff, upper_cutoff) in subset_dict.items():
                        
                        for sample in natsorted(dfTrim['Sample'].unique()):

                            sample_channel_data = dfTrim[dfTrim['Sample'] == sample][ch]

                            # drop cells < lower cutoff and > than upper cutoff
                            indices_to_drop = []

                            indices_to_drop.extend(
                                sample_channel_data.index[
                                    sample_channel_data < np.percentile(
                                        sample_channel_data,
                                        qc_report['pruneOutliers'][ch][0])])

                            indices_to_drop.extend(
                                sample_channel_data.index[
                                    sample_channel_data > np.percentile(
                                        sample_channel_data,
                                        qc_report['pruneOutliers'][ch][1])])

                            dfTrim.drop(
                                labels=set(indices_to_drop), axis=0, inplace=True,
                                errors='raise'
                            )

                            # rescale pruned antibody signal intensities
                            trimmed_data = dfTrim[dfTrim['Sample'] == sample][ch]

                            scaler = (
                                MinMaxScaler(feature_range=(0, 1), copy=True)
                                .fit(trimmed_data.values.reshape(-1, 1)))
                            rescaled_data = scaler.transform(
                                trimmed_data.values.reshape(-1, 1))
                            rescaled_data = pd.DataFrame(
                                data=rescaled_data, index=trimmed_data.index,
                            ).rename(columns={0: ch})

                            dfTrim.update(rescaled_data)
                else:
                    # select the first marker to pass to the callback function
                    channel = abx_channels[0]
                    
                    # rescale first channel's signal intensities 0-1 per sample
                    for sample in natsorted(dfTrim['Sample'].unique()):
                        raw_channel_data = dfTrim[dfTrim['Sample'] == sample][channel]
                        scaler = (
                            MinMaxScaler(feature_range=(0, 1), copy=True)
                            .fit(raw_channel_data.values.reshape(-1, 1)))
                        rescaled_data = scaler.transform(
                            raw_channel_data.values.reshape(-1, 1)
                        )
                        rescaled_data = pd.DataFrame(
                            data=rescaled_data, index=raw_channel_data.index
                        ).rename(columns={0: channel})
                        dfTrim.update(rescaled_data)
                
                dfTrim.to_parquet(os.path.join(pruning_dir, 'dfTrim.parquet'))
                
                initial_callback = False
                callback(
                    self, viewer, channel, dfTrim, data, initial_callback,
                    percentiles_widget, percentiles_layout, arbitrary_widget,
                    arbitrary_layout, plot_widget, plot_layout, pruning_dir,
                    plot_dir, qc_report, report_path
                )

                arbitrary_selection_toggle = True
                arbitrary_marker_idx = abx_channels.index(channel)
            
            else:
                print()
                napari.utils.notifications.show_info(
                    f'Cutoffs for {channel} have not been defined yet.'
                )

        else:
            napari.utils.notifications.show_warning('Marker name not in filtered data.')

    ##########################################################################################
    
    # delete and re-add percentile_selector to percentiles_layout if it exists 
    # so percentiles_widget appears above plot_widget in Napari and has (0, 100) defaults
    if not initial_callback:
        
        count = percentiles_layout.count()
        for i in range(count - 1, -1, -1):
            item = percentiles_layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.setParent(None)
        
        viewer.window.remove_dock_widget(percentiles_widget)
        percentiles_layout.addWidget(percentile_selector.native)
        viewer.window.add_dock_widget(
            percentiles_widget, name='Select Percentile Cutoffs', area='right'
        )

    # dock (or re-dock) plot_widget to Napari window 
    viewer.window.add_dock_widget(
        plot_widget, name=f'{channel} Intensity vs. Segmentation Area', area='right'
    )

    # remove and re-dock arbitrary_widget if it exists 
    # so plot_widget appears above arbitrary_widget in Napari window
    if not initial_callback:
        viewer.window.remove_dock_widget(arbitrary_widget)
        viewer.window.add_dock_widget(
            arbitrary_widget, name='Re-define Cutoff Series', area='right'
        )


def pruneOutliers(data, self, args):

    global marker_index
    
    print()

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # create pruning directory if it doesn't already exist
    pruning_dir = os.path.join(self.outDir, 'pruning')
    if not os.path.exists(pruning_dir):
        os.mkdir(pruning_dir)

    # create directory to save raw and trimmed intensity distributions if it doesn't already exist
    plot_dir = os.path.join(pruning_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # read QC report
    report_path = os.path.join(self.outDir, 'cylinter_report.yml')
    try:
        qc_report = yaml.safe_load(open(report_path))
        reload_report = False
        if qc_report is None:
            qc_report = {}
            reload_report = True
        if 'pruneOutliers' not in qc_report or qc_report['pruneOutliers'] is None:
            qc_report['pruneOutliers'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(qc_report, module='pruneOutliers', order=None)
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

    # generate Qt widget for specifying percentile cutoffs
    percentiles_widget = QtWidgets.QWidget()
    percentiles_layout = QtWidgets.QVBoxLayout(percentiles_widget)
    percentiles_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Preferred,
        QtWidgets.QSizePolicy.Fixed,
    )
    
    # generate Qt widget for plotting marker signal distributions
    plot_widget = QtWidgets.QWidget()
    plot_layout = QtWidgets.QVBoxLayout(plot_widget)
    plot_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Preferred,
        QtWidgets.QSizePolicy.Fixed
    )
    
    # generate Qt widget for specifying arbitrary marker selections
    arbitrary_widget = QtWidgets.QWidget()
    arbitrary_layout = QtWidgets.QVBoxLayout(arbitrary_widget)
    arbitrary_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Preferred,
        QtWidgets.QSizePolicy.Fixed,
    )

    # make a copy of the dataframe for serial redaction
    dfTrim = data.copy()

    # ensure that all marker keys in QC report have valid values
    valid_keys = []
    for ch in abx_channels:
        try:  
            lowerCutoff, upperCutoff = qc_report['pruneOutliers'][ch]
            if isinstance(lowerCutoff, float) and isinstance(upperCutoff, float):
                valid_keys.append(ch)
        except: 
            pass

    for ch in abx_channels:
        if ch not in valid_keys:
            first_missing_key = ch
            marker_index = abx_channels.index(first_missing_key) + 1
            
            # trim and rescale values of all pruneOutliers channels in QC report
            # whose index is less than that of first missing channel
            for ch in abx_channels:
                if abx_channels.index(ch) < abx_channels.index(first_missing_key):

                    for sample in natsorted(dfTrim['Sample'].unique()):

                        sample_channel_data = dfTrim[dfTrim['Sample'] == sample][ch]

                        # drop cells < lower cutoff and > than upper cutoff
                        indices_to_drop = []

                        indices_to_drop.extend(
                            sample_channel_data.index[
                                sample_channel_data < np.percentile(
                                    sample_channel_data,
                                    qc_report['pruneOutliers'][ch][0])])

                        indices_to_drop.extend(
                            sample_channel_data.index[
                                sample_channel_data > np.percentile(
                                    sample_channel_data,
                                    qc_report['pruneOutliers'][ch][1])])

                        dfTrim.drop(
                            labels=set(indices_to_drop), axis=0, inplace=True,
                            errors='raise'
                        )

                        # rescale pruned antibody signal intensities
                        trimmed_data = dfTrim[dfTrim['Sample'] == sample][ch]

                        scaler = (
                            MinMaxScaler(feature_range=(0, 1), copy=True)
                            .fit(trimmed_data.values.reshape(-1, 1)))
                        rescaled_data = scaler.transform(
                            trimmed_data.values.reshape(-1, 1))
                        rescaled_data = pd.DataFrame(
                            data=rescaled_data, index=trimmed_data.index,
                        ).rename(columns={0: ch})

                        dfTrim.update(rescaled_data)

            # save trimmed and rescaled dataframe
            dfTrim.to_parquet(os.path.join(pruning_dir, 'dfTrim.parquet'))

            viewer.window.add_dock_widget(
                percentiles_widget, name='Select Percentile Cutoffs', area='right'
            )
            
            initial_callback = True
            callback(
                self, viewer, first_missing_key, dfTrim, data, initial_callback,
                percentiles_widget, percentiles_layout, arbitrary_widget,
                arbitrary_layout, plot_widget, plot_layout, pruning_dir,
                plot_dir, qc_report, report_path
            )

            viewer.window.add_dock_widget(
                arbitrary_widget, name='Re-define Cutoff Series', area='right'
            )
            
            viewer.scale_bar.visible = True
            viewer.scale_bar.unit = 'um'

            napari.run()

            break
            # avoids, AttributeError: 'Window' object has no attribute '_qt_window'

    ##############################################################################################
    # prune outliers and rescale 0-1 per channel per sample

    for channel in abx_channels:

        # check that QC report entires are valid (i.e. have not been edited, deleted, etc.)
        try:
            # test for missing channel keys
            qc_report['pruneOutliers'][channel]

            # test for missing thresholds 
            lowerCutoff, upperCutoff = qc_report['pruneOutliers'][channel]
            
            # test for None-type thresholds, occurs when value is deleted in QC report 
            # but value placeholder remains (i.e. hyphans in cylinter_report.yml file)
            if lowerCutoff is None or upperCutoff is None:
                print()
                logger.info(
                    f'Aborting; QC report metadata for channel {channel} is invalid. '
                    'Removing this key from QC report; please re-run pruneOutliers module '
                    'to re-curate thresholds for this channel.'
                )
                del qc_report['pruneOutliers'][channel]
                f = open(report_path, 'w')
                yaml.dump(qc_report, f, sort_keys=False, allow_unicode=False)
                sys.exit()
            
        except (TypeError, KeyError, ValueError):
            print()
            logger.info(
                f'Aborting; QC report metadata for channel {channel} is invalid. '
                'Removing this key from QC report; please re-run pruneOutliers module '
                'to re-curate thresholds for this channel.'
            )
            del qc_report['pruneOutliers'][channel]
            f = open(report_path, 'w')
            yaml.dump(qc_report, f, sort_keys=False, allow_unicode=False)
            sys.exit()

        lowerCutoff, upperCutoff = qc_report['pruneOutliers'][channel]
        
        if (lowerCutoff == 0.0) and (upperCutoff == 100.0):
            logger.info(f'All data points selected for {channel} channel')
        else:
            logger.info(
                f'Applying percentile cutoffs ({lowerCutoff:.3f}, '
                f'{upperCutoff:.3f}) to {channel} channel.'
            )

        for sample in natsorted(data['Sample'].unique()):

            sample_channel_data = data[data['Sample'] == sample][channel]

            # drop cells < lower cutoff and > than upper cutoff
            indices_to_drop = []

            indices_to_drop.extend(
                sample_channel_data.index[
                    sample_channel_data < np.percentile(
                        sample_channel_data, lowerCutoff)])

            indices_to_drop.extend(
                sample_channel_data.index[
                    sample_channel_data > np.percentile(
                        sample_channel_data, upperCutoff)])

            data.drop(
                labels=set(indices_to_drop), axis=0, inplace=True, errors='raise'
            )

            # rescale trimmed antibody signal intensities
            trimmed_data = data[data['Sample'] == sample][channel]

            scaler = (
                MinMaxScaler(feature_range=(0, 1), copy=True)
                .fit(trimmed_data.values.reshape(-1, 1)))
            rescaled_data = scaler.transform(
                trimmed_data.values.reshape(-1, 1))
            rescaled_data = pd.DataFrame(
                data=rescaled_data, index=trimmed_data.index,
            ).rename(columns={0: channel})

            data.update(rescaled_data)

    ##############################################################################################
    # prune outliers and rescale 0-1 per channel (scale is across samples in this case)
    
    # for k, v in cutoffs_dict.items():
    #     print(f'Applying percentile cutoffs to the {k} channel.')

    #     marker_channel_data = data[k]

    #     scaler = (
    #         MinMaxScaler(feature_range=(0, 1), copy=True)
    #         .fit(marker_channel_data.values.reshape(-1, 1)))
    #     rescaled_data = scaler.transform(
    #         marker_channel_data.values.reshape(-1, 1))
    #     rescaled_data = pd.DataFrame(
    #         data=rescaled_data, index=marker_channel_data.index,
    #     ).rename(columns={0: k})

    #     data.update(rescaled_data)

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
