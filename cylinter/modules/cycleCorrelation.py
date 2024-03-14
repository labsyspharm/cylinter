import os
import re
import sys
import yaml
import logging

import numpy as np
import pandas as pd

from natsort import natsorted

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

import napari
from magicgui import magicgui

from ..utils import (
    input_check, read_markers, marker_channel_number, sort_qc_report,
    single_channel_pyramid, categorical_cmap, get_filepath, reorganize_dfcolumns, compute_gmm
)

logger = logging.getLogger(__name__)

arbitrary_selection_toggle = False
sample_index = 1


def callback(self, viewer, sample, samples_to_run, data, ratios_melt, initial_callback, selection_widget, selection_layout, hist_widget, hist_layout, cycles_dir, qc_report, report_path): 
    
    if sample in data['Sample'].unique():
        
        print()

        check, markers_filepath = input_check(self)

        # read marker metadata
        markers, abx_channels = read_markers( 
            markers_filepath=markers_filepath,
            counterstain_channel=self.counterstainChannel,
            markers_to_exclude=self.markersToExclude, data=None
        )

        # clear existing channels from Napari window if they exist
        viewer.layers.clear()

        # isolate sample data
        group = ratios_melt[ratios_melt['sample'] == sample]
        
        # loop over last cycle only
        cycle_ratio = group['cycle'].unique()[-1]

        # isolate ratio data
        cycle_data = group[group['cycle'] == cycle_ratio]
        first_dna = cycle_ratio.split('/')[0]
        last_dna = cycle_ratio.split('/')[1]

        # add cell segmentation outlines to Napari viewer
        file_path = get_filepath(self, check, sample, 'SEG')
        seg, min, max = single_channel_pyramid(file_path, channel=0)
        viewer.add_image(
            seg, rgb=False, blending='additive',
            opacity=0.5, colormap='gray', visible=False,
            name='segmentation', contrast_limits=(min, max)
        )
 
        # add last DNA image to Napari viewer
        file_path = get_filepath(self, check, sample, 'TIF')
        channel_number = marker_channel_number(self, markers, last_dna)
        dna_last, min, max = single_channel_pyramid(file_path, channel=channel_number)
        viewer.add_image(
            dna_last, rgb=False, blending='additive', colormap='magenta', name=last_dna,
            contrast_limits=(min, max)
        )

        # add first DNA image to Napari viewer
        file_path = get_filepath(self, check, sample, 'TIF')
        channel_number = marker_channel_number(self, markers, first_dna)
        dna_first, min, max = single_channel_pyramid(file_path, channel=channel_number)
        viewer.add_image(
            dna_first, rgb=False, blending='additive', colormap='green', name=first_dna,
            contrast_limits=(min, max)
        )

        # remove hist_widget dock and layout attributes from layout if they exist
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

        # set plot title
        log10 = '$Log_{10}$'
        fig.suptitle(f'Sample={sample}  {log10}({self.counterstainChannel}/{last_dna})', size=10)

        # get axis object from canvas
        ax = canvas.figure.subplots()
        ax2 = ax.twinx()

        # plot log(cycle 1/n) histogram for current sample
        counts, bins, patches = ax.hist(
            cycle_data['log10(ratio)'], bins=self.numBinsCorrelation,
            density=False, color='grey', ec='none', alpha=0.75,
            histtype='stepfilled', range=None, label='before'
        )

        ax.set_ylabel('Count', size=13)
        ax2.set_ylabel('GMM density')
        
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

        # add sliders to plot
        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = fig.add_axes(
            [0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = fig.add_axes(
            [0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        # specify data range
        rnge = [bins.min(), bins.max()]

        # convert histogram data into 2D numpy array with 1 column to pass to GMM
        gmm_data = cycle_data['log10(ratio)'].values.reshape(-1, 1)

        # compute GMM
        comp_lower_percentile, comp_upper_percentile = compute_gmm(
            data=gmm_data, x_min=rnge[0], x_max=rnge[1], ax=ax2
        )
        try:
            lowerCutoff, upperCutoff = qc_report['cycleCorrelation'][sample]
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

        # specify function for updating sliders
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
            lowerCutoff, upperCutoff = update(val=None)

            # get indices in log(ratio) series outside
            # lower and upper cutoffs
            idxs = list(
                cycle_data['index'][
                    (cycle_data['log10(ratio)'] < lowerCutoff) |
                    (cycle_data['log10(ratio)'] > upperCutoff)])

            # filter group data by selecting
            # indices NOT in idxs
            sample_data = data[data['Sample'] == sample]
            drop_df = sample_data.index.isin(idxs)
            centroids = sample_data[
                ['Y_centroid', 'X_centroid']][~drop_df]

            # remove existing centroids and plot new centroid selection in Napari window
            if not centroids.empty:
                if len(viewer.layers) == 4:
                    viewer.layers.pop(3)
                viewer.add_points(
                    centroids,
                    name='Selected Cells',
                    properties=None,
                    face_color='yellow',
                    edge_color='k',
                    edge_width=0.0, size=7.0)

        # add button functionality
        button.on_clicked(apply_cutoffs)

        # maintain reference to button after exiting callback()
        button_ax._button = button
        
        # dock (or re-dock) hist_widget to Napari window 
        viewer.window.add_dock_widget(
            hist_widget, name=f'Log({first_dna}/{last_dna}) histogram', area='right'
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
            call_button='Apply Gates and Move to Next Sample -->'
        )
        def next_sample(sample):

            global arbitrary_selection_toggle
            global sample_index
            
            # get current cutoffs
            lowerCutoff, upperCutoff = update(val=None)

            if lowerCutoff <= upperCutoff:
            
                # store cutoffs in QC report
                qc_report['cycleCorrelation'][sample] = [float(lowerCutoff), float(upperCutoff)]

                # sort and dump updated qc_report to YAML file
                qc_report_sorted = sort_qc_report(
                    qc_report, module='cycleCorrelation', order=None
                )
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
                        self, viewer, sample, samples_to_run, data, ratios_melt, initial_callback,
                        selection_widget, selection_layout, hist_widget, hist_layout,
                        cycles_dir, qc_report, report_path 
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
                self, viewer, sample, samples_to_run, data, ratios_melt, initial_callback,
                selection_widget, selection_layout, hist_widget, hist_layout,
                cycles_dir, qc_report, report_path
            )

            arbitrary_selection_toggle = True
        
        ###########################################################################
        napari.utils.notifications.show_info(f'Viewing sample {sample}')
        
    else:
        print()
        napari.utils.notifications.show_warning(
            'Sample name not in filtered data.'
        )


# main
def cycleCorrelation(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )

    # create cycles directory if it doesn't already exist
    cycles_dir = os.path.join(self.outDir, 'cycles')
    if not os.path.exists(cycles_dir):
        os.makedirs(cycles_dir)
    
    # read QC report
    report_path = os.path.join(self.outDir, 'cylinter_report.yml')
    try:
        qc_report = yaml.safe_load(open(report_path))
        reload_report = False
        if qc_report is None:
            qc_report = {}
            reload_report = True
        if 'cycleCorrelation' not in qc_report or qc_report['cycleCorrelation'] is None:
            qc_report['cycleCorrelation'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(qc_report, module='cycleCorrelation', order=None)
            f = open(report_path, 'w')
            yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
            qc_report = yaml.safe_load(open(report_path))
    except:
        logger.info(
            'Aborting; QC report missing from CyLinter output directory. Re-start pipeline '
            'from aggregateData module to initialize QC report.'
        )
        sys.exit()
    
    ##########################################################
    # get basename of DNA channels so they will not be read by Napari as immunomarker channels  
    dna_moniker = str(re.search(r'[^\W\d]+', self.counterstainChannel).group())  
    
    # get ordered list of DNA cycles
    dna_cycles = natsorted(data.columns[data.columns.str.contains(dna_moniker)])

    # compute DNA ratios
    ratios = pd.DataFrame(
        [np.log10((data[self.counterstainChannel] + 0.001) / (data[i] + 0.001))
         for i in dna_cycles]).T

    # computing ratios changes columns headers, creating new ratio column headers
    unnamed_headers = [
        i for i in ratios.columns if i.startswith('Unnamed')]
    ratio_headers = [f'{self.counterstainChannel}/{i}' for i in [j for j in dna_cycles[1:]]]
    ratio_columns = dict(zip(unnamed_headers, ratio_headers))
    ratio_columns[self.counterstainChannel] = (
        f'{self.counterstainChannel}/{self.counterstainChannel}'
    )
    ratios.rename(columns=ratio_columns, inplace=True)
    ratios['sample'] = data['Sample']

    # melt ratios dataframe
    ratios_melt = (
        ratios
        .reset_index()
        .melt(id_vars=['sample', 'index'], var_name='cycle',
              value_name='log10(ratio)'))

    # convert sample and cycle columns to ordered categoricals
    # and sort naturally on sample, cycle, and index
    ratios_melt['sample'] = pd.Categorical(
        ratios_melt['sample'], ordered=True, categories=natsorted(ratios_melt['sample'].unique())
    )
    ratios_melt['cycle'] = pd.Categorical(
        ratios_melt['cycle'], ordered=True, categories=natsorted(ratios_melt['cycle'].unique())
    )
    ratios_melt = ratios_melt.sort_values(['sample', 'cycle', 'index'])

    # convert columns back to strings
    ratios_melt['sample'] = ratios_melt['sample'].astype('str')
    ratios_melt['cycle'] = ratios_melt['cycle'].astype('str')

    ##########################################################
    
    # initialize Napari viewer
    viewer = napari.Viewer(title='CyLinter')

    # generate arbitrary sample selection Qt widget
    selection_widget = QtWidgets.QWidget()
    selection_layout = QtWidgets.QVBoxLayout(selection_widget)
    selection_widget.setSizePolicy(
        QtWidgets.QSizePolicy.Minimum,
        QtWidgets.QSizePolicy.Fixed,
    )
    
    # generate distribution Qt widget
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
                lowerCutoff, upperCutoff = qc_report['cycleCorrelation'][sample]
                if not isinstance(lowerCutoff, float) or not isinstance(upperCutoff, float):
                    samples_to_run.append(sample)
                
            except: 
                samples_to_run.append(sample)
        
    try:
        sample = samples_to_run[0]
    
        initial_callback = True
        callback(
            self, viewer, sample, samples_to_run, data, ratios_melt, initial_callback,
            selection_widget, selection_layout, hist_widget, hist_layout,
            cycles_dir, qc_report, report_path
        )
        
        viewer.window.add_dock_widget(
            selection_widget, name='Arbitrary Sample Selection', area='right'
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

    plot_dir = os.path.join(cycles_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    
    idxs_to_drop = {}
    for sample in samples_to_run:

        # isolate sample data
        group = ratios_melt[ratios_melt['sample'] == sample]
        
        # loop over last cycle only
        cycle_ratio = group['cycle'].unique()[-1]
        last_dna = cycle_ratio.split('/')[1]
        
        # isolate ratio data
        cycle_data = group[group['cycle'] == cycle_ratio]

        lowerCutoff, upperCutoff = qc_report['cycleCorrelation'][sample]

        # select all data points if sliders were not adjusted
        if lowerCutoff == upperCutoff:
            logger.info(f'All data points selected for sample {sample}.')
            # select all data points if sliders were not adjusted
            lowerCutoff = cycle_data['log10(ratio)'].min()
            upperCutoff = cycle_data['log10(ratio)'].max()
        else:
            logger.info(
                f'Applying cutoffs ({lowerCutoff:.3f}, '
                f'{upperCutoff:.3f}) to sample {sample}'
            )

        # plot DNA ratio histogram BEFORE filtering
        fig, ax = plt.subplots()

        # plot DNA ratio histogram BEFORE filtering
        counts, bins, patches = plt.hist(
            cycle_data['log10(ratio)'], bins=self.numBinsCorrelation,
            density=False, color='b', ec='none', alpha=0.5,
            histtype='stepfilled', range=None, label='before'
        )

        # apply lower and upper cutoffs
        idxs = list(
            cycle_data['index'][
                (cycle_data['log10(ratio)'] < lowerCutoff) |
                (cycle_data['log10(ratio)'] > upperCutoff)]
        )

        cycle_data_filtered = cycle_data[~cycle_data['index'].isin(idxs)]

        # plot DNA ratio histogram AFTER filtering
        counts, bins, patches = plt.hist(
            cycle_data_filtered['log10(ratio)'], bins=bins,
            density=False, color='r', ec='none', alpha=0.5,
            histtype='stepfilled', range=None, label='before'
        )
        
        log10 = '$Log_{10}$'
        plt.xlabel(f'{log10}({self.counterstainChannel}/{last_dna}')
        plt.ylabel('Count')
        plt.title(f'Sample = {sample})', size=10)

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

        idxs_to_drop[sample] = idxs

    # filter cells from all samples
    indices_to_drop = set()
    for k, v in idxs_to_drop.items():
        for idx in v:
            indices_to_drop.update(set([idx]))
    data = data.loc[~data.index.isin(indices_to_drop)]

    print()
    logger.info('Plotting cycle correlation graphs')
    
    # grab dna and sample columns of filtered dataframe
    facet_input = data.loc[:, data.columns.str.contains(f'{dna_moniker}|Sample')].copy()

    # melt filtered dataframe
    facet_per_cycle_melt = (
        facet_input
        .sample(frac=1.0)
        .reset_index()
        .melt(id_vars=['Sample', 'index'], var_name='cycle'))

    # convert sample and cycle columns to ordered categorical
    # with sorted catgories by natsorted
    facet_per_cycle_melt['Sample'] = pd.Categorical(
        facet_per_cycle_melt['Sample'], ordered=True,
        categories=natsorted(
            facet_per_cycle_melt['Sample'].unique()))
    facet_per_cycle_melt['cycle'] = pd.Categorical(
        facet_per_cycle_melt['cycle'], ordered=True,
        categories=natsorted(
            facet_per_cycle_melt['cycle'].unique()))

    # sort melt on sample, cycle, and index
    facet_per_cycle_melt = facet_per_cycle_melt.sort_values(
        ['Sample', 'cycle', 'index'])

    # plot dna intensity correlation per cycle
    fig, ax = plt.subplots(figsize=(5, 5))

    g = sns.FacetGrid(
        facet_per_cycle_melt, col='cycle', col_wrap=5,
        sharex=True, sharey=False
    )

    g.map(
        lambda y, color: plt.scatter(
            facet_per_cycle_melt['value'].loc[
                facet_per_cycle_melt['cycle'] == self.counterstainChannel], y, s=0.15,
            linewidth=0.0, marker='o', c='r'), 'value'
    )
    g.set(xlabel=None)
    g.set(ylabel=None)
    
    plt.savefig(os.path.join(plot_dir, 'correlation.png'), dpi=600, bbox_inches='tight')
    plt.close('all')

    # plot dna intensity correlation per cycle (color by sample)
    fig, ax = plt.subplots(figsize=(5, 5))

    # build cmap
    cmap = categorical_cmap(
        numUniqueSamples=len(facet_per_cycle_melt['Sample'].unique()),
        numCatagories=10, cmap='tab10', continuous=False
    )

    sample_color_dict = dict(zip(
        natsorted(facet_per_cycle_melt['Sample'].unique()), cmap.colors))

    g = sns.FacetGrid(
        facet_per_cycle_melt, col='cycle', hue='Sample',
        col_wrap=5, sharex=True, sharey=False
    )

    g.map(
        lambda sam, y, color, **kwargs: plt.scatter(
            facet_per_cycle_melt.loc[
                (facet_per_cycle_melt['Sample'] == sam.unique()[0])
                & (facet_per_cycle_melt['cycle'] == self.counterstainChannel), 'value'], y,
            c=np.reshape(sample_color_dict[sam.unique()[0]], (-1, 3)),
            s=0.15, linewidth=0.0, marker='o', **kwargs), 'Sample', 'value'
    )
    g.set(xlabel=None)
    g.set(ylabel=None)
    
    plt.legend(markerscale=10, bbox_to_anchor=(1.1, 1.05))

    plt.savefig(
        os.path.join(plot_dir, 'correlation_per_sample.png'), dpi=600, bbox_inches='tight'
    )
    plt.close('all')

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
