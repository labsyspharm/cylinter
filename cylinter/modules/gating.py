import os
import sys
import yaml
import math
import logging

import numpy as np
import pandas as pd
from itertools import combinations
from natsort import natsorted
import collections

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.lines import Line2D

from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas

from PyPDF2 import PdfReader, PdfWriter

from matplotlib.backends.qt_compat import QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
from matplotlib.widgets import Slider, Button
from qtpy.QtCore import QTimer

import napari
from magicgui import magicgui

from ..utils import (
    input_check, read_markers, single_channel_pyramid, marker_channel_number,
    log_banner, log_multiline, get_filepath, reorganize_dfcolumns, sort_qc_report
)

from ..config import BooleanTerm

logger = logging.getLogger(__name__)


def invert_bool(value):
    if value is None:
        return value
    return not value


def generate_pdf(data, marker, abx_channels, qc_report, gate_dir, dist_dir):

    logger.info(f'Writing PDF page for {marker}')
    
    ncols = 7
    nrows = math.ceil(len(data['Sample'].unique()) / ncols)
    shift = 1600
    pad = 50
    canvas_width = ncols * shift
    canvas_height = nrows * shift + pad

    # generate pdf canvas
    my_canvas = canvas.Canvas(
        os.path.join(dist_dir, f'{marker}.pdf'),
        pagesize=(canvas_width, canvas_height)
    )

    # specify x, y origin coords for first histogram on the page
    x_start, y_start = (pad + shift, canvas_height - shift)

    # add PDF figure label
    my_canvas.setFont('Helvetica', 250)
    my_canvas.drawString(pad, canvas_height - (shift / 2), marker)

    reportlab_graphics = {}
    for sample in natsorted(data['Sample'].unique()):

        cond = data['Condition'][data['Sample'] == sample].unique().item()
        rep = data['Replicate'][data['Sample'] == sample].unique().item()
        title = f'{sample}_{cond}_{rep}'

        hist_input = data[marker][(data['Sample'] == sample)]
        area_input = data['Area'][(data['Sample'] == sample)]

        #######################################################################
        # percentile filter for viz (no data filtration by default)
        percentile = np.percentile(hist_input, 0)
        hist_input = hist_input[hist_input >= percentile]
        area_input = area_input[area_input.index.isin(hist_input.index)]
        #######################################################################

        sns.set_style('whitegrid')
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.scatter(x=hist_input, y=area_input, s=2000 / len(hist_input), label='marker dist.')
        plt.xlim([0, 1])
        plt.ylim([area_input.min(), area_input.max()])

        # generate gate annotation and add to plot
        gate_val = qc_report['gating'][f'{marker}, {sample}, {cond}, {rep}']
        
        if gate_val is not None:
            ax.axvline(x=gate_val, color='orange', linewidth=4.0, linestyle='-')

        plt.title(title, weight='bold')
        plt.ylabel('Segmentation Area')

        plt.tight_layout()
        fig.savefig(os.path.join(dist_dir, 'temp.png'), dpi=400)
        plt.close('all')

        drawing = ImageReader(os.path.join(dist_dir, 'temp.png'))

        # delete temp.svg
        os.remove(os.path.join(dist_dir, 'temp.png'))

        # append to reportlab_graphics dict
        reportlab_graphics[title] = drawing

    # generate an ordered dict to ensure lexicographic order
    # and format pdf pages with marker histograms
    od = collections.OrderedDict(natsorted(reportlab_graphics.items()))
    ticker = 0

    # loop over reportlab_graphics dict
    for k, v in od.items():

        # advance tickers by 1
        ticker += 1

        my_canvas.drawImage(v, x_start, y_start)

        # wrap plots on canvas (every ncols-1 for first row)
        if ticker < ncols - 1:
            x_start += shift  # advance carriage to the right

        else:
            if ticker == ncols - 1:
                x_start = pad
                y_start -= shift  # carriage return
                ticker = ncols * 2  # wrap every ncols from here on

            elif ticker % ncols == 0:
                x_start = pad
                y_start -= shift  # carriage return

            else:
                x_start += shift  # advance carriage to the right

    # generate new pdf page
    my_canvas.showPage()

    # save pdf pages to disk
    my_canvas.save()


def multipage_pdf(abx_channels, dist_dir):
    
    logger.info('Writing multi-page PDF.')
    
    if os.path.exists(os.path.join(dist_dir, 'multipage.pdf')):
        os.remove(os.path.join(dist_dir, 'multipage.pdf'))

    file_list = [i for i in os.listdir(dist_dir) if not i.startswith('.')]
    file_order = [f'{i}.pdf' for i in abx_channels]
    sorted_files = sorted(
        file_list, key=lambda x: file_order.index(x) if
        x in file_order else len(file_order)
    )

    output_pdf = PdfWriter()
    for pdf_file in sorted_files:

        pdf = PdfReader(os.path.join(dist_dir, pdf_file))

        for page in pdf.pages:
            output_pdf.add_page(page)

    with open(os.path.join(dist_dir, 'multipage.pdf'), 'wb') as file:
        output_pdf.write(file)


def callback(self, viewer, data, hist_widget, hist_layout, selection_widget, selection_layout, gate_dir, sample, marker, initial_callback, dist_dir, abx_channels, markers, qc_report, report_path):

    # if valid sample and marker entered
    if (sample in data['Sample'].unique()) and (marker in abx_channels):

        check, markers_filepath = input_check(self)

        # clear existing channels from Napari window if they exist
        viewer.layers.clear()

        cond = data['Condition'][data['Sample'] == sample].unique().item()
        rep = data['Replicate'][data['Sample'] == sample].unique().item()

        sample_data = data[['X_centroid', 'Y_centroid', marker, 'Area']][data['Sample'] == sample]

        ###################################################################
        # percentile filter for viz (no data filtration by default)
        percentile = np.percentile(sample_data[marker], 0)
        sample_data = sample_data[sample_data[marker] >= percentile]
        ###################################################################

        # read DNA1 channel
        file_path = get_filepath(self, check, sample, 'TIF')
        channel_number = marker_channel_number(self, markers, self.counterstainChannel)
        dna, min, max = single_channel_pyramid(file_path, channel=channel_number)
        viewer.add_image(
            dna, rgb=False, blending='additive', colormap='gray', visible=False,
            name='DNA', contrast_limits=(min, max)
        )
            
        # read segmentation outlines
        file_path = get_filepath(self, check, sample, 'SEG')
        seg, min, max = single_channel_pyramid(file_path, channel=0)
        viewer.add_image(
            seg, rgb=False, blending='additive', opacity=1.0, colormap='gray',
            visible=False, name='segmentation', contrast_limits=(min, max)
        )

        # read marker channel
        reversed_abx_channels = abx_channels.copy()
        reversed_abx_channels.reverse()
        reversed_abx_channels = [marker]  # overriding to single channel
        
        for ch in reversed_abx_channels:

            channel_number = marker_channel_number(self, markers, ch)
            file_path = get_filepath(self, check, sample, 'TIF')
            img, min, max = single_channel_pyramid(file_path, channel=channel_number)

            if ch == marker:
                visible = True
            else:
                visible = False

            viewer.add_image(
                img, rgb=False, blending='additive', colormap='green', 
                visible=visible, name=ch, contrast_limits=(min, max)
            )

            ###################################################################
            # apply previously defined contrast limits if they exist

            try:
                viewer.layers[ch].contrast_limits = (
                    qc_report['setContrast'][ch][0], qc_report['setContrast'][ch][1]
                ) 
            except KeyError:
                pass

            ###################################################################

        # add centroids of currently gated cells
        current_gate = qc_report['gating'][f'{marker}, {sample}, {cond}, {rep}']
        
        current_gateOG = current_gate  # recording None status

        if current_gate is None:
            pass
        else:
            centroids = sample_data[['Y_centroid', 'X_centroid']][
                sample_data[marker] > current_gate]

            viewer.add_points(
                centroids, name='Reference gate', face_color='#00aaff', 
                edge_color='#00aaff', edge_width=0.0, size=7.0, opacity=1.0,
                blending='opaque', visible=True
            )

        try:
            # remove hist_widget and layout elements from Napari viewer if it exists
            viewer.window.remove_dock_widget(hist_widget)
        except LookupError:
            pass

        # remove hist_widget and layout attributes from Napari viewer if they exist
        count = hist_layout.count()
        for i in range(count - 1, -1, -1):
            item = hist_layout.itemAt(i)
            widget = item.widget()
            if widget:
                widget.setParent(None)

        # generate a blank figure canvas
        gate_canvas = FigureCanvas(Figure(figsize=(5.5, 4)))

        # add navigation tool bar and figure canvas to widget
        hist_layout.addWidget(NavigationToolbar(gate_canvas, hist_widget))
        hist_layout.addWidget(gate_canvas)

        # set plot style
        sns.set_style('white')

        # get figure object from canvas
        fig = gate_canvas.figure

        # adjust plot on canvas to taste
        fig.subplots_adjust(left=0.25, bottom=0.25)

        # set plot title
        fig.suptitle(
            f'Marker={marker}, Sample={sample}, Condition={cond}, Replicate={rep}', size=10
        )

        # get axis object from canvas
        ax = gate_canvas.figure.subplots()
        
        # plot scatter on canvas axis
        ax.scatter(
            x=sample_data[marker], y=sample_data['Area'], color='steelblue', 
            ec='white', lw=0.1, alpha=1.0, s=80000 / len(sample_data[marker]),
            label='--- reference gate'
        )

        # set y-axis label
        ax.set_ylabel('Segmentation Area')

        xlims = plt.xlim([sample_data[marker].min(), sample_data[marker].max()])
        plt.ylim([sample_data['Area'].min(), sample_data['Area'].max()])

        # add vertical line at current gate
        if current_gate is not None:
            ax.axvline(x=current_gate, color='k', linewidth=1.0, linestyle='--')

            ax.text(
                0.87, 1.04, '--- reference gate', fontsize=9, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes
            )

        # add sliders to plot
        axcolor = 'lightgoldenrodyellow'
        axGate = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

        rnge = [xlims[0], xlims[1]]  
        # small value added to avoid slider bar errors (line 350) if image has no signal

        # add slider functionality
        if current_gate is None:
            valinit = sample_data[marker].min()
            ax.axvline(x=valinit, c='orange', linewidth=2.5)
        else:
            valinit = current_gate

        sGate = Slider(
            axGate, 'Gate', rnge[0], rnge[1], valinit=valinit, valstep=(rnge[1] / 100000)
        )
        sGate.label.set_fontsize(11)
        sGate.label.set_color('k')
        sGate.poly.set_facecolor('orange')

        # specify function for updating sliders
        def update(val):

            # remove current gate
            [i.remove() for i in ax.get_lines()]

            # new gate
            adjusted_gate = np.float64(sGate.val).astype('float16')
            
            # update plot with adjusted gate
            ax.axvline(x=adjusted_gate, c='orange', linewidth=2.5)

            if current_gateOG is not None:

                # maintain current gate
                ax.axvline(x=current_gate, color='k', linewidth=1.0, linestyle='--')
            
            return adjusted_gate

        # update slider when moved
        sGate.on_changed(update)

        # add button to show selected centroids in Napari viewer
        button_ax = fig.add_axes([0.65, 0.025, 0.25, 0.06])
        button = Button(button_ax, 'Plot Points', color=axcolor, hovercolor='0.975')
        button.label.set_fontsize(11)

        def apply_cutoff(event):

            # get current gate
            adjusted_gate = update(val=None)

            # apply gate
            sample_update = sample_data[
                sample_data[marker] > adjusted_gate].copy()
            
            # isolate x, y coordinates of newly gated centroids
            centroids = sample_update[['Y_centroid', 'X_centroid']]

            # remove existing centroids and
            # plot new centroid selection in Napari window
            if current_gateOG is not None:
                layers = 5
            else:
                layers = 4

            # if not centroids.empty:
            if len(viewer.layers) == layers:
                viewer.layers.pop(layers - 1)

            viewer.add_points(
                centroids, name='adjusted gate', face_color='orange', edge_color='orange',
                edge_width=0.0, size=7.0, opacity=1.0, blending='opaque'
            )

        button.on_clicked(apply_cutoff)

        # maintain reference to button after exiting sample_selector_callback()
        button_ax._button = button
       
        #######################################################################
        viewer.window.add_dock_widget(
            hist_widget, name=f'{marker} intensity distribution for sample {sample}', area='right'
        )
        #######################################################################
        
        if not initial_callback:
            try:
                # remove selection_widget dock from layout if it exists
                viewer.window.remove_dock_widget(selection_widget)
            except LookupError:
                pass

            # re-dock selection_widget so hist_widget appears first in Napari
            viewer.window.add_dock_widget(
                selection_widget, name='Arbitrary Sample/Marker Selection', area='right'
            )

        #######################################################################

        @magicgui(layout='horizontal', call_button='Apply Gate and Move to Next Sample -->')
        def next_sample():
                
            print()
            
            adjusted_gate = update(val=None)

            # save gate value
            qc_report['gating'][f'{marker}, {sample}, {cond}, {rep}'] = float(adjusted_gate)

            # dump updated qc_report to YAML file
            qc_report_sorted = sort_qc_report(
                qc_report, module='gating', order=list(qc_report['gating'].keys())
            )
            f = open(report_path, 'w')
            yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
            
            napari.utils.notifications.show_info(f'Gate updated to {adjusted_gate:.4f}.')
            
            if any(value is None for value in qc_report['gating'].values()):
                for key, value in qc_report['gating'].items():
                    if qc_report['gating'][key] is None:
                        marker_name = key.split(', ')[0]
                        sample_name = key.split(', ')[1]
                        break

                print()
                initial_callback = False
                callback(
                    self, viewer, data, hist_widget,
                    hist_layout, selection_widget, selection_layout,
                    gate_dir, sample_name, marker_name, initial_callback,
                    dist_dir, abx_channels, markers, qc_report, report_path
                )
            else:
                print()
                
                logger.info('Gating complete!')

                ##################################################################################
                # ensure all PDFs are updated with current gates
                
                print()
                logger.info('Updating gates in PDF pages.')
                
                for ch in abx_channels:
                    generate_pdf(data, ch, abx_channels, qc_report, gate_dir, dist_dir)
                multipage_pdf(abx_channels, dist_dir)

                ##################################################################################
                
                print()
                
                QTimer().singleShot(0, viewer.close)

        next_sample.native.setSizePolicy(
            QtWidgets.QSizePolicy.Maximum,
            QtWidgets.QSizePolicy.Maximum,
        )
 
        hist_layout.addWidget(next_sample.native)

        #######################################################################
        # initialize "Enter" button
    
        @magicgui(
            layout='vertical',
            call_button='Enter',
            marker_name={'label': 'Marker Name'},
            sample_name={'label': 'Sample Name'}
        )
        def sample_selector(marker_name: str, sample_name: str):

            return marker_name, sample_name

        sample_selector.native.setSizePolicy(
            QtWidgets.QSizePolicy.Fixed,
            QtWidgets.QSizePolicy.Fixed
        )
        
        if initial_callback:   
            selection_layout.addWidget(sample_selector.native)

        #######################################################################

        @sample_selector.called.connect
        def sample_selector_callback(value: str):
            
            marker_name = value[0]
            sample_name = value[1]

            initial_callback = False
            callback(
                self, viewer, data, hist_widget,
                hist_layout, selection_widget, selection_layout,
                gate_dir, sample_name, marker_name, initial_callback,
                dist_dir, abx_channels, markers, qc_report, report_path
            )

        #######################################################################

        @magicgui(
            layout='vertical',
            call_button='Refresh PDF(s)',
            marker_name={'label': 'Marker Name (or "ALL")'},
        )
        def update_pdf(marker_name: str, gate_dir, dist_dir):

            return marker_name

        # give update_pdf access to zeros, gate_dir, and dist_dir
        # update_pdf.zeros.bind(zeros)
        update_pdf.gate_dir.bind(gate_dir)
        update_pdf.dist_dir.bind(dist_dir)

        update_pdf.native.setSizePolicy(
            QtWidgets.QSizePolicy.Maximum,
            QtWidgets.QSizePolicy.Maximum,
        )
        
        if initial_callback:
            selection_layout.addWidget(update_pdf.native)

        #######################################################################

        @update_pdf.called.connect
        def update_pdf_callback(marker: str):

            if marker not in abx_channels + ['ALL']:
                napari.utils.notifications.show_warning('Marker name not in filtered data.')
                pass
            else:
                print()
                if marker == 'ALL':
                    for marker in abx_channels:
                        generate_pdf(data, marker, abx_channels, qc_report, gate_dir, dist_dir)
                    multipage_pdf(abx_channels, dist_dir)
                else:
                    generate_pdf(data, marker, abx_channels, qc_report, gate_dir, dist_dir)
                    multipage_pdf(abx_channels, dist_dir)

                napari.utils.notifications.show_info('PDF(s) updated!')
        
        #######################################################################

        napari.utils.notifications.show_info(f'Gating {marker} in sample {sample}')
    
    else:
        print()
        if marker not in abx_channels and sample not in data['Sample'].unique():
            napari.utils.notifications.show_warning(
                'Marker and sample names not in filtered data.'
            )
            pass
        elif marker not in abx_channels and sample in data['Sample'].unique():
            napari.utils.notifications.show_warning('Marker name not in filtered data.')
            pass
        elif sample not in data['Sample'].unique() and marker in abx_channels:
            napari.utils.notifications.show_warning('Sample name not in filtered data.')
            pass


# main
def gating(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, abx_channels = read_markers( 
        markers_filepath=markers_filepath,
        counterstain_channel=self.counterstainChannel,
        markers_to_exclude=self.markersToExclude, data=None
    )
    
    if self.gating:
        
        print()

        # sox2 = pd.read_csv(
        #     '/Volumes/T7 Shield/cylinter_input/sandro_tma/output/gating/'
        #     'sox2_cells.csv')
        # sox2['Sample'] = sox2['Sample'].astype('str')
        # data = data.merge(sox2, how='inner', on=['CellID', 'Sample'])

        # create gating and distributions subdirectories if they haven't already
        gate_dir = os.path.join(self.outDir, 'gating')
        if not os.path.exists(gate_dir):
            os.makedirs(gate_dir)

        dist_dir = os.path.join(gate_dir, 'distributions')
        if not os.path.exists(dist_dir):
            os.makedirs(dist_dir)

        # read QC report
        report_path = os.path.join(self.outDir, 'cylinter_report.yml')
        try:
            qc_report = yaml.safe_load(open(report_path))
            reload_report = False
            if qc_report is None:
                qc_report = {}
                reload_report = True
            if 'gating' not in qc_report or qc_report['gating'] is None:
                qc_report['gating'] = {}
                reload_report = True
            if reload_report:
                qc_report_sorted = sort_qc_report(qc_report, module='gating', order=None)
                f = open(report_path, 'w')
                yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
                qc_report = yaml.safe_load(open(report_path))
        except:
            logger.info(
                'Aborting; QC report missing from CyLinter output directory. Re-start pipeline '
                'from aggregateData module to initialize QC report.'
            )
            sys.exit()
    
        data = data[~data['Sample'].isin(self.samplesToRemoveGating)]

        ##########################################################################################
        # initialize gating dictionary

        gating_keys = [
            f"{j}, {i}, {data['Condition'][data['Sample'] == i].unique().item()}, "
            f"{data['Replicate'][data['Sample'] == i].unique().item()}" 
            for j in abx_channels for i in natsorted(data['Sample'].unique())
        ]
        
        # pad gating dictionary with all possible gate keys
        # (ensures keys accidentally deleted from QC report are accounted for)
        for key in gating_keys:
            if key not in qc_report['gating']:
                qc_report['gating'][key] = None

        # sort and dump updated qc_report to YAML file
        qc_report_sorted = sort_qc_report(qc_report, module='gating', order=gating_keys)
        f = open(report_path, 'w')
        yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
        qc_report = yaml.safe_load(open(report_path))  # reload qc_report

        ##########################################################################################
        # generate initial PDFs of marker distributions
        
        for marker in abx_channels:
            if not os.path.exists(os.path.join(dist_dir, f'{marker}.pdf')):
                generate_pdf(data, marker, abx_channels, qc_report, gate_dir, dist_dir)
            
        if not os.path.exists(os.path.join(dist_dir, 'multipage.pdf')):
            logger.info('Generated multi-channel PDF.')
            multipage_pdf(abx_channels, dist_dir)
            print()
        
        ##########################################################################################
        # gate data

        viewer = napari.Viewer(title='CyLinter')
        
        # generate sample selection Qt widget
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

        if any(value is None for value in qc_report['gating'].values()):
            for key, value in qc_report['gating'].items():
                if qc_report['gating'][key] is None:
                    marker_name = key.split(', ')[0]
                    sample_name = key.split(', ')[1]
                    break

            initial_callback = True
            callback(
                self, viewer, data, hist_widget,
                hist_layout, selection_widget, selection_layout,
                gate_dir, sample_name, marker_name, initial_callback,
                dist_dir, abx_channels, markers, qc_report, report_path
            )
            
            viewer.window.add_dock_widget(
                selection_widget, name='Arbitrary Sample/Marker Selection', area='right'
            )

            viewer.scale_bar.visible = True
            viewer.scale_bar.unit = 'um'

            napari.run()

        else:
            pass

        ##########################################################################################
        # subtract gates and binarize data
        
        gated = pd.DataFrame()

        for sample in natsorted(data['Sample'].unique()):

            logger.info(f'Applying gates to sample {sample}')

            cond = data['Condition'][data['Sample'] == sample].unique().item()
            rep = data['Replicate'][data['Sample'] == sample].unique().item()
            
            # initialize dataframe to store zeroed sample data
            gated_temp = pd.DataFrame()

            gated_temp[['Sample', 'CellID']] = data[['Sample', 'CellID']][
                data['Sample'] == sample
            ]

            for marker in abx_channels:

                sample_data = data[marker][data['Sample'] == sample]

                gate = qc_report['gating'][f'{marker}, {sample}, {cond}, {rep}']

                if gate is None:
                    print()
                    logger.info(
                        'Aborting; some gates in QC report are missing (i.e. Null). '
                        'Re-run the gating module to ensure all sample/marker '
                        'combinations have a gate.'
                    )
                    sys.exit()

                gated_temp[marker] = sample_data - gate

            gated = pd.concat([gated, gated_temp], ignore_index=True)
        
        # include gate subtracted signal intensities in the output dataframe
        # gated[[f'{i}_gated' for i in abx_channels]] = gated.iloc[:, 2:]
        
        gated[abx_channels] = gated[abx_channels] > 0.0

        data = data.merge(
            gated, how='inner', on=['Sample', 'CellID'], suffixes=(None, '_bool')
        )

        bool_cols = [f'{i}_bool' for i in abx_channels]

        data['vector'] = data[bool_cols].apply(
            lambda row: ''.join('1' if cell else '0' for cell in row), axis=1
        )

        total_vector_counts = (
            data
            .groupby('vector')
            .size()
            .sort_values(ascending=False)
        )

        selected_vector_counts = total_vector_counts[total_vector_counts >= self.vectorThreshold]

        print()
        logger.info(
            '%d Boolean vectors with >= %d events.',
            len(selected_vector_counts), self.vectorThreshold
        )

        ##########################################################################################
        # plot Boolean vector counts

        sns.set_style('white')

        fig, ax = plt.subplots()

        plt.bar(
            x=list(range(len(total_vector_counts))),
            height=total_vector_counts, lw=0.0, color='grey'
        )
        ax.set_xlabel('Vector', fontsize=15, labelpad=10)
        ax.set_ylabel('Count', fontsize=15, labelpad=10)
        plt.yscale('log')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(os.path.join(gate_dir, 'total_vector_counts.pdf'))
        plt.close('all')

        ##########################################################################################
        # plot heatmap of Boolean vectors >= vectorThreshold
        font_scaler = 1.0
        
        unique_vectors = (
            data[['vector'] + bool_cols][
                data['vector'].isin(selected_vector_counts.index)]
            .drop_duplicates(subset=bool_cols)
            .set_index(keys='vector')
            .reindex(selected_vector_counts.index)
        )

        if not unique_vectors.empty:
            
            num_rows = unique_vectors.shape[0]
            num_cols = unique_vectors.shape[1]

            sns.set(font_scale=font_scaler)
            sns.set_style('white')
            fig, axs = plt.subplots(
                1, 2, figsize=(num_cols, (num_rows / 2) + 2), sharey=False,
                gridspec_kw=dict(width_ratios=[1, 0.25])
            )
            
            heatmap = sns.heatmap(
                unique_vectors, annot=False, lw=0.1, linecolor='k', xticklabels=True,
                cmap='Greys', vmax=1.75, cbar=False, square=False, ax=axs[0]
            )

            heatmap.set_xticklabels(
                [i.get_text().split('_bool')[0] for
                 i in heatmap.get_xticklabels()], rotation=90
            )

            for _, spine in heatmap.spines.items():
                spine.set_linewidth(2)
            axs[0].set_ylabel('')
            axs[0].set_yticks([])
            axs[0].set_xlabel('Marker', fontsize=18, labelpad=20, fontweight='bold')

            selected_vector_counts = pd.DataFrame(selected_vector_counts)
            selected_vector_counts.reset_index(inplace=True)
            selected_vector_counts.rename(columns={0: 'counts'}, inplace=True)
            
            # avoiding matplotlib info statement: "Using categorical units to plot a list
            # of strings that are all parsable as floats or dates." caused by the
            # below plotting function
            mlogger = logging.getLogger('matplotlib')
            mlogger.setLevel(logging.WARNING)

            sns.barplot(
                x=selected_vector_counts['counts'],
                y=selected_vector_counts['vector'],
                orient='horizontal', color='grey', ax=axs[1]
            )
            plt.xscale('log')
            axs[1].invert_yaxis()
            axs[1].spines['left'].set_visible(False)
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['top'].set_visible(False)
            axs[1].set_xlabel('Cell Count', fontsize=18, labelpad=10, fontweight='bold')
            axs[1].set_ylabel('')
            axs[1].set_yticks([])

            plt.subplots_adjust(wspace=-0.1)
            plt.tight_layout()
            plt.savefig(os.path.join(gate_dir, 'threshold_vectors.pdf'))
            plt.close('all')

        else:
            logger.info('No vector signatures to plot for threshold_vectors.pdf, skipping.')

        ##########################################################################################
        # classify cells in dataframe
        
        # add ignored markers to inner dictionary of self.classes
        for clss, inner_dict in self.classes.items():
            terms = []
            for name, channels in inner_dict.items():
                terms.extend(channels) 
            channel_list = [
                str(i).split('~')[1] if '~' in str(i) else str(i) for i in terms
            ]
            ignore = sorted(set(abx_channels).difference(set(channel_list)))
            boo = [BooleanTerm.parse_str(t) for t in ignore]
            inner_dict['ignore'] = boo
            self.classes[clss] = inner_dict

        # create expanded dictionary of immunomarker signatures
        signatures = {}
        for clss, inner_dict in self.classes.items():

            # generate all combinations of subset markers
            combos = []
            for length in range(0, len(inner_dict['subsets']) + 1):
                combos.extend(combinations(inner_dict['subsets'], length))
            combos = [sorted(i) for i in combos]

            # name signatures
            signature_names = [f"{'_'.join(i)}_{clss}" for i in combos]
            signature_names[0] = clss  # removing "_" from base phenotype
            
            # generate complementary lists of negated subset markers
            negated = [
                list(set(inner_dict['subsets']).difference(set(i))) for i in combos
            ]
            negated = [sorted(i) for i in negated]
            
            # prepend minus signs to negated subset markers 
            negated = [[f'-{element}' for element in inner_list] for inner_list in negated]

            # prepend plus signs to expressed subset markers 
            combos = [[f'+{element}' for element in inner_list] for inner_list in combos]
            
            # combine expressed and negated subset markers
            subsets = [
                expressed + not_expressed for expressed, not_expressed in zip(combos, negated)
            ]
            
            # write cell type subsets to signatures dictionary
            for name, i in zip(signature_names, subsets):
                boo = [BooleanTerm.parse_str(t) for t in i]
                boo = inner_dict['definition'] + boo
                boo = boo + inner_dict['ignore']
                signatures[name] = boo

        # write signatures to yaml file 
        f = open(os.path.join(gate_dir, 'signatures.yml'), 'w')
        sigs = {
            key: [f'{value}/~{value}' if value.negated is None else
                  str(value) for value in values] for key, values in signatures.items()
        }
        yaml.dump(sigs, f)

        if signatures: 
            data['class'] = None
            for class_name, terms in signatures.items():
                
                positives = pd.Series(
                    np.ones(len(data), dtype=bool), index=data.index
                )  # row indices to pair down in IDing cells of each class
                
                for term in terms:
                    if term.negated is not None:
                        
                        try:
                            col = data[f'{term.name}_bool']
                        except KeyError:
                            logger.warning(
                                f'Aborting; classes dictionary term {term.name} '
                                'not a marker in dataframe.'
                            )
                            sys.exit()
                        
                        if term.negated:
                            col = ~col
                        positives = positives[col]

                indexer = (positives.index, 'class')
                conflicts = set(data.loc[indexer][data.loc[indexer].notna()])
                if conflicts:
                    raise ValueError(
                        f"Boolean class '{class_name}' overlaps with {conflicts}"
                    )
                data.loc[indexer] = class_name

            data['class'] = data['class'].fillna('unclassified')
            data['class'] = data['class'].astype('str')

        else:
            print()
            logger.info(
                'No cell state classifications have been made. Please '
                'update "classes" dictionary in config.yml.'
            )
            sys.exit()

        classified_counts = data.groupby('class').size()
        classified_counts.sort_values(ascending=False, inplace=True)
        classified_counts.name = 'cell_count'
        classified_counts.index.name = 'class'
        classified_counts = classified_counts.reset_index()

        log_banner(logger.info, 'Boolean classifications')
        log_multiline(logger.info, classified_counts.to_string(index=False))

        pct_classified = classified_counts['cell_count'][
            classified_counts['class'] != 'unclassified'].sum() / len(data) * 100

        unclassified = data[data['class'] == 'unclassified']

        unclassified_vector_counts = (
            unclassified
            .groupby('vector')
            .size()
            .sort_values(ascending=False)
        )

        unclassified_vector_counts = unclassified_vector_counts[
            unclassified_vector_counts >= self.vectorThreshold
        ]

        unique_unclassified_vectors = (
            unclassified[['vector'] + bool_cols][
                unclassified['vector'].isin(unclassified_vector_counts.index)]
            .drop_duplicates(subset=bool_cols)
            .set_index(keys='vector')
        )
        
        logger.info('')
        logger.info(
            'Current classification accounts for %.2f%% of data; '
            '%d Boolean vectors with >= %d events remain unclassified.',
            pct_classified, len(unique_unclassified_vectors), self.vectorThreshold
        )

        ##########################################################################################
        # plot heatmap of Boolean vectors >= vectorThreshold left unclassified

        if not unique_unclassified_vectors.empty:

            sns.set(font_scale=1)
            sns.set_style('white')

            num_rows = unique_unclassified_vectors.shape[0]
            num_cols = unique_unclassified_vectors.shape[1]

            fig, axs = plt.subplots(
                1, 2, figsize=(num_cols, (num_rows / 2) + 2), sharey=False,
                gridspec_kw=dict(width_ratios=[1, 0.25])
            )

            heatmap = sns.heatmap(
                unique_unclassified_vectors, annot=False, lw=0.1, linecolor='k',
                xticklabels=True, cmap='Greys', vmax=1.75, cbar=False,
                square=False, ax=axs[0]
            )

            heatmap.set_xticklabels(
                [i.get_text().split('_bool')[0] for i in heatmap.get_xticklabels()], rotation=90
            )

            for _, spine in heatmap.spines.items():
                spine.set_linewidth(2)

            axs[0].set_ylabel('')
            axs[0].set_yticks([])
            axs[0].set_xlabel('Marker', fontsize=18, labelpad=20, fontweight='bold')

            unclassified_vector_counts = pd.DataFrame(unclassified_vector_counts)
            unclassified_vector_counts.rename(columns={0: 'counts'}, inplace=True)
            
            sns.barplot(
                data=unclassified_vector_counts, x='counts',
                y=unclassified_vector_counts.index, orient='horizontal',
                color='grey', ax=axs[1]
            )

            plt.xscale('log')
            axs[1].invert_yaxis()
            axs[1].spines['left'].set_visible(False)
            axs[1].spines['right'].set_visible(False)
            axs[1].spines['top'].set_visible(False)
            axs[1].set_xlabel('Cell Count', fontsize=18, labelpad=10, fontweight='bold')
            axs[1].set_ylabel('')
            axs[1].set_yticks([])
            plt.subplots_adjust(wspace=-0.1)
            plt.tight_layout()
            plt.savefig(os.path.join(gate_dir, 'unclassified_vectors.pdf'))
            plt.close('all')

        ##########################################################################################
        # plot heatmap of classified cell type signatures
        
        classes = pd.DataFrame.from_dict(signatures)
        table = pd.DataFrame(columns=abx_channels)

        for cls in classes.columns:
            class_declarations = classes[cls].copy()
            class_declarations.index = [i.name for i in class_declarations]
            channel_negations = [class_declarations[i].negated for i in abx_channels]
            row = dict(zip(abx_channels, channel_negations))
            table = pd.concat(
                [table, pd.DataFrame(index=[cls], data=[row])], ignore_index=False
            )

        # apply the custom function to invert Boolean calls in the dataFrame
        table = table.map(invert_bool)
        
        # create two Boolean heatmaps; one in which "don't cares" are filled
        # with 0s, and one in which they are filled with 1s
        black = table.copy()
        black[black.isna()] = 0
        black = black.astype('int')
        
        white = table.copy()
        white[white.isna()] = 3
        white = white.astype('int')
        white[white == 0] = 2
        white[white == 3] = 1

        num_classes = table.shape[0]
        num_markers = table.shape[1]
        
        fig, ax = plt.subplots(figsize=(num_markers / 2, num_classes / 2))

        x = np.arange(num_markers + 1)
        y = np.arange(num_classes + 1)
        xs, ys = np.meshgrid(x, y[::-1])

        triangles1 = [
            (i + j * (num_markers + 1), i + 1 + j * (num_markers + 1),
             i + (j + 1) * (num_markers + 1)) for j in range(num_classes) for
            i in range(num_markers)
        ]
        triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1)

        triangles2 = [
            (i + 1 + j * (num_markers + 1), i + 1 + (j + 1) * (num_markers + 1),
             i + (j + 1) * (num_markers + 1)) for j in range(num_classes) for
            i in range(num_markers)
        ]
        triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2)
        
        ax.tripcolor(triang1, white.values.ravel(), lw=0.0, cmap='Greys_r', vmin=0, vmax=2) 
        ax.tripcolor(triang2, 1 - black.values.ravel(), lw=0.0, cmap='Greys_r', vmin=-1, vmax=1)
        # vmin and vmax values set to create grey color

        for i in range(num_markers + 1):
            ax.axvline(x=i, ymin=0, ymax=num_classes, c='k', lw=0.7)
        plt.xlim([0, num_markers])

        for i in range(num_classes + 1):
            ax.axhline(y=i, xmin=0, xmax=num_markers, c='k', lw=0.7)
        plt.ylim([0, num_classes])

        custom_xtick_locations = list(range(num_markers))
        custom_xtick_labels = abx_channels
        plt.xticks(
            [i + 0.5 for i in custom_xtick_locations], custom_xtick_labels,
            fontsize=10, rotation=90
        )

        custom_ytick_locations = list(range(num_classes))
        custom_ytick_labels = table.index
        plt.yticks(
            [i + 0.5 for i in custom_ytick_locations[::-1]], custom_ytick_labels,
            fontsize=10, rotation=0
        )

        legend_elements = []

        legend_elements.append(
            Line2D([0], [0], marker='s', color='none', label='True',
                   markerfacecolor='gray', markeredgecolor='k', lw=0.01, markersize=12)
        )

        legend_elements.append(
            Line2D([0], [0], marker='s', color='none', label='False',
                   markerfacecolor='white', markeredgecolor='k', lw=0.01, markersize=12)
        )

        ax.legend(
            handles=legend_elements, prop={'size': 8}, loc='upper left',
            bbox_to_anchor=[1.01, 1.1], labelspacing=1.0, frameon=False
        )

        plt.title(
            'Classification accounts for {:.2f}% of data'.format(pct_classified), size=8
        )
        plt.savefig(os.path.join(gate_dir, 'class_signatures.pdf'), bbox_inches='tight')
        plt.close('all')

    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
