import os
import sys
import yaml
import math
import logging
import pickle

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
    input_check, read_markers, single_channel_pyramid, marker_channel_number, napari_notification,
    log_banner, log_multiline, get_filepath, reorganize_dfcolumns
)

from ..config import BooleanTerm

logger = logging.getLogger(__name__)


def invert_bool(value):
    if value is None:
        return value
    return not value


def generate_pdf(data, marker, abx_channels, zeros, gate_dir, dist_dir):

    napari_notification(f'Writing PDF page for {marker}.')
    
    ncols = 7
    nrows = math.ceil(len(data['Sample'].unique()) / ncols)
    shift = 1600
    pad = 50
    canvas_width = ncols * shift
    canvas_height = nrows * shift + pad

    zeros = pd.read_csv(os.path.join(gate_dir, 'zeros.csv'))
    zeros['sample'] = zeros['sample'].astype(str)

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
        hist_input = hist_input[hist_input > percentile]
        area_input = area_input[area_input.index.isin(hist_input.index)]
        #######################################################################

        sns.set_style('whitegrid')
        fig, ax = plt.subplots(figsize=(4, 4))

        ax.scatter(x=hist_input, y=area_input, s=2000 / len(hist_input), label='marker dist.')
        plt.xlim([0, 1])
        plt.ylim([area_input.min(), area_input.max()])

        # generate gate annoatation and add to plot
        gate_val = zeros['gate'][
            (zeros['marker'] == marker) & (zeros['sample'] == sample)
        ].iloc[0]

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

    zeros.to_csv(os.path.join(gate_dir, 'zeros.csv'), index=False)


def multipage_pdf(abx_channels, dist_dir):
    
    napari_notification('Writing multi-page PDF.')
    
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


def callback(self, viewer, data, zeros, hist_widget, hist_layout, selection_widget,selection_layout, gate_dir, sample, marker, initial_callback, dist_dir, abx_channels, markers):

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
        sample_data = sample_data[sample_data[marker] > percentile]
        ###################################################################

        reversed_abx_channels = abx_channels.copy()
        reversed_abx_channels.reverse()
        reversed_abx_channels = [marker]  # overriding to single channel
        
        for ch in reversed_abx_channels:

            channel_number = marker_channel_number(markers, ch)
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

            if os.path.exists(os.path.join(f'{self.outDir}/contrast/contrast_limits.yml')):

                contrast_limits = yaml.safe_load(open(
                    f'{self.outDir}/contrast/contrast_limits.yml')
                )

                viewer.layers[ch].contrast_limits = contrast_limits[ch][0], contrast_limits[ch][1]

            ###################################################################

        # read DNA1 channel
        file_path = get_filepath(self, check, sample, 'TIF')
        dna, min, max = single_channel_pyramid(file_path, channel=0)
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

        # add centroids of currently gated cells
        current_gate = zeros['gate'][
            (zeros['marker'] == marker) & (zeros['sample'] == sample)
        ].iloc[0]
        
        current_gateOG = current_gate  # recording NaN status

        if math.isnan(current_gate):
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
        gate_canvas = FigureCanvas(Figure())

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
            f'Sample={sample}, Condition={cond}, Replicate={rep}, Marker={marker}', size=10
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

        plt.xlim([sample_data[marker].min(), sample_data[marker].max()])
        plt.ylim([sample_data['Area'].min(), sample_data['Area'].max()])

        # add vertical line at current gate
        if not math.isnan(current_gate):
            ax.axvline(x=current_gate, color='k', linewidth=1.0, linestyle='--')

            ax.text(
                0.87, 1.04, '--- reference gate', fontsize=9, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes
            )

        # add sliders to plot
        axcolor = 'lightgoldenrodyellow'
        axGate = fig.add_axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

        rnge = [sample_data[marker].min(), sample_data[marker].max()]

        # add slider functionality
        if math.isnan(current_gate):
            valinit = sample_data[marker].min()
        else:
            valinit = current_gate
        
        sGate = Slider(
            axGate, 'Adjusted gate', rnge[0], rnge[1], valinit=valinit, valstep=(rnge[1] / 100000)
        )
        
        sGate.label.set_fontsize(11)
        sGate.label.set_color('k')
        sGate.poly.set_facecolor('orange')

        # avoid rounding errors in equivalence between
        # current and adjusted gates
        current_gate = sGate.val.astype('float16')

        # specify function for updating sliders
        def update(val):

            # remove current gate
            [i.remove() for i in ax.get_lines()]

            # new gate
            adjusted_gate = sGate.val.astype('float16')

            # update plot with adjusted gate
            ax.axvline(x=adjusted_gate, c='orange', linewidth=2.5)

            if not math.isnan(current_gateOG):

                # maintain current gate
                ax.axvline(x=current_gate, color='k', linewidth=1.0, linestyle='--')

            napari_notification(f'Gate updated to {adjusted_gate:.4f}.')
            
            return adjusted_gate

        # update slider when moved
        sGate.on_changed(update)

        # add button to show selected centroids in Napari viewer
        button_ax = fig.add_axes([0.65, 0.025, 0.25, 0.06])
        button = Button(button_ax, 'Plot Points', color=axcolor, hovercolor='0.975')
        button.label.set_fontsize(11)

        def apply_cutoff(event):

            # get current gate
            gate = update(val=None)

            # apply gate
            sample_update = sample_data[
                sample_data[marker] > gate].copy()

            # isolate x, y coordinates of newly gated centroids
            centroids = sample_update[['Y_centroid', 'X_centroid']]

            # remove existing centroids and
            # plot new centroid selection in Napari window
            if not math.isnan(current_gateOG):
                layers = 5
            else:
                layers = 4

            if not centroids.empty:
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

            # get current gate
            adjusted_gate = update(val=None)

            if adjusted_gate != current_gate:
            
                zeros = pd.read_csv(os.path.join(gate_dir, 'zeros.csv'))
                zeros['sample'] = zeros['sample'].astype(str)

                row_indexer = zeros.index[
                    (zeros['marker'] == marker) & (zeros['sample'] == sample)
                ].item()
                zeros.loc[row_indexer, 'gate'] = adjusted_gate
                zeros.to_csv(os.path.join(gate_dir, 'zeros.csv'), index=False)

            zeros = pd.read_csv(os.path.join(gate_dir, 'zeros.csv'))
            zeros['sample'] = zeros['sample'].astype(str)
            
            # check gate has been selected for current sample
            row = zeros.index[
                (zeros['marker'] == marker) & (zeros['sample'] == sample)
            ].item()
            
            if math.isnan(zeros.loc[row]['gate']):
                napari_notification(
                    "Select a gate then click 'Apply Gate and Move to Next Sample -->'."
                )
            else:
                if zeros['gate'].isnull().any():
                    for row in range(len(zeros)):
                        if math.isnan(zeros.loc[row]['gate']):
                            marker_name = zeros.loc[row]['marker']
                            sample_name = zeros.loc[row]['sample']
                            break

                    print()
                    initial_callback = False
                    callback(
                        self, viewer, data, zeros, hist_widget,
                        hist_layout, selection_widget, selection_layout,
                        gate_dir, sample_name, marker_name, initial_callback,
                        dist_dir, abx_channels, markers
                    )
                else:
                    print()
                    napari_notification('Gating complete!')
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
            sample_name={'label': 'Sample Name'},
            marker_name={'label': 'Marker Name'},
        )
        def sample_selector(sample_name: str, marker_name: str):

            return sample_name, marker_name

        sample_selector.native.setSizePolicy(
            QtWidgets.QSizePolicy.Minimum,
            QtWidgets.QSizePolicy.Maximum
        )
        
        if initial_callback:   
            selection_layout.addWidget(sample_selector.native)

        #######################################################################

        @sample_selector.called.connect
        def sample_selector_callback(value: str):

            sample_name = value[0]
            marker_name = value[1]

            zeros = pd.read_csv(os.path.join(gate_dir, 'zeros.csv'))
            zeros['sample'] = zeros['sample'].astype(str)

            print()
            initial_callback = False
            callback(
                self, viewer, data, zeros, hist_widget,
                hist_layout, selection_widget, selection_layout,
                gate_dir, sample_name, marker_name, initial_callback,
                dist_dir, abx_channels, markers
            )

        #######################################################################

        @magicgui(
            layout='vertical',
            call_button='Refresh PDF(s)',
            marker_name={'label': 'Marker Name (or ALL)'},
        )
        def update_pdf(marker_name: str, zeros, gate_dir, dist_dir):

            return marker_name

        # give update_pdf access to zeros, gate_dir, and dist_dir
        update_pdf.zeros.bind(zeros)
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

            print()
            if marker == 'ALL':
                for marker in abx_channels:
                    generate_pdf(data, marker, abx_channels, zeros, gate_dir, dist_dir)
                multipage_pdf(abx_channels, dist_dir)
            else:
                generate_pdf(data, marker, abx_channels, zeros, gate_dir, dist_dir)
                multipage_pdf(abx_channels, dist_dir)

            napari_notification('PDF(s) updated!')
            print()
        
        #######################################################################

        napari_notification(f'Gating {marker} in sample {sample}')
    
    else:
        print()
        napari_notification('Invalid entry.')
        pass


# main
def gating(data, self, args):

    check, markers_filepath = input_check(self)

    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
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

        data = data[~data['Sample'].isin(self.samplesToRemoveGating)]

        ##########################################################################################
        # initialize zeros tables

        if not os.path.exists(os.path.join(gate_dir, 'zeros.csv')):

            mylist = [
                f"{j},{i},{data['Condition'][data['Sample'] == i].unique().item()},"
                f"{data['Replicate'][data['Sample'] == i].unique().item()},"
                for j in abx_channels for i in natsorted(data['Sample'].unique())
            ]
            rows = [row.split(',') for row in mylist]
            zeros = pd.DataFrame(
                rows, columns=['marker', 'sample', 'condition', 'replicate', 'gate']
            )
            zeros['sample'] = zeros['sample'].astype(str)
            zeros['gate'] = np.nan
            zeros.to_csv(os.path.join(gate_dir, 'zeros.csv'), index=False)
        
        else:
            zeros = pd.read_csv(os.path.join(gate_dir, 'zeros.csv'))
            zeros['sample'] = zeros['sample'].astype(str)

            mylist = [
                f"{j},{i},{data['Condition'][data['Sample'] == i].unique().item()},"
                f"{data['Replicate'][data['Sample'] == i].unique().item()},"
                for j in abx_channels for i in natsorted(data['Sample'].unique())
            ]
            rows = [row.split(',') for row in mylist]

            zeros_new = pd.DataFrame(
                rows, columns=['marker', 'sample', 'condition', 'replicate', 'gate']
            )
            zeros_new['replicate'] = zeros_new['replicate'].astype('int')
            zeros_new['gate'] = np.nan
            
            zeros = zeros.merge(
                zeros_new, how='outer', suffixes=(None, '_y'),
                on=['marker', 'sample', 'condition', 'replicate']
            )
            zeros.drop(columns=['gate_y'], inplace=True)
            
            zeros.to_csv(os.path.join(gate_dir, 'zeros.csv'), index=False)

        ##########################################################################################
        # generate initial PDFs of marker distributions
        for marker in abx_channels:
            if not os.path.exists(os.path.join(dist_dir, f'{marker}.pdf')):
                generate_pdf(data, marker, abx_channels, zeros, gate_dir, dist_dir)

        if not os.path.exists(os.path.join(dist_dir, 'multipage.pdf')):
            multipage_pdf(abx_channels, dist_dir)

            napari_notification('PDF(s) generated!')
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

        if zeros['gate'].isnull().any():
            for row in range(len(zeros)):
                if math.isnan(zeros.loc[row]['gate']):
                    marker_name = zeros.loc[row]['marker']
                    sample_name = zeros.loc[row]['sample']
                    break
            
            initial_callback = True
            callback(
                self, viewer, data, zeros, hist_widget,
                hist_layout, selection_widget, selection_layout,
                gate_dir, sample_name, marker_name, initial_callback,
                dist_dir, abx_channels, markers
            )
            
            viewer.window.add_dock_widget(
                selection_widget, name='Arbitrary Sample/Marker Selection', area='right'
            )

            viewer.scale_bar.visible = True
            viewer.scale_bar.unit = 'um'

            napari.run()

            print()

        else:
            napari_notification('Gating complete!')
            print()
        
        ##########################################################################################
        # ensure all PDFs are updated with current gates
        
        logger.info('Updating PDFs...')
        
        for marker in abx_channels:
            # if not os.path.exists(os.path.join(dist_dir, f'{marker}.pdf')):
            generate_pdf(data, marker, abx_channels, zeros, gate_dir, dist_dir)
        
        # if not os.path.exists(os.path.join(dist_dir, 'multipage.pdf')):
        multipage_pdf(abx_channels, dist_dir)

        ##########################################################################################
        # subtract gates and binarize data
        
        print()
        
        zeros = pd.read_csv(os.path.join(gate_dir, 'zeros.csv'))
        zeros['sample'] = zeros['sample'].astype(str)
        
        gated = pd.DataFrame()

        for sample in natsorted(data['Sample'].unique()):

            logger.info(f'Applying gates to sample {sample}.')

            # initialize dataframe to store zeroed sample data
            gated_temp = pd.DataFrame()

            gated_temp[['Sample', 'CellID']] = data[['Sample', 'CellID']][
                data['Sample'] == sample
            ]

            for marker in abx_channels:

                sample_data = data[marker][data['Sample'] == sample]

                gate = zeros['gate'][
                    (zeros['marker'] == marker) & (zeros['sample'] == sample)
                ].iloc[0]

                if math.isnan(gate):
                    print()
                    logger.info(
                        'Aborting; zeros.csv contains NaNs. '
                        'Ensure all sample/marker combinations have a gate.'
                    )
                    sys.exit()

                gated_temp[marker] = sample_data - gate

            gated = pd.concat([gated, gated_temp], ignore_index=True)
        
        # include gate subtracted signal intensities in the output dataframe
        # gated[[f'{i}_gated' for i in abx_channels]] = gated.iloc[:, 2:]
        
        gated.loc[:, abx_channels] = gated[abx_channels] > 0

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

        sns.barplot(
            data=selected_vector_counts, x=selected_vector_counts, y=selected_vector_counts.index,
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

        f = open(os.path.join(gate_dir, 'signatures.pkl'), 'wb')
        pickle.dump(signatures, f)
        f.close()

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

            sns.barplot(
                data=unclassified_vector_counts, x=unclassified_vector_counts,
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
            channel_negations = [i.negated for i in classes[cls]]
            row = dict(zip(abx_channels, channel_negations))
            table = pd.concat(
                [table, pd.DataFrame(index=[cls], data=[row])], ignore_index=False
            )

        # apply the custom function to invert Boolean calls in the dataFrame
        table = table.applymap(invert_bool)

        # create two Boolean heatmaps; one in which "don't cares" are filled
        # with 1s, and one in which they are filled with 0s
        black = table.fillna(value=1)
        black = black.astype('int')
        white = table.fillna(value=0)
        white = white.astype('int')

        num_classes = table.shape[1]
        num_markers = table.shape[0]
        fig, ax = plt.subplots(figsize=(num_classes / 2, num_markers / 2))

        x = np.arange(num_classes + 1)
        y = np.arange(num_markers + 1)
        xs, ys = np.meshgrid(x, y[::-1])

        triangles1 = [
            (i + j * (num_classes + 1), i + 1 + j * (num_classes + 1),
             i + (j + 1) * (num_classes + 1)) for j in range(num_markers) for
            i in range(num_classes)
        ]
        triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1)

        triangles2 = [
            (i + 1 + j * (num_classes + 1), i + 1 + (j + 1) * (num_classes + 1),
             i + (j + 1) * (num_classes + 1)) for j in range(num_markers) for
            i in range(num_classes)
        ]
        triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2)

        ax.tripcolor(triang1, white.values.ravel(), lw=0.0, cmap='Greys', vmax=1.75)
        ax.tripcolor(triang2, black.values.ravel(), lw=0.0, cmap='Greys', vmax=1.75)

        for i in range(num_classes + 1):
            ax.axvline(x=i, ymin=0, ymax=num_markers, c='k', lw=0.4)
        plt.xlim([0, num_classes])

        for i in range(num_markers + 1):
            ax.axhline(y=i, xmin=0, xmax=num_classes, c='k', lw=0.4)
        plt.ylim([0, num_markers])

        custom_xtick_locations = list(range(num_classes))
        custom_xtick_labels = abx_channels
        plt.xticks(
            [i + 0.5 for i in custom_xtick_locations], custom_xtick_labels,
            fontsize=10, rotation=90
        )

        custom_ytick_locations = list(range(num_markers))
        custom_ytick_labels = table.index
        plt.yticks(
            [i + 0.5 for i in custom_ytick_locations[::-1]], custom_ytick_labels,
            fontsize=10, rotation=0
        )

        legend_elements = []

        legend_elements.append(
            Line2D([0], [0], marker='s', color='none', label='True',
                   markerfacecolor='grey', markeredgecolor='k', lw=0.01, markersize=12)
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
