import os
import sys
import yaml
import pickle
import logging

import numpy as np

from tifffile import imread
from PIL import Image, ImageDraw

import matplotlib.pyplot as plt
from matplotlib.backends.qt_compat import QtWidgets
from qtpy.QtCore import QTimer

import napari
from magicgui import magicgui
from magicgui import widgets

from ..utils import (
    input_check, read_markers, get_filepath, marker_channel_number, napari_notification,
    single_channel_pyramid, triangulate_ellipse, reorganize_dfcolumns
)

logger = logging.getLogger(__name__)

arbitrary_selection_toggle = False
sample_index = 1


def callback(self, viewer, sample, data, initial_callback, next_widget, next_layout, arbitrary_widget, arbitrary_layout, roi_dir):

    if sample in data['Sample'].unique():

        print()

        check, markers_filepath = input_check(self)

        # read marker metadata
        markers, dna1, dna_moniker, abx_channels = read_markers(
            markers_filepath=markers_filepath, markers_to_exclude=self.markersToExclude, data=data
        )

        # clear existing channels from Napari window if they exist
        viewer.layers.clear()

        # remove next_widget dock and layout attributes from Napari viewer if they exist
        if not initial_callback:
            viewer.window.remove_dock_widget(next_widget)

            count = next_layout.count()
            for i in range(count - 1, -1, -1):
                item = next_layout.itemAt(i)
                widget = item.widget()
                if widget:
                    widget.setParent(None)

        # load polygon dictionary if it exists
        if os.path.exists(os.path.join(roi_dir, 'polygons.pkl')):
            f = open(os.path.join(roi_dir, 'polygons.pkl'), 'rb')
            polygon_dict = pickle.load(f)
            
        else:
            # create polygon dictionary
            polygon_dict = {}
        
        if self.showAbChannels:
            for e, ch in enumerate(reversed(abx_channels)):

                channel_number = marker_channel_number(markers, ch)

                # add immunomarker channels to Napari viewer
                file_path = get_filepath(self, check, sample, 'TIF')
                img, min, max = single_channel_pyramid(file_path, channel=channel_number)
                viewer.add_image(
                    img, rgb=False, blending='additive',
                    colormap='green', visible=False, name=ch,
                    contrast_limits=(min, max)
                )

        ###############################################################

        # add H&E image (SINGLE CHANNEL IMAGE)
        # for file_path in glob.glob(
        #   f'{self.inDir}/he/{sample}.*tif'):
        #     tiff = TiffFile(file_path, is_ome=False)
        #     rgb_pyramid = [
        #         zarr.open(tiff.series[0].levels[0].aszarr())[i]
        #         for i in list(range(len(tiff.series[0].levels)))
        #         ]
        #     rgb_pyramid = [
        #         DatasetView(i).lazy_transpose([1, 2, 0]) for
        #         i in rgb_pyramid
        #         ]
        #     rgb_pyramid = [da.from_zarr(z) for z in rgb_pyramid]
        #     viewer.add_image(
        #         rgb_pyramid, rgb=True, blending='additive',
        #         opacity=1.0, visible=False, name='H&E')
        # OR...

        # add H&E image (SHOW IMAGE AS SEPARATE RGB CHANNELS)
        # for file_path in glob.glob(
        #   f'{self.inDir}/he/{sample}.*tif'):
        #
        #     for color, channel in zip(
        #       ['red', 'green', 'blue'], [0, 1, 2]
        #       ):
        #
        #         img, min, max = single_channel_pyramid(
        #             file_path, channel=channel
        #             )
        #
        #         viewer.add_image(
        #             img, rgb=False, colormap=color,
        #             blending='additive', visible=False,
        #             name=f'H&E_{color}'
        #             )

        ###############################################################

        # add cell segmentation outlines to Napari viewer
        file_path = get_filepath(self, check, sample, 'SEG')
        seg, min, max = single_channel_pyramid(file_path, channel=0)
        viewer.add_image(
            seg, rgb=False, blending='additive', opacity=0.5, 
            colormap='gray', visible=False, name='segmentation'
        )

        # add DNA1 channel to Napari viewer
        file_path = get_filepath(self, check, sample, 'TIF')
        dna, min, max = single_channel_pyramid(file_path, channel=0)
        viewer.add_image(
            dna, rgb=False, blending='additive', colormap='gray', visible=True,
            name=f'{dna1}: {sample}', contrast_limits=(min, max)
        )

        # make lists of existing polygon shapes and vertices
        try:
            polygon_dict[sample]
            shapes = [polygon_dict[sample][i][0] for i in range(0, len(polygon_dict[sample]))]
            polygons = [polygon_dict[sample][i][1] for i in range(0, len(polygon_dict[sample]))]
        
        # or create empty list to store polygon (shapes, vertices)
        except KeyError:
            polygons = []

        if polygons:
            viewer.add_shapes(
                data=polygons, shape_type=shapes, ndim=2,
                face_color=[1.0, 1.0, 1.0, 0.2], edge_color=[0.0, 0.66, 1.0, 1.0],
                edge_width=10.0, name='ROI(s)'
            )
        else:
            viewer.add_shapes(
                data=None, shape_type='polygon', ndim=2,
                face_color=[1.0, 1.0, 1.0, 0.2], edge_color=[0.0, 0.66, 1.0, 1.0],
                edge_width=10.0, name='ROI(s)'
            )
        
        # apply previously defined contrast limits if they exist
        if os.path.exists(os.path.join(f'{self.outDir}/contrast/contrast_limits.yml')):

            logger.info('Reading existing contrast settings.')

            contrast_limits = yaml.safe_load(open(f'{self.outDir}/contrast/contrast_limits.yml'))

            viewer.layers[
                f'{dna1}: {sample}'].contrast_limits = (
                contrast_limits[dna1][0], contrast_limits[dna1][1])

            for ch in reversed(abx_channels):
                viewer.layers[ch].contrast_limits = (
                    contrast_limits[ch][0], contrast_limits[ch][1])
        
        # dock (or re-dock) next_widget to Napari window
        viewer.window.add_dock_widget(next_widget, name=f'Sample: {sample}', area='right')

        # remove and re-dock arbitrary_widget if it exists 
        # so next_widget appears first in Napari window
        if not initial_callback:
            viewer.window.remove_dock_widget(arbitrary_widget)
            viewer.window.add_dock_widget(
                arbitrary_widget, name='Arbitrary Sample Selection', area='right'
            )

        #######################################################################
        
        @magicgui(
            layout='horizontal',
            call_button='Apply ROI(s) and Move to Next Sample'
        )
        def next_sample(sample):

            global arbitrary_selection_toggle
            global sample_index

            # update ROI(s)
            updated_polygons = []

            # append updated polygons to polygon_dict and store
            for shape_type, roi in zip(viewer.layers[-1].shape_type, viewer.layers[-1].data):
                updated_polygons.append((shape_type, roi))

            polygon_dict[sample] = updated_polygons

            f = open(os.path.join(roi_dir, 'polygons.pkl'), 'wb')
            pickle.dump(polygon_dict, f)
            f.close()

            napari_notification(f'ROI(s) applied to sample {sample}.')

            # go to next sample
            try:
                if arbitrary_selection_toggle:
                    sample_index -= 1 

                sample = self.samplesForROISelection[sample_index]
                
                initial_callback = False
                callback(
                    self, viewer, sample, data, initial_callback, 
                    next_widget, next_layout, arbitrary_widget, arbitrary_layout,
                    roi_dir
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

        # if initial_callback:
        next_layout.addWidget(next_sample.native)
        
        #######################################################################
    
def selectROIs(data, self, args):

    check, markers_filepath = input_check(self)
    
    # read marker metadata
    markers, dna1, dna_moniker, abx_channels = read_markers(
        markers_filepath=markers_filepath, 
        markers_to_exclude=self.markersToExclude, 
        data=data
    )

    # print(markers_filepath, self.markersToExclude, data)
    # print("########")
    # print(markers)

    if self.samplesForROISelection:
        
        # check that all samples in samplesForROISelection are in dataframe
        if set(self.samplesForROISelection).issubset(data['Sample'].unique()):

            # create ROIs directory if it hasn't already
            roi_dir = os.path.join(self.outDir, 'ROIs')
            if not os.path.exists(roi_dir):
                os.makedirs(roi_dir)

            samples = self.samplesForROISelection
            sample = samples[0]
            
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

            initial_callback = True
            callback(
                self, viewer, sample, data, initial_callback, 
                next_widget, next_layout, arbitrary_widget, arbitrary_layout,
                roi_dir
            )

            viewer.window.add_dock_widget(
                arbitrary_widget, name='Arbitrary Sample Selection', area='right'
            )
            
            viewer.scale_bar.visible = True
            viewer.scale_bar.unit = 'um'

            napari.run()

            print()

            ###################################################################
            # load polygon dictionary if it exists

            if os.path.exists(os.path.join(roi_dir, 'polygons.pkl')):
                f = open(os.path.join(roi_dir, 'polygons.pkl'), 'rb')
                polygon_dict = pickle.load(f)

            else:
                logger.info(
                    'Aborting; ROI polygon dictionary does not exist. '
                    'Please select ROIs.'
                )
                sys.exit()
            
            idxs_to_drop = {}
            for sample in samples:
                try:
                    if polygon_dict[sample]:

                        logger.info(f'Generating ROI mask(s) for sample: {sample}')

                        sample_data = data[['X_centroid', 'Y_centroid', 'CellID']][
                            data['Sample'] == sample].astype(int)

                        sample_data['tuple'] = list(
                            zip(sample_data['X_centroid'],
                                sample_data['Y_centroid'])
                        )
 
                        file_path = get_filepath(self, check, sample, 'TIF')
                        dna = imread(file_path, key=0)

                        # create pillow image to convert into boolean mask
                        img = Image.new('L', (dna.shape[1], dna.shape[0]))

                        for shape_type, verts in polygon_dict[sample]:

                            selection_verts = np.round(verts).astype(int)

                            if shape_type == 'ellipse':

                                vertices, triangles = triangulate_ellipse(
                                    selection_verts
                                )

                                # flip 2-tuple coordinates returned by
                                # triangulate_ellipse() to draw image mask
                                vertices = [tuple(reversed(tuple(i))) for i in vertices]

                                # update pillow image with polygon
                                ImageDraw.Draw(img).polygon(
                                    vertices, outline=1, fill=1
                                )

                            else:
                                vertices = list(tuple(
                                    zip(selection_verts[:, 1],
                                        selection_verts[:, 0])
                                ))

                                # update pillow image with polygon
                                ImageDraw.Draw(img).polygon(vertices, outline=1, fill=1)

                        # convert pillow image into boolean numpy array
                        mask = np.array(img, dtype=bool)

                        # use numpy fancy indexing to get centroids
                        # where boolean mask is True
                        xs, ys = zip(*sample_data['tuple'])

                        inter = mask[ys, xs]

                        # update sample_data with boolean calls per cell
                        sample_data['inter'] = inter

                        if self.delintMode is True:
                            idxs_to_drop[sample] = list(
                                sample_data['CellID'][sample_data['inter']]
                            )
                        else:
                            idxs_to_drop[sample] = list(
                                sample_data['CellID'][~sample_data['inter']]
                            )
                    else:
                        logger.info(f'No ROIs selected for sample: {sample}')
                        idxs_to_drop[sample] = []

                except KeyError:
                    logger.info(
                        f'Aborting; ROIs have not been ' 
                        f'selected for sample {sample}. '
                        'Please re-run selectROIs module to select '
                        'ROIs for this sample.'
                    )
                    sys.exit()
                
            print()

            # drop cells from samples
            if not all([True if not value else False for value in idxs_to_drop.values()]):
                for sample, cell_ids in idxs_to_drop.items():
                    if cell_ids:
                        logger.info(f'Dropping cells from sample: {sample}')
                        global_idxs_to_drop = data[
                            (data['Sample'] == sample) 
                            & (data['CellID'].isin(set(cell_ids)))].index
                        data.drop(global_idxs_to_drop, inplace=True)
                    else:
                        pass
                print()

            # save plots of selected data points
            plot_dir = os.path.join(roi_dir, 'plots')
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

            for sample in samples:

                logger.info(f'Plotting ROI selections for sample: {sample}')

                file_path = get_filepath(self, check, sample, 'TIF')
                dna = imread(file_path, key=0)

                fig, ax = plt.subplots()
                ax.imshow(dna, cmap='gray')
                ax.grid(False)
                ax.set_axis_off()
                coords = data[['X_centroid', 'Y_centroid', 'Area']][
                    data['Sample'] == sample]
                ax.scatter(
                    coords['X_centroid'], coords['Y_centroid'], s=0.35, lw=0.0, c='yellow'
                )
                plt.title(f'Sample {sample}', size=10)
                plt.tight_layout()
                plt.savefig(os.path.join(plot_dir, f'{sample}.png'), dpi=800)
                plt.close('all')

        else:
            logger.info(
                'Aborting; one or more samples in "samplesForROISelection" configuration ' 
                'parameter are not in the dataframe.'
            )
            sys.exit()
            
        # create ROIs directory if it hasn't already
        roi_dir = os.path.join(self.outDir, 'ROIs')
        if not os.path.exists(roi_dir):
            os.makedirs(roi_dir)

        samples = iter(self.samplesForROISelection)
        sample = next(samples)
        global last_sample
        last_sample = None
        
        viewer = napari.Viewer(title='CyLinter')
        viewer.scale_bar.visible = True
        viewer.scale_bar.unit = 'um'

        ### Load data for the data layer(s), i.e., polygon dicts, and potentially
        ### the points corresponding to centroids of cells classified as artifacts.
        varname_filename_lst = [('ROI', 'polygons.pkl'),
                            ('ROI2', 'polygons2.pkl'),
                            ('Detected Artifacts', 'points.pkl')]
        layer_type = {'ROI': 'shape', 
                      'ROI2': 'shape',
                      'Detected Artifacts': 'point'}
        extra_layers = {}
        for varname, fname in varname_filename_lst:
            if os.path.exists(os.path.join(roi_dir, fname)):
                f = open(os.path.join(roi_dir, fname), 'rb')
                extra_layers[varname] = pickle.load(f)
            else:
                extra_layers[varname] = {}
        
        ###################################################################

        def add_layers(sample):
            # antibody/immunomarker channels
            if self.showAbChannels:
                for _, ch in enumerate(reversed(abx_channels)):
                    channel_number = marker_channel_number(markers, ch)
                    file_path = get_filepath(self, check, sample, 'TIF')
                    img, min, max = single_channel_pyramid(file_path, channel=channel_number)
                    viewer.add_image(
                        img, rgb=False, blending='additive',
                        colormap='green', visible=False, name=ch,
                        contrast_limits=(min, max)
                    )
            
            # H&E channel, as single image or using separate RGB channels, to be implemented 

            # cell segmentation outlines channel
            file_path = get_filepath(self, check, sample, 'SEG')
            seg, min, max = single_channel_pyramid(file_path, channel=0)
            viewer.add_image(
                seg, rgb=False, blending='additive', opacity=0.5, 
                colormap='gray', visible=False, name='segmentation'
            )

            # DNA1 channel
            file_path = get_filepath(self, check, sample, 'TIF')
            dna, min, max = single_channel_pyramid(file_path, channel=0)
            viewer.add_image(
                dna, rgb=False, blending='additive', colormap='gray', visible=True,
                name=f'{dna1}: {sample}', contrast_limits=(min, max)
            )

            # ROI selection channel, as well as ROI2 and labeled artifacts
            for varname, layer_data in extra_layers.items():
                if layer_type[varname] == 'shape':
                    try:
                        shapes = [shape_data[0] for shape_data in layer_data[sample]]
                        polygons = [shape_data[1] for shape_data in layer_data[sample]]
                    except KeyError:
                        shapes = 'polygon'
                        polygons = None
                    
                    viewer.add_shapes(
                            data=polygons, shape_type=shapes, ndim=2,
                            face_color=[1.0, 1.0, 1.0, 0.2], edge_color=[0.0, 0.66, 1.0, 1.0],
                            edge_width=10.0, name=varname
                    )
                elif layer_type[varname] == 'point':
                    try:
                        points = layer_data[sample]
                    except:
                        points = None
                    viewer.add_points(points, ndim=2, face_color=[1.0, 0, 0, 1.0], 
                                      edge_color=[0.0, 0.0, 0.0, 0.0],
                                      edge_width=0.0, name=varname, size=10.0)
            

            ### Apply previously defined contrast limits if they exist
            if os.path.exists(os.path.join(f'{self.outDir}/contrast/contrast_limits.yml')):

                logger.info('Reading existing contrast settings.')

                contrast_limits = yaml.safe_load(open(f'{self.outDir}/contrast/contrast_limits.yml'))

                viewer.layers[
                    f'{dna1}: {sample}'].contrast_limits = (
                    contrast_limits[dna1][0], contrast_limits[dna1][1])

                for ch in reversed(abx_channels):
                    viewer.layers[ch].contrast_limits = (
                        contrast_limits[ch][0], contrast_limits[ch][1])
                    
        def clear_widgets(widgets):
            for widget in widgets:
                viewer.window.remove_dock_widget(widget)

        def add_widgets(last_sample, sample):
            @magicgui(
                call_button="Apply ROI(s) and Move to Next Sample"
            )
            def next_sample(last_sample, sample, widgets):
                try: 
                    # update ROI and other shape layers
                    for i, tupl in enumerate(varname_filename_lst[::-1]):
                        varname, filename = tupl
                        updated_layer_data = []
                        layer = viewer.layers[-(i+1)]
                        if layer_type[varname] == 'shape':
                            for shape_type, roi in zip(layer.shape_type, layer.data):
                                updated_layer_data.append((shape_type, roi))
                            extra_layers[varname][sample] = updated_layer_data
                        elif layer_type[varname] == 'point':
                            updated_layer_data = layer.data
                            extra_layers[varname][sample] = updated_layer_data

                        f = open(os.path.join(roi_dir, filename), 'wb')
                        pickle.dump(extra_layers[varname], f)
                        f.close()

                    napari_notification(f'ROI(s) applied to sample {sample}.')

                    # move to next sample and re-render the GUI
                    if last_sample is None:
                        sample = next(samples)
                    else: 
                        sample = last_sample
                        last_sample = None
                    next_sample.sample.bind(sample)
                    viewer.layers.clear()
                    clear_widgets(widgets)
                    add_layers(sample)
                    add_widgets(last_sample, sample)
                except StopIteration:
                    print()
                    napari_notification('Gating complete!')
                    QTimer().singleShot(0, viewer.close)

            next_sample.sample.bind(sample) #should bind first sample_id
            next_sample.last_sample.bind(last_sample)

            @magicgui(
                call_button='Enter',
                next_sample={'label': 'Sample Name'}
            )
            def arbitrary_sample(sample, next_sample: str, widgets):
                global last_sample
                if last_sample is None:
                    last_sample = sample
                viewer.layers.clear()
                clear_widgets(widgets)
                add_layers(next_sample)
                add_widgets(last_sample, next_sample)
            arbitrary_sample.sample.bind(sample)


            @magicgui(
                call_button="Auto label artifacts"
            )
            def label_artifacts():
                viewer.layers[-1].data = None
                global artifact_mask
                artifact_mask = [True] * len(data)
                centroids = data[['Y_centroid', 'X_centroid']][artifact_mask]
                viewer.layers[-1].add(centroids)
            


            widgets = [next_sample, arbitrary_sample, label_artifacts]
            widget_names = [f'Sample {sample}', 'Arbitrary Sample Selection', 'Automation Module']
            napari_widgets = []
            for widget, widget_name in zip(widgets, widget_names):
                napari_widgets.append(viewer.window.add_dock_widget(widget=widget,
                                                                    add_vertical_stretch=False,
                                                                    name=widget_name))
            next_sample.widgets.bind(napari_widgets)  
            arbitrary_sample.widgets.bind(napari_widgets) 

            for napari_widget in napari_widgets:
                napari_widget.widget().setSizePolicy(
                    QtWidgets.QSizePolicy.Minimum,
                    QtWidgets.QSizePolicy.Maximum,
                )

        ###################################################################

        add_layers(sample)
        add_widgets(last_sample, sample)

        napari.run() # blocks until window is closed

        ###################################################################

        idxs_to_drop = {}
        for sample in samples:
            try:
                if extra_layers['ROI'][sample]:

                    logger.info(f'Generating ROI mask(s) for sample: {sample}')

                    sample_data = data[['X_centroid', 'Y_centroid', 'CellID']][
                        data['Sample'] == sample].astype(int)

                    sample_data['tuple'] = list(
                        zip(sample_data['X_centroid'],
                            sample_data['Y_centroid'])
                    )

                    file_path = get_filepath(self, check, sample, 'TIF')
                    dna = imread(file_path, key=0)

                    # create pillow image to convert into boolean mask
                    def shape_layer_to_mask(layer_data):
                        img = Image.new('L', (dna.shape[1], dna.shape[0]))

                        for shape_type, verts in layer_data:

                            selection_verts = np.round(verts).astype(int)

                            if shape_type == 'ellipse':

                                vertices, triangles = triangulate_ellipse(
                                    selection_verts
                                )

                                # flip 2-tuple coordinates returned by
                                # triangulate_ellipse() to draw image mask
                                vertices = [tuple(reversed(tuple(i))) for i in vertices]

                                # update pillow image with polygon
                                ImageDraw.Draw(img).polygon(
                                    vertices, outline=1, fill=1
                                )

                            else:
                                vertices = list(tuple(
                                    zip(selection_verts[:, 1],
                                        selection_verts[:, 0])
                                ))

                                # update pillow image with polygon
                                ImageDraw.Draw(img).polygon(vertices, outline=1, fill=1)

                        # convert pillow image into boolean numpy array
                        mask = np.array(img, dtype=bool)
                        return mask
                    
                    ROI_mask = shape_layer_to_mask(extra_layers['ROI'][sample])
                    ROI2_mask = shape_layer_to_mask(extra_layers['ROI2'][sample])
                    # use numpy fancy indexing to get centroids
                    # where boolean mask is True
                    xs, ys = zip(*sample_data['tuple'])

                    inter1 = ROI_mask[ys, xs]
                    global artifact_mask
                    inter2 = ~ROI2_mask[ys, xs] * artifact_mask[ys, xs]

                    # update sample_data with boolean calls per cell
                    sample_data['inter1'] = inter1
                    sample_data['inter2'] = inter2

                    if self.delintMode is True:
                        idxs_to_drop[sample] = list(
                            sample_data['CellID'][sample_data['inter1'] * sample_data['inter2']]
                        )
                    else:
                        idxs_to_drop[sample] = list(
                            sample_data['CellID'][~sample_data['inter1'] * sample_data['inter2']]
                        )
                else:
                    logger.info(f'No ROIs selected for sample: {sample}')
                    idxs_to_drop[sample] = []

            except KeyError:
                logger.info(
                    f'Aborting; ROIs have not been ' 
                    f'selected for sample {sample}. '
                    'Please re-run selectROIs module to select '
                    'ROIs for this sample.'
                )
                sys.exit()
        print()

        # drop cells from samples
        for sample, cell_ids in idxs_to_drop.items():
            if cell_ids:
                logger.info(f'Dropping cells from sample: {sample}')
                global_idxs_to_drop = data[
                    (data['Sample'] == sample) 
                    & (data['CellID'].isin(set(cell_ids)))].index
                data.drop(global_idxs_to_drop, inplace=True)
            else:
                pass
        print()

        # save plots of selected data points
        plot_dir = os.path.join(roi_dir, 'plots')
        if not os.path.exists(plot_dir):
            os.mkdir(plot_dir)

        for sample in samples:

            logger.info(f'Plotting ROI selections for sample: {sample}')

            for file_path in glob.glob(f'{self.inDir}/tif/{sample}.*tif'):
                dna = imread(file_path, key=0)

            fig, ax = plt.subplots()
            ax.imshow(dna, cmap='gray')
            ax.grid(False)
            ax.set_axis_off()
            coords = data[['X_centroid', 'Y_centroid', 'Area']][
                data['Sample'] == sample]
            ax.scatter(
                coords['X_centroid'], coords['Y_centroid'], s=0.35, lw=0.0, c='yellow'
            )
            plt.title(f'Sample {sample}', size=10)
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f'{sample}.png'), dpi=800)
            plt.close('all')

    else:
        logger.info(
            'Skipping ROI selection, no samples for ROI selection specified in config.yml'
        )
    
    data = reorganize_dfcolumns(data, markers, self.dimensionEmbedding)

    print()
    print()
    return data
