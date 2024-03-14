import os
import sys
import yaml
import pickle
import logging
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

import cv2
from PIL import Image, ImageDraw

import matplotlib.pyplot as plt
from qtpy.QtCore import QTimer
from matplotlib.backends.qt_compat import QtWidgets

import napari
from magicgui import magicgui
from skimage.morphology import flood
from magicgui.widgets import ComboBox, SpinBox, Container, Button, CheckBox

from ..utils import (
    input_check, read_markers, get_filepath, marker_channel_number, sort_qc_report,
    single_channel_pyramid, triangulate_ellipse, reorganize_dfcolumns,
    upscale, ArtifactInfo, artifact_detector_v3
)

logger = logging.getLogger(__name__)


@dataclass
class GlobalState():
    loaded_ims = {}
    abx_layers = {}
    base_clf = None
    last_sample = None
    artifacts = {}
    _artifact_mask: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=int))
    binarized_artifact_mask: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=bool))
    artifact_proba: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=float))
    artifact_detection_threshold: float = 0.5
 
    current_layer: napari.layers.Layer = None
    current_point: np.ndarray = None
    current_tol: int = None

    @property
    def artifact_mask(self):
        return self._artifact_mask

    @artifact_mask.setter
    def artifact_mask(self, data: np.ndarray):
        self._artifact_mask = data
        # may need to remove this hard-coded assumption later
        self.binarized_artifact_mask = self.artifact_mask != 1


def selectROIs(data, self, args):

    print()
    
    global_state = GlobalState()

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
        if 'selectROIs' not in qc_report or qc_report['selectROIs'] is None:
            qc_report['selectROIs'] = {}
            reload_report = True
        if reload_report:
            qc_report_sorted = sort_qc_report(qc_report, module='selectROIs', order=None)
            f = open(report_path, 'w')
            yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
            qc_report = yaml.safe_load(open(report_path))
    except:
        logger.info(
            'Aborting; QC report missing from CyLinter output directory. Re-start pipeline '
            'from aggregateData module to initialize QC report.'
        )
        sys.exit()
    
    if self.samplesForROISelection:

        # check that all samples in samplesForROISelection are in dataframe
        if not (set(self.samplesForROISelection).issubset(data['Sample'].unique())):
            logger.info(
                'Aborting; one or more samples in "samplesForROISelection" configuration '
                'parameter are not in the dataframe.'
            )
            sys.exit()

        # create manual ROIs directory if it hasn't already
        roi_dir = os.path.join(self.outDir, 'ROIs')
        if not os.path.exists(roi_dir):
            os.makedirs(roi_dir)
        
        # create automated ROIs directory if it hasn't already
        art_dir = os.path.join(roi_dir, 'masks', self.artifactDetectionMethod)
        if not os.path.exists(art_dir):
            os.makedirs(art_dir)

        samples = iter(self.samplesForROISelection)
        sample = next(samples)

        viewer = napari.Viewer(title='CyLinter')
        viewer.scale_bar.visible = True
        viewer.scale_bar.unit = 'um'

        ###################################################################
        # Load data for the data layer(s), i.e., polygon dicts, and potentially
        # the points corresponding to centroids of cells classified as artifacts.
        def load_extra_layers():
            
            varname_filename_lst = [
                ('ROI', 'Manual ROI Selections (neg.)' if self.delintMode else
                 'Manual ROI Selections (pos.)')
            ]
            layer_type = [('ROI', 'shape')]
            layer_name = [
                ('ROI', 'Manual ROI Selections (neg.)' if self.delintMode else
                 'Manual ROI Selections (pos.)')
            ]
            
            if self.autoArtifactDetection:
                
                if self.artifactDetectionMethod == 'MLP':
                    varname_filename_lst += [
                        ('ROI2', 'artifact_pred_selection.pkl'),
                        ('Detected Artifacts', 'points.pkl')
                    ]
                    layer_type += [('ROI2', 'shape'), ('Detected Artifacts', 'point')]
                    layer_name += [
                        ('ROI2', 'Ignore Artifact Prediction Selection'), 
                        ('Detected Artifacts', 'Predicted Artifacts')
                    ]
                
                elif self.artifactDetectionMethod == 'classical':
                    
                    # to avoid complication, we implement the most straightforward method
                    # of registering every abx channel, even if some may not have artifacts
                    for abx_channel in abx_channels:
                        varname_filename_lst += [
                            (f'{abx_channel}_seeds', f'{abx_channel}_artifact_seeds.pkl'),
                            (f'{abx_channel}_mask', f'{abx_channel}_artifact_mask.pkl')
                        ]
                        layer_type += [
                            (f'{abx_channel}_seeds', 'point'), (f'{abx_channel}_mask', 'image') 
                        ]
                        layer_name += [
                            (f'{abx_channel}_seeds', f'{abx_channel} Artifact Seeds'),
                            (f'{abx_channel}_mask', f'{abx_channel} Artifact Mask')  
                        ]
            
            layer_type = dict(layer_type)
            layer_name = dict(layer_name)

            extra_layers = {}
            for varname, fname in varname_filename_lst: 

                if varname == 'ROI':
                    try: 
                        extra_layers[varname] = qc_report['selectROIs'][fname]
                    except KeyError:
                        extra_layers[varname] = {}
                else:
                    try:
                        fdir = os.path.join(art_dir, fname)
                        f = open(fdir, 'rb')
                        extra_layers[varname] = pickle.load(f)
                    except FileNotFoundError:
                        extra_layers[varname] = {}
            
            return extra_layers, layer_type, layer_name, varname_filename_lst

        def float_roi_layer_to_top():
            roi_layer_id = viewer.layers.index(viewer.layers[layer_name['ROI']])
            top_layer_id = len(viewer.layers)
            viewer.layers.move(roi_layer_id, top_layer_id)
            viewer.layers.select_next()  # ensuring manual ROI layer is highlighted

        def add_layers(sample):
            
            # reset global state if using classical artifact detection
            if self.artifactDetectionMethod == 'classical':
                global_state.artifacts = {}
                global_state.loaded_ims = {}
                global_state.abx_layers = {}
                for abx_channel in abx_channels:
                    try:
                        global_state.artifacts[abx_channel] = (
                            extra_layers[f'{abx_channel}_mask'][sample]
                        )
                    except:
                        pass

            # antibody channels
            if self.showAbChannels:
                for ch in reversed(abx_channels):
                    channel_number = marker_channel_number(self, markers, ch)
                    file_path = get_filepath(self, check, sample, 'TIF')
                    img, min, max = single_channel_pyramid(
                        file_path, channel=channel_number)
                    layer = viewer.add_image(
                        img, rgb=False, blending='additive',
                        colormap='green', visible=False, name=ch,
                        contrast_limits=(min, max)
                    )
                    global_state.loaded_ims[ch] = img
                    global_state.abx_layers[ch] = layer

            # H&E channel (single image or separate RGB channels), to be implemented here

            # cell segmentation outlines channel
            file_path = get_filepath(self, check, sample, 'SEG')
            seg, min, max = single_channel_pyramid(file_path, channel=0)
            viewer.add_image(
                seg, rgb=False, blending='additive', opacity=1.0,
                colormap='gray', visible=False, name='segmentation',
                contrast_limits=(min, max)
            )

            # DNA1 channel
            file_path = get_filepath(self, check, sample, 'TIF')
            channel_number = marker_channel_number(self, markers, self.counterstainChannel)
            dna, min, max = single_channel_pyramid(file_path, channel=channel_number)
            viewer.add_image(
                dna, rgb=False, blending='additive', colormap='gray', visible=True,
                name=self.counterstainChannel, contrast_limits=(min, max)
            )

            # ROI selection channel, as well as ROI2 and labeled artifacts
            for varname, layer_data in extra_layers.items():
                
                if layer_type[varname] == 'shape':
                    try:
                        shapes = [shape_data[0] for shape_data in layer_data[sample]]
                        polygons = [np.array(shape_data[-1]) for shape_data in layer_data[sample]]
                    except KeyError:
                        shapes = 'polygon'
                        polygons = None

                    if varname == 'ROI':
                        if self.delintMode:
                            edge_color = [1.0, 0.00, 0.0, 1.0]
                        else:
                            edge_color = [0.0, 0.66, 1.0, 1.0]
                    
                    elif varname == 'ROI2':
                        edge_color = [0.0, 0.66, 1.0, 1.0]

                    viewer.add_shapes(
                        data=polygons, shape_type=shapes, ndim=2,
                        face_color=[1.0, 1.0, 1.0, 0.3], edge_color=edge_color,
                        edge_width=5.0, name=layer_name[varname]
                    )
                
                elif layer_type[varname] == 'image':
                    if self.autoArtifactDetection and self.artifactDetectionMethod == 'classical': 
                        abx_channel = varname.split('_')[0]
                        
                        try:
                            extra_layers[varname][sample].render_seeds(
                                viewer, global_state.loaded_ims, layer_name, abx_channel
                            )
                        except:
                            pass 
                        
                        try:
                            extra_layers[varname][sample].render_mask(
                                viewer, global_state.loaded_ims, layer_name, abx_channel
                            )
                        except:
                            pass

                    elif self.autoArtifactDetection and self.artifactDetectionMethod == 'MLP':
                        try:
                            (points, global_state.artifact_mask,
                             global_state.artifact_proba) = layer_data[sample]
                            artifact_mask_ = (
                                global_state.artifact_mask[
                                    global_state.binarized_artifact_mask]
                            )
                            artifact_proba_ = (
                                global_state.artifact_proba[
                                    global_state.binarized_artifact_mask]
                            )
                        except:  # may want to be explicit here about the expected exception.
                            points = None
                            artifact_mask_ = []
                            artifact_proba_ = []
                        
                        viewer.add_points(
                            points, ndim=2, edge_color=[0.0, 0.0, 0.0, 0.0],
                            edge_width=0.0, name=layer_name[varname], size=10.0,
                            face_color_cycle={
                                1: 'white', 2: 'red', 3: 'blue',
                                4: 'green', 5: 'cyan', 6: 'magenta'},
                            face_color='artifact_class',
                            features={'artifact_class': np.array(
                                artifact_mask_, dtype=int)}
                        )  # face_color=[1.0, 0, 0, 0.2],
                        points_layer = viewer.layers[-1]
                        points_layer.face_color_mode = 'direct'
                        points_layer.face_color[:, -1] * \
                            np.array(artifact_proba_)
                        points_layer.refresh()
            
            # apply previously defined contrast limits if they exist
            try:
                viewer.layers[
                    f'{self.counterstainChannel}'].contrast_limits = (
                    qc_report['setContrast'][self.counterstainChannel][0],
                    qc_report['setContrast'][self.counterstainChannel][1])

                for ch in reversed(abx_channels):
                    viewer.layers[ch].contrast_limits = (
                        qc_report['setContrast'][ch][0], qc_report['setContrast'][ch][1]
                    ) 
                logger.info('Existing contrast settings applied.')
            except:
                pass
            
            float_roi_layer_to_top()

        def save_shapes(sample):
            
            if self.artifactDetectionMethod == 'MLP':
                for i, tupl in enumerate(varname_filename_lst[::-1]):
                    varname, filename = tupl
                    updated_layer_data = []
                    layer = viewer.layers[-(i + 1)]
                    if layer_type[varname] == 'shape':
                        for shape_type, roi in zip(layer.shape_type, layer.data):
                            updated_layer_data.append([shape_type, round(roi.tolist())])
                        extra_layers[varname][sample] = updated_layer_data
                    elif layer_type[varname] == 'point':
                        updated_layer_data = (
                            layer.data, global_state.artifact_mask, global_state.artifact_proba
                        )
                        extra_layers[varname][sample] = updated_layer_data
                    f = open(os.path.join(art_dir, filename), 'wb')
                    pickle.dump(extra_layers[varname], f)
                    f.close()
            
            elif self.artifactDetectionMethod == 'classical':
                for varname, filename in varname_filename_lst:
                    save_layer = False
                    
                    if layer_type[varname] == 'shape':
                        layer = viewer.layers[layer_name[varname]]
                        updated_layer_data = []
                        for shape_type, roi in zip(layer.shape_type, layer.data):
                            updated_layer_data.append(
                                [shape_type, np.round(roi).astype('int').tolist()]
                            )
                        extra_layers[varname][sample] = updated_layer_data
                        save_layer = True
                    
                    elif layer_type[varname] == 'image':
                        abx_channel = varname.split('_')[0]
                        artifact_info = global_state.artifacts.get(abx_channel)
                        try:
                            viewer.layers.index(layer_name[varname])
                        except:
                            if artifact_info is not None:  # napari layer was deleted 
                                extra_layers[varname].pop(sample, None)
                                artifact_info = None
                                try:
                                    os.remove(os.path.join(art_dir, filename))
                                except:
                                    pass  # file does not exist yet
                                # continue
                        
                        if artifact_info is not None:
                            extra_layers[varname][sample] = artifact_info
                            artifact_info.features = None
                            artifact_info.artifact_layer = None
                            artifact_info.seed_layer = None
                        if len(extra_layers[varname]) > 0:
                            save_layer = True
                    
                    if save_layer:
                        if layer_type[varname] == 'shape':

                            # empty ROI lists can occur when ROIs are deleted in Napari, delete 
                            # key entries in extra_layers[varname] with no associated ROIs
                            extra_layers[varname] = {
                                key: value for key, value in extra_layers[varname].items()
                                if value
                            }

                            qc_report['selectROIs'][layer_name[varname]] = extra_layers[varname]
                            if not bool(extra_layers[varname]):
                                del qc_report['selectROIs'][layer_name[varname]]

                        elif layer_type[varname] == 'image':
                            f = open(os.path.join(art_dir, filename), 'wb')
                            pickle.dump(extra_layers[varname], f)
                            f.close()

        def add_widgets(last_sample, sample):
            
            @magicgui(call_button="Apply ROI(s) and Move to Next Sample")
            def next_sample(qc_report, last_sample, sample, widgets):
                
                try:
                    # update ROI and other shape layers
                    save_shapes(sample)
                    
                    # compute vertices for automatically generated artifact channel masks
                    channel_mask_verts = []
                    for varname, filename in varname_filename_lst: 
                        if layer_type[varname] == 'image':
                            
                            abx_channel = varname.split('_')[0]
                            artifact_info = extra_layers[varname].get(sample)

                            if artifact_info is None:
                                try:
                                    if bool(qc_report['selectROIs'][
                                            'Automated ROI Selections (neg.)'][sample]):
                                        del qc_report['selectROIs'][
                                            'Automated ROI Selections (neg.)'][sample]
                                except KeyError:
                                    pass
                                else:
                                    continue
                            else:
                                upscaled_mask = upscale(
                                    artifact_info.mask > 0, 
                                    global_state.loaded_ims[abx_channel][0]
                                )

                                # save polygon vertices for automated ROI masks
                                # find contours in the binary mask
                                contours, _ = cv2.findContours(
                                    upscaled_mask.astype(np.uint8),
                                    cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
                                )

                                # compute mask vertices for current sample
                                for contour in contours:
                                    contour = np.squeeze(contour, axis=1)
                                    contour = np.flip(contour, axis=-1)  # inverting x/y coords
                                    channel_mask_verts.append(
                                        ['polygon', abx_channel, contour.tolist()]
                                    )

                    if 'Automated ROI Selections (neg.)' not in qc_report['selectROIs'].keys():
                        qc_report['selectROIs']['Automated ROI Selections (neg.)'] = {}

                    if len(channel_mask_verts) > 0:
                        qc_report['selectROIs'][
                            'Automated ROI Selections (neg.)'][sample] = channel_mask_verts
                    
                    if not bool(qc_report['selectROIs']['Automated ROI Selections (neg.)']):
                        del qc_report['selectROIs']['Automated ROI Selections (neg.)']

                    # if all ROIs were deleted
                    if not bool(qc_report['selectROIs']):
                        del qc_report['selectROIs']

                    # store automated and manually-defined polygon vertices (added to
                    # qc_report dict by save_shapes function) to CyLinter QC report
                    qc_report_sorted = sort_qc_report(qc_report, module='selectROIs', order=None)
                    f = open(report_path, 'w')
                    yaml.dump(qc_report_sorted, f, sort_keys=False, allow_unicode=False)
                    
                    # move to next sample and re-render the GUI
                    if last_sample is None:
                        sample = next(samples)
                    else:
                        sample = last_sample
                        global_state.last_sample = None
                    next_sample.sample.bind(sample)
                    viewer.layers.clear()
                    clear_widgets(widgets)
                    add_layers(sample)
                    add_widgets(last_sample, sample)
                    logger.info(f'ROI(s) applied to sample {sample}.')
                    print()
                
                except StopIteration:
                    logger.info('ROI selection complete')
                    QTimer().singleShot(0, viewer.close)
                    print()

            next_sample.sample.bind(sample)  # should bind first sample_id
            next_sample.last_sample.bind(global_state.last_sample)
            next_sample.qc_report.bind(qc_report)

            @magicgui(
                call_button='Enter',
                next_sample={'label': 'Sample Name'}
            )
            def arbitrary_sample(sample, next_sample: str, widgets):
                
                if next_sample not in self.samplesForROISelection:
                    napari.utils.notifications.show_error(
                        'Sample ID not in samplesForROISelection.'
                    )
                    print()
                    return

                if global_state.last_sample is None:
                    global_state.last_sample = sample
                save_shapes(sample)
                viewer.layers.clear()
                clear_widgets(widgets)
                add_layers(next_sample)
                add_widgets(global_state.last_sample, next_sample)
            arbitrary_sample.sample.bind(sample)

            @magicgui(
                proba_threshold={'label': 'Threshold',
                                 'widget_type': 'FloatSlider',
                                 'min': 0.001,
                                 'max': 0.999,
                                 'step': 0.001},
                call_button="Auto label artifacts"
            )
            def label_artifacts_MLP(proba_threshold: float = 0.5):
                
                viewer.layers.selection.active = viewer.layers[-1]
                viewer.layers[-1].data = None
                global_state.artifact_detection_threshold = proba_threshold
                global_state.artifact_mask, class_probas = artifact_detection_model_MLP(data)
                global_state.artifact_proba = 1 - class_probas[:, 0]
                artifact_mask, binarized_artifact_mask, artifact_proba = (
                    global_state.artifact_mask, 
                    global_state.binarized_artifact_mask,
                    global_state.artifact_proba
                )
                centroids = data[['Y_centroid', 'X_centroid']][binarized_artifact_mask]
                points_layer = viewer.layers[-1]
                points_layer.face_color_mode = 'cycle'
                points_layer.add(centroids)
                points_layer.features.loc[-len(centroids):, 'artifact_class'] = (
                    artifact_mask[binarized_artifact_mask]
                )
                points_layer.refresh_colors()
                points_layer.face_color_mode = 'direct'
                points_layer.face_color[:, -1] = np.array(
                    artifact_proba[binarized_artifact_mask])
                points_layer.refresh()
                viewer.layers.selection.active = viewer.layers[-2]

            @label_artifacts_MLP.proba_threshold.changed.connect
            def on_slider_changed(threshold):
                
                # consider using resource management "with...as" here
                viewer.layers.selection.active = viewer.layers[-1]
                points_layer = viewer.layers[-1]
                global_state.artifact_detection_threshold = threshold
                points_layer.face_color_mode = 'direct'
                proba = np.array(
                    global_state.artifact_proba[global_state.binarized_artifact_mask])
                points_layer.face_color[:, -1] = np.piecewise(
                    proba, [proba < threshold, proba > threshold],
                    [lambda v: 0, lambda v: (v - threshold) / (1 - threshold) * (1 - 0.3) + 0.3]
                )
                # np.clip((proba - threshold) / (np.max(proba) - threshold), 0.001, 0.999)
                points_layer.refresh()
                viewer.layers.selection.active = viewer.layers[-2]
            
            # might want to rewrite the above code for MLP model to reduce
            # overhead of defining these widgets if we are never going to use them
            def generate_widgets_for_classical_method():
                
                channel_selector_dropdown = ComboBox(choices=abx_channels, label='Abx channel')
                sensitivity_spinbox = SpinBox(value=0, label='Sensitivity', min=0, max=255)
                sensitivity_auto_checkbox = CheckBox(value=True, label='Auto')
                compute_mask_button = Button(label='Compute Artifact Mask')
                widget_combo_1 = Container(widgets=[
                    channel_selector_dropdown, 
                    sensitivity_spinbox,
                    sensitivity_auto_checkbox,
                    compute_mask_button
                ])

                tolerance_spinbox = SpinBox(value=0, label='Tolerance')
                widget_combo_2 = Container(widgets=[
                    tolerance_spinbox, 
                ])

                # bind callbacks
                artifacts = global_state.artifacts
                loaded_ims = global_state.loaded_ims

                @compute_mask_button.clicked.connect
                def compute_artifact_mask():
                    
                    # first remove any existing layers if exist
                    abx_channel = channel_selector_dropdown.value
                    if abx_channel in artifacts.keys():
                        try:
                            viewer.layers.remove(
                                viewer.layers[layer_name[f'{abx_channel}_seeds']]
                            )
                            viewer.layers.remove(
                                viewer.layers[layer_name[f'{abx_channel}_mask']]
                            )
                        except:
                            pass
                    
                    # next, compute
                    params = {'downscale': 2}
                    if sensitivity_auto_checkbox.value:
                        h = None
                    else:
                        h = 255 - sensitivity_spinbox.value
                    
                    artifact_mask, im_transformed, seeds, tols, opt_h = \
                        artifact_detector_v3(
                            loaded_ims[abx_channel], downscale=params['downscale'], h=h
                        )
                    sensitivity_spinbox.value = 255 - opt_h
                    sensitivity_auto_checkbox.value = True
                    artifact_info = ArtifactInfo(
                        params, artifact_mask, im_transformed, 
                        dict(zip(range(len(seeds)), seeds)), tols
                    )
                    artifacts[abx_channel] = artifact_info
                    artifact_info.render(viewer, loaded_ims, layer_name, abx_channel)
                    float_roi_layer_to_top()

                    seed_layer, artifact_layer = (
                        artifact_info.seed_layer, artifact_info.artifact_layer
                    )
                    artifact_info.bind_listener_seeds(viewer, global_state, tolerance_spinbox)
                    for layer in viewer.layers:
                        layer.visible = False
                    global_state.abx_layers[abx_channel].visible = True
                    artifact_layer.visible = True
                    seed_layer.visible = False
                    
                    # this does not seem to work to highlight a specific layer
                    # artifact_layer.selected = True 
                    
                    # this does work to move up and down from the current highlighted layer 
                    # viewer.layers.select_previous()  

                @sensitivity_spinbox.changed.connect
                def deselect_auto():
                    
                    sensitivity_auto_checkbox.value = False

                @tolerance_spinbox.changed.connect
                def update_flood_mask():
                    
                    current_layer = global_state.current_layer
                    if current_layer is None:
                        napari.utils.notifications.show_error(
                            "Please select a seed point. Note, you may need to"
                            " re-compute the artifact mask."
                        )
                        return
                    
                    abx_channel = current_layer.metadata['abx_channel']
                    artifacts = global_state.artifacts
                    pt_id = current_layer.current_properties['id'][0]
                    df = artifacts[abx_channel].seed_layer.features
                    pt_idx = df[df['id'] == pt_id]['tol'].index.to_list()[0]
                    old_tol = df.loc[pt_idx, 'tol']
                    im_transformed = artifacts[abx_channel].transformed
                    old_fill = flood(
                        im_transformed, seed_point=tuple(artifacts[abx_channel].seeds[pt_id]), 
                        tolerance=old_tol
                    ) 
                    new_tol = tolerance_spinbox.value

                    new_fill = flood(
                        im_transformed, seed_point=tuple(artifacts[abx_channel].seeds[pt_id]), 
                        tolerance=new_tol
                    )
                    artifacts[abx_channel].update_mask(
                        artifacts[abx_channel].mask + new_fill - old_fill
                    )
                    artifacts[abx_channel].tols[pt_idx] = new_tol
                    artifacts[abx_channel].seed_layer.features.loc[pt_idx, 'tol'] = new_tol
                
                ################################################################
                
                widget_lst = [widget_combo_1, widget_combo_2]
                widget_names = ['Automated Artifact Detection', 'Fine-tune Seed Points']

                for artifact_info in global_state.artifacts.values():
                    artifact_info.bind_listener_seeds(viewer, global_state, tolerance_spinbox)

                return widget_lst, widget_names

            widgets = [next_sample, arbitrary_sample]
            widget_names = [f'Sample {sample}', 'Sample Selector']
            if self.autoArtifactDetection:
                if self.artifactDetectionMethod == 'MLP':
                    widgets.append(label_artifacts_MLP)  # better function name?
                    widget_names.append(['Automation Module'])
                elif self.artifactDetectionMethod == 'classical':
                    new_widget_lst, new_widget_names = generate_widgets_for_classical_method()
                    widgets.extend(new_widget_lst)  # need to implement this
                    widget_names.extend(new_widget_names)
            napari_widgets = []
            for widget, widget_name in zip(widgets, widget_names):
                napari_widgets.append(
                    viewer.window.add_dock_widget(widget=widget, add_vertical_stretch=False,
                                                  name=widget_name)
                )
            next_sample.widgets.bind(napari_widgets)
            arbitrary_sample.widgets.bind(napari_widgets)

            for napari_widget in napari_widgets:
                napari_widget.widget().setSizePolicy(
                    QtWidgets.QSizePolicy.Preferred,
                    QtWidgets.QSizePolicy.Fixed,
                )

        def clear_widgets(widgets):
            
            for widget in widgets:
                viewer.window.remove_dock_widget(widget)

        def artifact_detection_model_MLP(data):
            
            if global_state.base_clf is None:
                model_path = os.path.join(
                    Path(__file__).absolute().parent,
                    '../pretrained_models/pretrained_model.pkl'
                )
                with open(model_path, 'rb') as f:
                    global_state.base_clf = pickle.load(f)

            # may consider putting this into the sklearn pipeline?!
            # make sure to strip off cell_id first
            if ('CellID' in data.columns):
                model_input = data.loc[:, data.columns != 'CellID']

            # keep just relevant columns for model
            marker_list = [
                'Hoechst0', 'anti_CD3', 'anti_CD45RO', 'Keratin_570', 'aSMA_660', 'CD4_488',
                'CD45_PE', 'PD1_647', 'CD20_488', 'CD68_555', 'CD8a_660', 'CD163_488',
                'FOXP3_570', 'PDL1_647', 'Ecad_488', 'Vimentin_555', 'CDX2_647', 'LaminABC_488',
                'Desmin_555', 'CD31_647', 'PCNA_488', 'CollagenIV_647', 'Area', 'MajorAxisLength',
                'MinorAxisLength', 'Eccentricity', 'Solidity', 'Extent', 'Orientation'
            ]

            model_input = model_input[marker_list]

            # clf_probs = base_clf.predict_proba(data.values)
            clf_preds = global_state.base_clf.predict(model_input)
            clf_proba = global_state.base_clf.predict_proba(model_input)

            return clf_preds, clf_proba

        ###################################################################
        
        extra_layers, layer_type, layer_name, varname_filename_lst = load_extra_layers()
        add_layers(sample)
        add_widgets(global_state.last_sample, sample)

        napari.run()  # blocks until window is closed

        ###################################################################

        idxs_to_drop = {}
        samples = self.samplesForROISelection
        for sample in samples:

            try:
                extra_layers['ROI'][sample]

                logger.info(f'Generating ROI mask(s) for sample: {sample}')

                sample_data = data[['X_centroid', 'Y_centroid', 'CellID']][
                    data['Sample'] == sample].astype(int)

                sample_data['tuple'] = list(
                    zip(sample_data['X_centroid'],
                        sample_data['Y_centroid'])
                )

                file_path = get_filepath(self, check, sample, 'TIF')
                channel_number = marker_channel_number(self, markers, self.counterstainChannel)
                dna, min, max = single_channel_pyramid(file_path, channel=channel_number)

                # create pillow image to convert into boolean mask
                def shape_layer_to_mask(layer_data):
                    
                    img = Image.new(
                        'L', (dna[0].shape[1], dna[0].shape[0]))

                    for shape_type, verts in layer_data:

                        # snap any floating point verts to array
                        selection_verts = np.round(verts).astype(int)  

                        if shape_type == 'ellipse':

                            vertices, triangles = triangulate_ellipse(
                                selection_verts
                            )

                            # flip 2-tuple coordinates returned by
                            # triangulate_ellipse() to draw image mask
                            vertices = [tuple(reversed(tuple(i)))
                                        for i in vertices]

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
                            ImageDraw.Draw(img).polygon(
                                vertices, outline=1, fill=1)

                    # convert pillow image into boolean numpy array
                    mask = np.array(img, dtype=bool)
                    
                    return mask

                # use numpy fancy indexing to get centroids
                # where boolean mask is True
                xs, ys = zip(*sample_data['tuple'])
                ROI_mask = shape_layer_to_mask(extra_layers['ROI'][sample])
                inter1 = ROI_mask[ys, xs]
                sample_data['inter1'] = inter1
            
            except KeyError:
                logger.info(f'No ROIs selected for sample: {sample}')
                sample_data = data[['X_centroid', 'Y_centroid', 'CellID']][
                    data['Sample'] == sample].astype(int)
                sample_data['tuple'] = list(
                    zip(sample_data['X_centroid'],
                        sample_data['Y_centroid'])
                )
                xs, ys = zip(*sample_data['tuple'])
                sample_data['inter1'] = False

            # may want to keep track if auto artifact detection is used separately in global state
            if self.artifactDetectionMethod == 'MLP':
                autoArtifactDetectionUsed = len(global_state.artifact_proba) > 0
            elif self.artifactDetectionMethod == 'classical':
                autoArtifactDetectionUsed = len(global_state.artifacts) > 0
                
            if self.autoArtifactDetection:
                if self.artifactDetectionMethod == 'MLP':
                    autoArtifactDetectionUsed = len(global_state.artifact_proba) > 0
                elif self.artifactDetectionMethod == 'classical':
                    autoArtifactDetectionUsed = sum(
                        [len(d) for d in list(extra_layers.values())[1:]]
                    ) > 0

                if autoArtifactDetectionUsed:
                    if self.artifactDetectionMethod == 'MLP':
                        ROI2_mask = shape_layer_to_mask(extra_layers['ROI2'][sample])
                        inter2 = ~ROI2_mask[ys, xs] & (global_state.artifact_mask != 1) & \
                            (global_state.artifact_proba > 
                             global_state.artifact_detection_threshold)
                    elif self.artifactDetectionMethod == 'classical':
                        masks = []
                        for abx_channel in abx_channels:
                            artifact_info = extra_layers[f'{abx_channel}_mask'].get(sample)
                            if artifact_info is None:
                                continue
                            else:
                                upscaled_mask = upscale(
                                    artifact_info.mask > 0, 
                                    global_state.loaded_ims[abx_channel][0]
                                )
                                masks.append(upscaled_mask)

                        if len(masks) == 0:
                            inter2 = False
                        else:
                            ROI2_mask = np.logical_or.reduce(masks)
                            inter2 = ROI2_mask[ys, xs]

                    sample_data['inter2'] = inter2
            else:
                logger.info(f'Artifact Detection not enabled for sample: {sample}')

            drop_artifact_ids = (
                sample_data['inter2'] if self.autoArtifactDetection and 
                autoArtifactDetectionUsed else False
            )
            
            if self.delintMode is False:
                sample_data['inter1'] = ~sample_data['inter1']

            idxs_to_drop[sample] = list(
                sample_data['CellID'][sample_data['inter1'] | drop_artifact_ids]
            )
       
        print()

        # drop cells from samples
        for sample, cell_ids in idxs_to_drop.items():
            if cell_ids:
                logger.info(f'Dropping cells from sample: {sample}')
                global_idxs_to_drop = data[
                    (data['Sample'] == sample) &
                    (data['CellID'].isin(set(cell_ids)))].index
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
            channel_number = marker_channel_number(self, markers, self.counterstainChannel)
            dna, min, max = single_channel_pyramid(file_path, channel=channel_number)

            fig, ax = plt.subplots()
            ax.imshow(dna[0], cmap='gray')
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
