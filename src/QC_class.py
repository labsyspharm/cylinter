import os
import re
import yaml
import math
import pickle
import subprocess

import gc
import zarr
import hdbscan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# from skimage import io
from skimage.color import gray2rgb
from skimage.filters import gaussian
from skimage.util.dtype import img_as_float
from skimage.util.dtype import img_as_uint

from matplotlib.lines import Line2D
from matplotlib.widgets import Slider, Button
from matplotlib.widgets import TextBox
from matplotlib.colors import ListedColormap
from matplotlib.path import Path

import napari
from tifffile import imread

from sklearn.preprocessing import MinMaxScaler
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize as norm

from natsort import natsorted, order_by_index, index_natsorted
from numcodecs import Blosc
from datetime import datetime
from joblib import Memory
from scipy.stats import ttest_ind
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from decimal import Decimal
from bridson import poisson_disc_samples

from QC_utils import (save_dataframe, read_dataframe, read_markers,
                      loadZarrs, cluster_expression, categorical_cmap,
                      SelectFromCollection)

# map matplotlib color codes to the default seaborn palette
sns.set()
sns.set_color_codes()
_ = plt.plot([0, 1], color='r')
sns.set_color_codes()
_ = plt.plot([0, 2], color='b')
sns.set_color_codes()
_ = plt.plot([0, 3], color='g')
sns.set_color_codes()
_ = plt.plot([0, 4], color='m')
sns.set_color_codes()
_ = plt.plot([0, 5], color='y')
plt.close('all')

SUPPORTED_EXTENSIONS = ['.csv']


def dataset_files(root):
    """
    Returns a list of all supported extension
    files in the specified directory
    """
    total = os.listdir(root)
    pruned = [i for i in total if i.endswith(tuple(SUPPORTED_EXTENSIONS))]
    return pruned


class QC(object):
    def __init__(self,

                 inDir=None,
                 outDir=None,
                 markers_filepath=None,
                 mask_object=None,
                 sample_metadata=None,
                 randomSampleSize=None,

                 numPCAComponents=2,
                 pointSize=125.0,
                 normalize=True,
                 labelPoints=True,
                 condHueDict={
                     'cd': (0.5, 0.5, 0.5, 1.0),
                     'hfd': (1.0, 1.0, 0.0, 1.0)
                     },
                 numTSNEComponents=2,
                 perplexity=50.0,
                 earlyExaggeration=12.0,
                 learningRate=200.0,
                 metric='euclidean',
                 random_state=5,

                 denominator_cluster=2,
                 FDRCorrection=False,

                 bonferroniCorrection=False,

                 numFingernails=10,

                 cropDict={
                     'cd13': ('top', 10000),
                     'hfd3': ('bottom', 11000),
                     'hfd4': ('top', 8000),
                     'hfd8': ('bottom', 7500),
                     'hfd9': ('top', 9500),
                     'hfd11': ('bottom', 6600),
                     'hfd12': ('top', 9000),
                     'hfd13': ('top', 7000),
                     },
                 spatialDict1={
                    'aco2': 0.07, 'glut1': 0.25
                    },
                 spatialDict2={
                    'TCF1': 0, 'CD8T': 1
                    },
                 radiusRange=[40, 600],

                 dfSaveCount=1
                 ):
        """
        Args:
          config.yaml —
            input_dir: path to csv tables and t-CyCIF images (Zarr arrays)
            output_dir: path to save directory
            markers: probes to be used in the analysis
            samples: map of sample names to experimental condition
            replicates: map of sample names to replicate number
            randomSampleSize: analyze a random subset of data; float (0-1)

          performPCA module —
            numPCAComponents: number of PCs
            pointSize: scatter point size
            normalize: scale input vectors individually to
            unit norm (vector length)
            labelPoints: annotate scatter points
            condHueDict: color scatter points by experimental condition

          performTSNE module —
            numTSNEComponents: dimension of the TSNE embedding
            perplexity: related to the number of nearest neighbors
            used in other manifold learning algorithms. Larger datasets
            usually require larger perplexity. Different values can
            result in significanlty different results.
            earlyExaggeration: for larger values, the space between
            natural clusters will be larger in the embedded space.
            learningRate: t-SNE learning rate: usually [10.0, 1000.0]
            metric: string allowed by scipy.spatial.distance.pdist for
            its metric parameter, or a metric listed in
            pairwise.PAIRWISE_DISTANCE_FUNCTIONS.
            random_state: integer, determines the random number generator.
            For reproducible results across multiple function calls.

          frequencyStats —
            denominator_cluster: HDBSCAN cluster to use as the
            denominator when computing cell frequency ratios
            FDRCorrection: true, report p-vals corrected for
            multiple comparisions by False Discovery Rate method (q-vals);
            false, report uncorrected p-vals

          clusterBoxplots —
            bonferroniCorrection: true, report p-vals corrected for
            multiple comparisions by Bonferroni method (q-vals);
            false, report uncorrected p-vals

          curateFingernails —
            numFingernails: number of example cells to randomly draw from
            each HDBSCAN cluster

          spatialAnalysis —
            cropDict: vertical crop coordinate (numpy row) and
            sub-image to use for t-CyCIF images containing more
            than one tissue section
            spatialDict1: cutoff for pixel-level protein signal instensities
            spatialDict2: map of cell state call to HDBSCAN cluster
            for cell states of interest
            radiusRange: range of radii (in pixels) for Poisson-disc sampling

            dfSaveCount: integer (typically 1) from which start counting
            each time the dataframe is updated and saved to disc
        """

        # assert(SOMETHING)  # placeholder for now

        self.inDir = inDir
        self.outDir = outDir
        self.markers_filepath = markers_filepath
        self.mask_object = mask_object
        self.sample_metadata = sample_metadata
        self.randomSampleSize = randomSampleSize

        self.numPCAComponents = numPCAComponents
        self.pointSize = pointSize
        self.normalize = normalize
        self.labelPoints = labelPoints
        self.condHueDict = condHueDict

        self.numTSNEComponents = numTSNEComponents
        self.perplexity = perplexity
        self.earlyExaggeration = earlyExaggeration
        self.learningRate = learningRate
        self.metric = metric
        self.random_stats = random_state

        self.denominator_cluster = denominator_cluster
        self.FDRCorrection = FDRCorrection

        self.bonferroniCorrection = bonferroniCorrection

        self.numFingernails = numFingernails

        self.cropDict = cropDict
        self.spatialDict1 = spatialDict1
        self.spatialDict2 = spatialDict2
        self.radiusRange = radiusRange

        self.dfSaveCount = dfSaveCount

    def getSingleCellData(self, args):

        files = dataset_files(f'{self.inDir}/csv')

        sample_metadata = self.sample_metadata

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        df_list = []
        for file in files:

            sample = file.split('.')[0]

            print(f'Importing single-cell data for sample {sample}.')

            data = pd.read_csv(os.path.join(f'{self.inDir}/csv', file))

            data['Sample'] = sample

            # append dataframe to list
            df_list.append(data)
            del data

        # stack dataframes row-wise
        df = pd.concat(df_list, axis=0)
        del df_list

        # assign global index
        df.reset_index(drop=True, inplace=True)

        # add condition column
        df['Condition'] = [
            sample_metadata[s][0] for s in df['Sample']]

        # add replicate column
        df['Replicate'] = [
            sample_metadata[s][1] for s in df['Sample']]

        # organize columns
        cols = (
            ['CellID', 'Sample', 'Condition', 'Replicate', 'X_centroid',
             'Y_centroid', 'column_centroid', 'row_centroid', 'Area',
             'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Solidity',
             'Extent', 'Orientation'] +
            [f'{i}_cellMask' for i in markers['marker_name']]
             )
        df = df[cols]

        # handle data subsetting
        df = df.sample(frac=self.randomSampleSize, random_state=1)
        df.sort_values(by=['Sample', 'CellID'], inplace=True)
        df.reset_index(drop=True, inplace=True)

        save_dataframe(
            df=df, outDir=self.outDir, moduleName='getSingleCellData')
        self.data = df

        return self.data

    def lassoROIs(self, args):

        df = self.data

        if not os.path.exists(os.path.join(self.outDir, 'lasso_dict.pkl')):

            lasso_dict = {}

            for sample_name, group in df.groupby('Sample'):

                sample_data = group[['X_centroid', 'Y_centroid']].astype(int)

                if '-' in sample_name:
                    img_name = sample_name.split('-')[1]
                else:
                    img_name = sample_name
                dna = imread(
                    f'{self.inDir}/images/{img_name}.tif',
                    key=0  # key 0 is first dna cycle in markers.csv file
                    )

                columns, rows = np.meshgrid(
                    np.arange(dna.shape[1]),
                    np.arange(dna.shape[0])
                    )
                columns, rows = columns.flatten(), rows.flatten()
                pixel_coords = np.vstack((rows, columns)).T

                with napari.gui_qt():
                    viewer = napari.view_image(
                        dna, rgb=False, name=f'{sample_name}_dna'
                        )

                    selection_layer = viewer.add_shapes(
                        shape_type='polygon',
                        face_color=[1.0, 1.0, 1.0, 0.2],
                        edge_color=[0.0, 0.66, 1.0, 1.0],
                        edge_width=10.0,
                        name='ROI(s)'
                        )

                    # @viewer.mouse_drag_callbacks.append
                    # def get_cell_indices(viewer, event):
                    #
                    #     # on mouse press
                    #     yield
                    #
                    #     # on mouse move
                    #     while event.type == 'mouse_move':
                    #         yield
                    #
                    #     # on mouse release
                    #     selection_layer = viewer.layers['ROI_1']
                    #     yield

                sample_data['tuple'] = list(
                    zip(sample_data['Y_centroid'], sample_data['X_centroid'])
                    )

                idxs_to_drop = []
                mask_coords = set()
                cell_coords = set(
                    [tuple(i) for i in np.array(
                        sample_data[['Y_centroid', 'X_centroid']])]
                    )

                # take all cells if no ROIs are drawn
                if len(selection_layer.data) == 0:
                    continue

                else:
                    for roi in selection_layer.data:

                        selection_verts = np.round(roi).astype(int)
                        polygon = Path(selection_verts)

                        grid = polygon.contains_points(pixel_coords)
                        mask = grid.reshape(
                            dna.shape[0], dna.shape[1])

                        mask_coords.update(
                            [tuple(i) for i in np.argwhere(mask)]
                            )

                    inter = mask_coords.intersection(cell_coords)

                    idxs_to_drop.extend(
                        [i[0] for i in sample_data.iterrows() if
                         i[1]['tuple'] not in inter]
                         )

                    lasso_dict[sample_name] = idxs_to_drop

            os.chdir(self.outDir)
            f = open('lasso_dict.pkl', 'wb')
            pickle.dump(lasso_dict, f)
            f.close()

            for key, value in lasso_dict.items():
                df.drop(lasso_dict[key], inplace=True)

            save_dataframe(
                df=df, outDir=self.outDir,  moduleName='lassoROIs')
            self.data = df

            return self.data

        else:

            os.chdir(self.outDir)
            pickle_in = open('lasso_dict.pkl', 'rb')
            lasso_dict = pickle.load(pickle_in)

            for key, value in lasso_dict.items():
                df.drop(lasso_dict[key], inplace=True)

            save_dataframe(df=df, outDir=self.outDir, moduleName='lassoROIs')
            self.data = df

            return self.data

    def dnaIntensityCutoff(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        bins = 100

        histtype = 'stepfilled'

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        n, bins, patches = plt.hist(
            df[f'{dna1}_cellMask'], bins=bins,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype=histtype,
            range=None, label='before'
            )

        plt.title('median DNA intensity')
        plt.ylabel('count')

        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        rnge = [0, bins.max()]

        sLower = Slider(
            axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sLower.label.set_color('b')

        sUpper = Slider(
            axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sUpper.label.set_color('r')

        def update(val):
            [i.remove() for i in ax.get_lines()]
            lowerCutoff = sLower.val
            upperCutoff = sUpper.val
            blueLine = ax.axvline(
                x=lowerCutoff, c='b', linewidth=2.5)
            redLine = ax.axvline(
                x=upperCutoff, c='r', linewidth=2.5)
            return lowerCutoff, upperCutoff

        sLower.on_changed(update)
        sUpper.on_changed(update)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sLower.reset()
            sUpper.reset()
        button.on_clicked(reset)

        def submit(text, df):
            lowerCutoff, upperCutoff = update(val=None)

            # apply lower and upper cutoffs
            df_test = df[
                (df[f'{dna1}_cellMask'] > lowerCutoff) &
                (df[f'{dna1}_cellMask'] < upperCutoff)
                ]

            if text in df['Sample'].unique():

                if '-' in text:
                    img_name = text.split('-')[1]
                else:
                    img_name = text
                dna = imread(
                    f'{self.inDir}/images/{img_name}.tif',
                    key=0
                    )

                fig, ax = plt.subplots()
                ax.imshow(dna, cmap='gray')
                ax.grid(False)
                coords = df_test[
                    ['X_centroid', 'Y_centroid', f'{dna1}_cellMask']][
                        df_test['Sample'] == text]
                sp = ax.scatter(
                    coords['X_centroid'], coords['Y_centroid'], s=3.5,
                    c=coords[f'{dna1}_cellMask'], cmap='viridis')
                plt.title(
                    f'Sample {text}. '
                    f'Selected cells colored by DNA intensity.'
                    )
                plt.colorbar(sp)
                plt.show(block=True)

        axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
        text_box = TextBox(
            axbox, 'evaluation sample name', initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.on_submit(lambda val: submit(val, df))
        plt.show(block=True)

        lowerCutoff, upperCutoff = update(val=None)

        fig, ax = plt.subplots()
        # plot DNA intensity histogram BEFORE filtering
        plt.hist(
            df[f'{dna1}_cellMask'], bins=bins,
            density=False, color='b', ec='none',
            alpha=0.5, histtype=histtype,
            range=None, label='before'
            )

        # apply lower and upper cutoffs
        df = df[
            (df[f'{dna1}_cellMask'] > lowerCutoff) &
            (df[f'{dna1}_cellMask'] < upperCutoff)
            ]

        # plot DNA intensity histogram AFTER filtering
        plt.hist(
            df[f'{dna1}_cellMask'], bins=bins, color='r', ec='none',
            alpha=0.5, histtype=histtype, range=None, label='after')
        plt.xlabel('median DNA intensity')
        plt.ylabel('count')

        legend_elements = []
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='excluded data',
                   markerfacecolor='b', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='included data',
                   markerfacecolor='r', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        plt.legend(
            handles=legend_elements, prop={'size': 10},
            loc='best'
            )

        plt.tight_layout()
        plt.savefig(os.path.join(self.outDir, 'histogram_dna.pdf'))
        plt.close()

        save_dataframe(
            df=df, outDir=self.outDir, moduleName='dnaIntensityCutoff'
            )
        self.data = df

        return self.data

    def nuclearAreaCutoff(self, args):

        df = self.data

        bins = 100

        histtype = 'stepfilled'

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        n, bins, patches = plt.hist(
            df['Area'], bins=bins,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype=histtype,
            range=None, label='before'
            )

        plt.title('nuclear area')
        plt.ylabel('count')

        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

        rnge = [0, bins.max()]

        sLower = Slider(
            axLowerCutoff, 'lowerCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sLower.label.set_color('b')

        sUpper = Slider(
            axUpperCutoff, 'upperCutoff', rnge[0], rnge[1],
            valinit=0.00, valstep=(rnge[1]/100000)
            )
        sUpper.label.set_color('r')

        def update(val):
            [i.remove() for i in ax.get_lines()]
            lowerCutoff = sLower.val
            upperCutoff = sUpper.val
            blueLine = ax.axvline(
                x=lowerCutoff, c='b', linewidth=2.5)
            redLine = ax.axvline(
                x=upperCutoff, c='r', linewidth=2.5)
            return lowerCutoff, upperCutoff

        sLower.on_changed(update)
        sUpper.on_changed(update)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

        def reset(event):
            sLower.reset()
            sUpper.reset()
        button.on_clicked(reset)

        def submit(text, df):
            lowerCutoff, upperCutoff = update(val=None)

            # apply lower and upper cutoffs
            df_test = df[
                (df['Area'] > lowerCutoff) &
                (df['Area'] < upperCutoff)
                ]

            if text in df['Sample'].unique():

                if '-' in text:
                    img_name = text.split('-')[1]
                else:
                    img_name = text
                dna = imread(
                    f'{self.inDir}/images/{img_name}.tif',
                    key=0
                    )

                fig, ax = plt.subplots()
                ax.imshow(dna, cmap='gray')
                ax.grid(False)
                coords = df_test[
                    ['X_centroid', 'Y_centroid', 'Area']][
                        df_test['Sample'] == text]
                sp = ax.scatter(
                    coords['X_centroid'], coords['Y_centroid'], s=3.5,
                    c=coords['Area'], cmap='viridis')
                plt.title(
                    f'Sample {text}. '
                    f'Selected cells colored by nuclear area.'
                    )
                plt.colorbar(sp)
                plt.show(block=True)

        axbox = plt.axes([0.4, 0.025, 0.35, 0.045])
        text_box = TextBox(
            axbox, 'evaluation sample name', initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.on_submit(lambda val: submit(val, df))
        plt.show(block=True)

        lowerCutoff, upperCutoff = update(val=None)

        fig, ax = plt.subplots()
        # plot nuclear area histogram BEFORE filtering
        plt.hist(
            df['Area'], bins=bins,
            density=False, color='b', ec='none',
            alpha=0.5, histtype=histtype,
            range=None, label='before'
            )

        # apply lower and upper cutoffs
        df = df[
            (df['Area'] > lowerCutoff) &
            (df['Area'] < upperCutoff)
            ]

        # plot nuclear area histogram AFTER filtering
        plt.hist(
            df['Area'], bins=bins, color='r', ec='none',
            alpha=0.5, histtype=histtype, range=None, label='after')
        plt.xlabel('nuclear area')
        plt.ylabel('count')

        legend_elements = []
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='excluded data',
                   markerfacecolor='b', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        legend_elements.append(
            Line2D([0], [0], marker='o', color='none',
                   label='included data',
                   markerfacecolor='r', alpha=0.5,
                   markeredgecolor='none', lw=0.001,
                   markersize=8)
                   )
        plt.legend(
            handles=legend_elements, prop={'size': 10},
            loc='best'
            )

        plt.tight_layout()
        plt.savefig(os.path.join(self.outDir, 'histogram_area.pdf'))
        plt.close()

        # save images
        lasso_dir = os.path.join(self.outDir, 'lassos')
        if not os.path.exists(lasso_dir):
            os.mkdir(lasso_dir)

        for sample_name in df['Sample'].unique():

            if '-' in sample_name:
                img_name = sample_name.split('-')[1]
            else:
                img_name = sample_name
            dna = imread(
                f'{self.inDir}/images/{img_name}.tif',
                key=0
                )

            fig, ax = plt.subplots()
            ax.imshow(dna, cmap='gray')
            ax.grid(False)
            coords = df[['X_centroid', 'Y_centroid', 'Area']][
                df['Sample'] == sample_name]
            sp = ax.scatter(
                coords['X_centroid'], coords['Y_centroid'],
                s=coords['Area']/100, c=coords['Area'], cmap='viridis'
                )
            plt.title(
                f'Sample {sample_name}. Lasso colored by nuclear area')
            plt.colorbar(sp)
            plt.savefig(
                os.path.join(
                    lasso_dir, f'{sample_name}.png'), dpi=1000)
            plt.close('all')

        save_dataframe(
            df=df, outDir=self.outDir, moduleName='nuclearAreaCutoff'
            )
        self.data = df

        return self.data

    def crossCyleCorrelation(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        ratios = pd.DataFrame(
            [np.log10(
                (df[f'{dna1}_cellMask'] + 0.00001) /
                (df[i] + 0.00001)) for i in
                df.columns[df.columns.str.contains(dna_prefix)]]).T

        list1 = [i for i in ratios.columns if i.startswith('Unnamed')]
        list2 = [f'1/{i+1}' for i in range(1, len(list1)+1)]
        ratio_columns = dict(zip(list1, list2))
        ratios.rename(columns=ratio_columns, inplace=True)
        ratios.drop(f'{dna1}_cellMask', axis=1, inplace=True)
        ratios['sample'] = df['Sample']

        ratios = ratios.melt(id_vars=['sample'])
        ratios.rename(columns={'variable': 'cycles'}, inplace=True)
        ratios.rename(columns={'value': 'log10(ratio)'}, inplace=True)

        sns.set(font_scale=0.5)
        sns.set_style('whitegrid')

        g = sns.FacetGrid(
            ratios, row='sample',
            col='cycles', sharey=False
            )
        g = g.map(
            plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
            ec='none', range=(-1, 1), bins=200, density=True
            )

        plt.savefig(
            os.path.join(self.outDir, 'cycle_correlation(logRatio).pdf'))
        plt.close()

        subprocess.call(
            ['open', '-a', 'Preview', os.path.join(
                self.outDir, 'cycle_correlation(logRatio).pdf')])

        def submit(text, g):

            count_cutoff = float(text)

            sns.set(font_scale=0.5)
            sns.set_style('whitegrid')

            g = sns.FacetGrid(
                ratios, row='sample',
                col='cycles', sharey=False
                )
            g = g.map(
                plt.hist, 'log10(ratio)', color='r', histtype='stepfilled',
                ec='none', range=(-1, 1), bins=200, density=True
                )

            for ax in g.axes.ravel():
                ax.axhline(y=count_cutoff, c='k', linewidth=0.5)

            plt.savefig(
                os.path.join(self.outDir, 'cycle_correlation(logRatio).pdf'))

            subprocess.call(
                ['open', '-a', 'Preview',
                 os.path.join(self.outDir, 'cycle_correlation(logRatio).pdf')]
                 )

            os.chdir(self.outDir)
            f = open('count_cutoff.pkl', 'wb')
            pickle.dump(count_cutoff, f)
            f.close()

            plt.show(block=False)
            plt.close()

        plt.rcParams['figure.figsize'] = (6, 3)
        axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
        text_box = TextBox(
            axbox, 'countCutoff', initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.label.set_size(15)

        text_box.on_submit(lambda val: submit(val, g))

        plt.show(block=True)

        os.chdir(self.outDir)
        pickle_in = open('count_cutoff.pkl', 'rb')
        count_cutoff = pickle.load(pickle_in)

        # grab dna and sample columns of dataframe
        facet_input = df.loc[
            :, df.columns.str.contains(f'{dna_prefix}|Sample')
            ]

        # initialize a set to append indices to drop
        indices_to_drop = set()

        # loop over samples
        for sample in df['Sample'].unique():

            # slice sample-specific data
            sample_df = df[df['Sample'] == sample]

            # loop over dna cycle columns
            for col_name in sample_df.columns:
                if dna_prefix in col_name:

                    # get cycle number
                    cycle_num = str(re.search(r'\d+', col_name).group())

                    if cycle_num != '1':

                        # get ratios
                        ratios = np.log10(
                            (sample_df[f'{dna1}_cellMask'] + 0.00001) /
                            (sample_df[f'{dna_prefix}{cycle_num}_cellMask'] +
                             0.00001))

                        # get histogram elements
                        sns.set_style('whitegrid')
                        fig, ax = plt.subplots()
                        counts, bins, patches = plt.hist(
                            ratios, color='r', histtype='stepfilled',
                            ec='none', range=(-1, 1),
                            bins=200, density=True
                            )
                        plt.close('all')

                        # plot histogram of log(ratios)
                        # with a horizontal line at cutoff point
                        # ax.axhline(y=count_cutoff, c='k', linewidth=0.5)
                        # plt.title(sample + '_' + col_name)
                        # plt.show(block=True)

                        # get bin values (i.e. ratios) where cell
                        # counts are > than count_cutoff
                        count_indices = np.where(counts > count_cutoff)
                        bin_values = [bins[i] for i in count_indices[0]]
                        min_bin_val = min(bin_values)
                        max_bin_val = max(bin_values)

                        # get indices in log(ratio) series outside
                        # min_bin_val and max_bin_val
                        idx = list(
                            ratios.index[
                                (ratios < min_bin_val) |
                                (ratios > max_bin_val)]
                                )

                        # append indices of uncorrelated
                        # log(ratios) to idx_list
                        indices_to_drop |= set(idx)

        # filter dataframe by selecting indices NOT in the indices_to_drop list
        df = df.loc[~df.index.isin(indices_to_drop)]

        # grab dna and sample columns
        facet_input = df.loc[
            :, df.columns.str.contains(f'{dna_prefix}|Sample')].copy()

        # plot cell dropout facet (per cycle)
        facet_per_cycle_melt = (
            facet_input.drop(['Sample'], axis=1)
            .melt(var_name='cycle', value_name='signal'))

        fig, ax = plt.subplots(figsize=(5, 5))
        g = sns.FacetGrid(facet_per_cycle_melt, col='cycle', col_wrap=4)

        g.map(
            lambda y, color: plt.scatter(
                y, facet_input[f'{dna1}_cellMask'],
                s=0.25, alpha=0.1, linewidth=None,
                marker='o', c='r'), 'signal')
        plt.savefig(
            os.path.join(
                self.outDir, 'cycle_correlation(perCycle).pdf'))
        plt.close('all')

        # plot cell dropout facet (per sample per cycle)
        # take a fraction of the total dataframe for plotting
        facet_per_sample_per_cycle_melt = (
            facet_input.sample(frac=0.1)
            .melt(id_vars=['Sample'],
                  var_name='cycle', value_name='signal')
            )

        facet_per_sample_per_cycle_melt['Sample'] = pd.Categorical(
            facet_per_sample_per_cycle_melt['Sample'], ordered=True,
            categories=natsorted(
                facet_per_sample_per_cycle_melt['Sample'].unique()))

        facet_per_sample_per_cycle_melt['cycle'] = pd.Categorical(
            facet_per_sample_per_cycle_melt['cycle'], ordered=True,
            categories=natsorted(
                facet_per_sample_per_cycle_melt['cycle'].unique()))

        facet_per_sample_per_cycle_melt = (
            facet_per_sample_per_cycle_melt.sort_values(['Sample', 'cycle'])
            )

        fig, ax = plt.subplots(figsize=(5, 5))

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(
                facet_per_sample_per_cycle_melt['Sample'].unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False
            )

        sample_color_dict = dict(
            zip(
                natsorted(facet_per_sample_per_cycle_melt['Sample'].unique()),
                cmap.colors)
                )

        g = sns.FacetGrid(
            facet_per_sample_per_cycle_melt, col='cycle', hue='Sample'
            )

        g.map(
            lambda sam, y, color, **kwargs: plt.scatter(
                y, facet_per_sample_per_cycle_melt.loc[
                    (facet_per_sample_per_cycle_melt['Sample'] ==
                     sam.unique()[0])
                    & (facet_per_sample_per_cycle_melt['cycle'] ==
                       f'{dna1}_cellMask'), 'signal'],
                c=np.reshape(sample_color_dict[sam.unique()[0]], (-1, 3)),
                s=0.25, linewidth=None, marker='o', **kwargs),
            'Sample', 'signal'
            )

        plt.legend(markerscale=10, bbox_to_anchor=(1.1, 1.05))

        plt.savefig(
            os.path.join(
                self.outDir, 'cycle_correlation(perSample).pdf'),
            bbox_inches='tight')
        plt.close('all')

        save_dataframe(
            df=df, outDir=self.outDir, moduleName='crossCyleCorrelation'
            )
        self.data = df

        return self.data

    def log10transform(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        # rescale abx intensities between 0 and 1
        # min_max_scaler = MinMaxScaler()
        # df[abx_channels] = min_max_scaler.fit_transform(df[abx_channels])

        # log10 transform
        df[abx_channels] += 0.00000000001
        df[abx_channels] = np.log10(df[abx_channels])

        save_dataframe(
            df=df, outDir=self.outDir, moduleName='log10transform'
            )
        self.data = df

        return self.data

    def pruneOutliers(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        sample_metadata = self.sample_metadata

        # save images
        correlation_dir = os.path.join(self.outDir, 'lassos')
        if not os.path.exists(correlation_dir):
            os.mkdir(correlation_dir)

        # plot raw signal intensity histrograms for inspection
        for ab in abx_channels:

            print(ab)

            hist_facet = (
                df[['Sample', 'Condition', 'Area'] + [ab]]
                .sample(frac=1.0)
                .melt(
                    id_vars=['Sample', 'Condition', 'Area'],
                    var_name='channel', value_name='signal')
                )

            # sort and create labels column for plotting
            hist_facet.sort_values(by=['Condition', 'Sample'], inplace=True)
            hist_facet['for_plot'] = (
                hist_facet['Condition'] + '_' +
                hist_facet['Sample']
                )

            sns.set_style('white')
            g = sns.FacetGrid(
                hist_facet, col='for_plot', col_wrap=15,
                height=0.5, aspect=1.0, sharex=True, sharey=False
                )

            g.map(
                lambda x, y, color: plt.hexbin(
                    x, y, gridsize=20, linewidths=0.02, color='dimgrey'),
                'signal', 'Area')

            g.set_titles(
                col_template="{col_name}", fontweight='bold',
                size=1.5, pad=0.75
                )

            g.fig.suptitle(ab, y=0.97)

            for ax in g.axes.flatten():
                ax.tick_params(
                    axis='both', which='major',
                    labelsize=0.2, pad=-5
                    )
                ax.xaxis.label.set_size(2.0)
                ax.yaxis.label.set_size(2.0)
                ax.spines['left'].set_visible(False)
                ax.spines['bottom'].set_visible(False)

            plt.subplots_adjust(
                left=0.01, bottom=0.01, right=0.99,
                top=0.90, hspace=0.4, wspace=0.4
                )

            plt.savefig(
                os.path.join(
                    correlation_dir, f'{ab}_hexbins(raw).pdf'))

            plt.close()

            subprocess.call(
                ['open', '-a', 'Preview', os.path.join(
                    correlation_dir, f'{ab}_hexbins(raw).pdf')])

            def submit(text, df):

                lowerPercentileCutoff = float(text.split(',')[0])
                upperPercentileCutoff = float(text.split(',')[1])

                ab = str(text.split(',')[2][1:])

                # for s in df['Sample'].unique():

                data = df[ab]

                indices_to_drop = []

                # add row index to list if column value < lower bound
                indices_to_drop.extend(
                    data.index[
                        data < np.percentile(data, lowerPercentileCutoff)]
                        )

                # add row index to list if column value > upper bound
                indices_to_drop.extend(
                    data.index[
                        data > np.percentile(data, upperPercentileCutoff)])

                df = df.drop(
                    labels=set(indices_to_drop), axis=0,
                    inplace=False, errors='raise')

                pruned_data = df[ab]

                # rescale channel signal intensities
                scaler = (
                    MinMaxScaler(feature_range=(0, 1), copy=True)
                    .fit(pruned_data.values.reshape(-1, 1))
                        )
                rescaled_data = scaler.transform(
                    pruned_data.values.reshape(-1, 1)
                    )

                rescaled_data = pd.DataFrame(
                    data=rescaled_data,
                    index=pruned_data.index,
                    ).rename(columns={0: ab})

                df.update(rescaled_data)

                # plot pruned and rescaled signal intensity histrograms
                hist_facet = (
                    df[['Sample', 'Condition', 'Area'] + [ab]]
                    .sample(frac=1.0)
                    .melt(
                        id_vars=['Sample', 'Condition', 'Area'],
                        var_name='channel', value_name='signal')
                    )

                # sort and create labels column for plotting
                hist_facet.sort_values(
                    by=['Condition', 'Sample'], inplace=True
                    )
                hist_facet['for_plot'] = (
                    hist_facet['Condition'] + '_' +
                    hist_facet['Sample']
                    )

                g = sns.FacetGrid(
                    hist_facet, col='for_plot', col_wrap=15,
                    height=0.5, aspect=1.0, sharex=True, sharey=False
                    )

                g.map(
                    lambda x, y, color: plt.hexbin(
                        x, y, gridsize=20, linewidths=0.02, color='dimgrey'),
                    'signal', 'Area')

                g.set_titles(
                    col_template="{col_name}", fontweight='bold',
                    size=1.5, pad=0.75
                    )

                g.fig.suptitle(ab, y=0.97)

                for ax in g.axes.flatten():
                    ax.tick_params(
                        axis='both', which='major',
                        labelsize=0.2, pad=-5
                        )
                    ax.xaxis.label.set_size(2.0)
                    ax.yaxis.label.set_size(2.0)
                    ax.spines['left'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)

                plt.subplots_adjust(
                    left=0.01, bottom=0.01, right=0.99,
                    top=0.90, hspace=0.4, wspace=0.4
                    )

                plt.savefig(
                    os.path.join(
                        correlation_dir,
                        f'{ab}_hexbins(pruned_rescaled).pdf')
                        )

                plt.close()

                subprocess.call(
                    ['open', '-a', 'Preview', os.path.join(
                        correlation_dir,
                        f'{ab}_hexbins(pruned_rescaled).pdf')]
                        )
                ##################

            plt.rcParams['figure.figsize'] = (7, 3)
            axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
            text_box = TextBox(
                axbox,
                'lowerCutoff, upperCutoff, abxName',
                initial='',
                color='0.95',
                hovercolor='1.0',
                label_pad=0.05
                )
            text_box.label.set_size(10)

            text_box.on_submit(lambda text: submit(text, df))

            plt.show(block=True)

        save_dataframe(
            df=df, outDir=self.outDir, moduleName='pruneOutliers'
            )
        self.data = df

        return self.data

    def makeZarrs(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        # make directory to store zarr arrays
        zarrs_dir = os.path.join(self.outDir, 'zarrs')
        if not os.path.exists(zarrs_dir):
            os.mkdir(zarrs_dir)

            # initialize a dictionary to index and access zarr arrays
            zs = {}

            # select data compression algorithm for zarr arrays
            compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE)

            # loop over tissue samples
            for sample_name in df['Sample'].unique():

                # read segmentation mask for sample to crop
                # ashlar_output images to the same size
                if '-' in sample_name:
                    img_name = sample_name.split('-')[1]
                else:
                    img_name = sample_name

                segmentation = imread(
                    os.path.join(
                        f'{self.inDir}/masks',
                        f'{img_name}-{self.mask_object}.tif')
                        )

                # read DNA ashlar_output image
                img = imread(
                    f'{self.inDir}/images/{img_name}.tif',
                    key=0  # key 0 is first dna cycle in markers.csv file
                    )

                # crop DNA image to size of U-Net segmentation mask
                img = img[0:segmentation.shape[0], 0:segmentation.shape[1]]

                # initialize zarr array for DNA image, append to a dictionary
                zs[f'{img_name}_{dna1}'] = zarr.open(
                    f'{zarrs_dir}/{img_name}_{dna1}.zarr', mode='w',
                    shape=(
                        img.shape[0], img.shape[1]),
                    chunks=(1000, 1000), dtype='uint16', compressor=compressor
                    )

                # update zarr array with DNA image data
                zs[f'{img_name}_{dna1}'][:] = img

                # update zarr array with all antibody images
                # for k, v in cycle_map.items():
                abx_channels = [
                    i for i in markers['marker_name'] if
                    dna_prefix not in i
                    ]

                for ab in abx_channels:

                    print(f'Storing sample {img_name} {ab}.tif')

                    channel_number = markers['channel_number'][
                        markers['marker_name'] == ab].iloc[0]

                    # read antibody image
                    img = imread(
                        f'{self.inDir}/images/{img_name}.tif',
                        key=(channel_number-1)
                        )

                    # crop image to size of U-Net segmentation mask
                    img = img[0:segmentation.shape[0], 0:segmentation.shape[1]]

                    # initialize zarr array for image, append to a dictionary
                    zs[f'{img_name}_{ab}'] = zarr.open(
                        f'{zarrs_dir}/{img_name}_{ab}.zarr', mode='w',
                        shape=(
                            img.shape[0], img.shape[1]),
                        chunks=(1000, 1000), dtype='uint16',
                        compressor=compressor
                        )

                    # update zarr array with image data
                    zs[f'{img_name}_{ab}'][:] = img
            print()

            # apply OMERO image rendering settings if available and desired
            if os.path.exists(f'{self.inDir}/images/omero_settings.yml'):

                # load settings
                img_rend_settings = yaml.safe_load(
                    open(f'{self.inDir}/images/omero_settings.yml'))

                # set channel start to 0.0 to see all background,
                # or use: img_rend_settings['channels'][channel]['start']
                settings_dict = {}
                for channel in img_rend_settings['channels']:
                    settings_dict[
                        img_rend_settings['channels'][channel]['label']] = ((
                            img_rend_settings['channels'][channel]['start'],
                            img_rend_settings['channels'][channel]['end']),
                            img_rend_settings['channels'][channel]['color'])

                for k, v in zs.items():

                    print(f'Applying OMERO signal intensity cutoffs for {k}')

                    param_map_key = k.split('_')[1]

                    bottom_omero_cutoff = settings_dict[param_map_key][0][0]
                    top_omero_cutoff = settings_dict[param_map_key][0][1]

                    temp = img_as_float(zs[k])
                    temp -= (bottom_omero_cutoff/65535)
                    temp /= (
                        (top_omero_cutoff/65535)-(bottom_omero_cutoff/65535)
                        )
                    temp = np.clip(temp, 0, 1)

                    zs[k][:] = img_as_uint(temp)

        return self.data

    def performPCA(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        sample_metadata = self.sample_metadata

        # print(zs.keys())
        medians = df.groupby(
            ['Sample']).median()[abx_channels]

        # set grid background style
        sns.set_style('whitegrid')

        # specify PCA parameters
        pca = PCA(self.numPCAComponents)

        idx = medians.index

        # normalize signal intensities across samples (axis=0)
        if self.normalize is True:
            medians = norm(
                medians, norm='l2', axis=0, copy=True, return_norm=False
                )
        else:
            medians = medians.values

        # apply PCA parameters to data
        projected = pca.fit_transform(medians)

        # generate dataframe for plot input
        scatter_input = pd.DataFrame(data=projected, index=idx)
        scatter_input.rename(columns={0: 'PC1', 1: 'PC2'}, inplace=True)

        # switch sample names for condition names
        scatter_input.index = [
            sample_metadata[i][0] for i in scatter_input.index
            ]

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(scatter_input.index.unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False
            )

        sample_color_dict = dict(
            zip(
                natsorted(scatter_input.index.unique()),
                cmap.colors)
                )

        # plot scores plot for first 2 PCs
        sns.scatterplot(
            data=scatter_input, x='PC1', y='PC2',
            hue=scatter_input.index,
            palette=sample_color_dict,
            edgecolor='k', s=self.pointSize, alpha=1.0, legend=False)

        # annotate data points
        if self.labelPoints is True:
            for label, x, y in zip(
              scatter_input.index,
              scatter_input['PC1'], scatter_input['PC2']):

                plt.annotate(
                    label, xy=(x, y), xytext=(3, 5), size=4.0,
                    textcoords='offset points', ha='left', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.1', fc='yellow',
                              alpha=0.0))

        plt.xlabel(
            'PC1 ' + '(' + str(
                round((pca.explained_variance_ratio_[0] * 100), 2)) + '%)',
            fontsize=12)
        plt.ylabel(
            'PC2 ' + '(' + str(
                round((pca.explained_variance_ratio_[1] * 100), 2)) + '%)',
            fontsize=12)

        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.tight_layout()
        plt.savefig(os.path.join(self.outDir, 'sample_scores.pdf'))
        plt.close('all')

        return self.data

    def performTSNE(self, args):

        # df = self.data

        df = pd.read_csv(
            '/Users/greg/projects/cycif-qc/output/tma22/dataframe_archive/' +
            'pruneOutliers.csv', index_col=0)

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        if os.path.exists(os.path.join(self.outDir, 'embedding.npy')):
            embedded = np.load(os.path.join(self.outDir, 'embedding.npy'))
            df['emb1'] = embedded[:, 0]
            df['emb2'] = embedded[:, 1]

        else:
            startTime = datetime.now()
            embedded = TSNE(
                n_components=self.numTSNEComponents,
                init='pca',
                perplexity=self.perplexity,
                early_exaggeration=self.earlyExaggeration,
                learning_rate=self.learningRate,
                metric=self.metric,
                random_state=5,
                n_jobs=-1).fit_transform(df[abx_channels])
            print('Embedding completed in ' + str(datetime.now() - startTime))

            np.save(os.path.join(self.outDir, 'embedding'), embedded)
            df['emb1'] = embedded[:, 0]
            df['emb2'] = embedded[:, 1]

        sns.set_style('white')

        def submit(text, df):

            markers, dna1, dna_prefix = read_markers(
                markers_filepath=self.markers_filepath
                )

            abx_channels = [
                f'{i}_cellMask' for i in markers['marker_name'] if
                dna_prefix not in i
                ]

            numerical_input = text.split('.')[0].strip()
            tup = tuple(map(int, numerical_input.split('-')))

            if len(tup) == 1:

                mylist = [tup[0]]

                for i in mylist:

                    min_cluster_size = i

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=min_cluster_size, min_samples=None,
                        metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,
                        memory=Memory(
                            location=None),
                        approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False).fit(
                            df[['emb1', 'emb2']]
                            )
                    df['cluster'] = clustering.labels_

                    print(
                        f'min_cluster_size={i}', np.unique(clustering.labels_)
                        )

                    fig, (ax1, ax2) = plt.subplots(1, 2)
                    fig.suptitle(
                        f'min cluster size = {min_cluster_size}',
                        fontsize=6
                        )

                    plt.subplots_adjust(
                        wspace=0.7,
                        left=0.04
                        )

                    # PLOT TSNE
                    for color_by in ['cluster', 'Condition']:

                        highlight = 'none'

                        if color_by == 'cluster':

                            num_colors_required = len(
                                df[color_by].unique()) - 1
                            catagories = math.ceil(
                                num_colors_required/10
                                )

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(df[color_by].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            # make black the first color to specify
                            # cluster outliers (i.e. cluster -1 cells)
                            # colors_list.insert(0, (0.0, 0.0, 0.0))
                            cmap = ListedColormap(
                                np.insert(
                                    arr=cmap.colors, obj=0,
                                    values=[0, 0, 0], axis=0)
                                    )

                            # trim qualitative cmap to number of unique samples
                            trim = (
                                len(cmap.colors) - len(df[color_by].unique())
                                )
                            cmap = ListedColormap(
                                cmap.colors[:-trim]
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df[color_by].unique()),
                                    list(range(len(df[color_by].unique()))))
                                    )

                            c = [sample_dict[i] for i in df[color_by]]

                            ax1.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=135000/len(df),
                                ec=[
                                    'k' if i == highlight else 'none' for
                                    i in df[color_by]
                                    ],
                                linewidth=0.1
                                )

                            ax1.axis('equal')
                            ax1.tick_params(labelsize=5)
                            ax1.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(df[color_by].unique())
                              ):

                                hi_markers = cluster_expression(
                                    df=df, markers=abx_channels,
                                    cluster=i, num_proteins=3
                                    )

                                legend_elements.append(
                                    Line2D([0], [0], marker='o', color='none',
                                           label=f'Cluster {i} {hi_markers}',
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor='none', lw=0.001,
                                           markersize=4)
                                           )

                            ax1.legend(
                                handles=legend_elements,
                                prop={'size': 4},
                                bbox_to_anchor=[1.02, 1.0]
                                )

                        elif color_by == 'Condition':

                            num_colors_required = len(df[color_by].unique())
                            numSubcatagories = math.ceil(
                                num_colors_required/10
                                )

                            # build cmap
                            cmap = categorical_cmap(
                                numUniqueSamples=len(df[color_by].unique()),
                                numCatagories=10,
                                cmap='tab10',
                                continuous=False
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df[color_by].unique()),
                                    list(range(len(df[color_by].unique()))))
                                    )

                            c = [sample_dict[i] for i in df[color_by]]

                            ax2.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=135000/len(df),
                                ec=[
                                    'k' if i == highlight else 'none' for
                                    i in df[color_by]
                                    ],
                                linewidth=0.1
                                )

                            ax2.axis('equal')
                            ax2.tick_params(labelsize=5)
                            ax2.grid(False)

                            legend_elements = []
                            for e, i in enumerate(
                                natsorted(df[color_by].unique())
                              ):

                                if i == highlight:
                                    markeredgecolor = 'k'
                                else:
                                    markeredgecolor = 'none'

                                legend_elements.append(
                                    Line2D([0], [0], marker='o', color='none',
                                           label=f'Sample {i}',
                                           markerfacecolor=cmap.colors[e],
                                           markeredgecolor=markeredgecolor,
                                           lw=0.001,
                                           markersize=4)
                                           )

                            ax2.legend(
                                handles=legend_elements,
                                prop={'size': 4},
                                bbox_to_anchor=[1.02, 1.0]
                                )

                    plt.tight_layout()

                    if '.save' in text:

                        plt.savefig(
                            os.path.join(
                                self.outDir,
                                f'tsne_{color_by}.pdf')
                            )

                        save_dataframe(
                            df=df, outDir=self.outDir,
                            moduleName='performTSNE'
                            )
                        self.data = df

                        plt.close('all')

                    plt.show(block=False)

            else:

                # df = df.sample(frac=0.01, random_state=22)

                mylist = list(range(tup[0], tup[1] + 1, 1))
                mylist.reverse()  # run higher sizes first for plot order

                for i in mylist:

                    min_cluster_size = i

                    clustering = hdbscan.HDBSCAN(
                        min_cluster_size=min_cluster_size, min_samples=None,
                        metric='euclidean', alpha=1.0, p=None,
                        algorithm='best', leaf_size=40,
                        memory=Memory(
                            location=None),
                        approx_min_span_tree=True,
                        gen_min_span_tree=False, core_dist_n_jobs=4,
                        cluster_selection_method='eom',
                        allow_single_cluster=False,
                        prediction_data=False,
                        match_reference_implementation=False).fit(
                            df[['emb1', 'emb2']]
                            )
                    df['cluster'] = clustering.labels_

                    print(
                        f'min_cluster_size={i}', np.unique(clustering.labels_)
                        )

        plt.rcParams['figure.figsize'] = (7, 3)
        axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
        text_box = TextBox(
            axbox,
            'min_cluster_size (single # or range #-# .save)',
            initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.label.set_size(10)
        text_box.on_submit(lambda text: submit(text, df))
        plt.show(block=True)

        return self.data

    def getClustermap(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        clustermap_input = df[df['cluster'] != -1]

        cluster_heatmap_input = clustermap_input[
            abx_channels + ['cluster']].groupby('cluster').mean()

        sns.set(font_scale=0.4)
        g = sns.clustermap(
            cluster_heatmap_input, cmap='viridis', standard_scale=1,
            square=False, yticklabels=1, linewidth=0.1, cbar=True
            )

        plt.savefig(
            os.path.join(
                self.outDir, 'clustermap.pdf'), bbox_inches='tight')

        plt.show(block=True)

        return self.data

    def lassoClusters(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        subplot_kw = dict(
            xlim=(df['emb1'].min(), df['emb1'].max()),
            ylim=(df['emb2'].min(), df['emb2'].max()),
            autoscale_on=False)

        subplot_kw = dict()

        fig, ax = plt.subplots(subplot_kw=subplot_kw)

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(df['cluster'].unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False
            )
        # add black as first element to represent HDBSCAN outliers
        cmap.colors = np.append([[0.0, 0.0, 0.0]], cmap.colors, axis=0)
        cmap.colors = np.delete(cmap.colors, -1, axis=0)

        pts = ax.scatter(
            df['emb1'],
            df['emb2'],
            c=df['cluster'],
            cmap=cmap,
            s=2.0,
            ec='none'
            )
        # pts = ax.scatter(
        #     lasso_X['emb1'], lasso_X['emb2'], c='lime', s=0.0005
        #     )
        selector = SelectFromCollection(ax, pts)

        def accept(event):
            if event.key == "enter":
                print("Selected points:")
                print(selector.xys[selector.ind])
                selector.disconnect()
                ax.set_title("")
                fig.canvas.draw()

        fig.canvas.mpl_connect("key_press_event", accept)
        ax.set_title('TSNE embedding. Press enter to accept selected points.')
        ax.set_aspect('equal')
        plt.show(block=True)

        # filter dataframe
        drop_idx = (
            set(np.array(list(range(selector.xys.data.shape[0]))))
            - set(selector.ind))
        idx = df.iloc[list(drop_idx)].index

        # show highest expression channels
        markers = df.copy()
        markers.loc[~markers.index.isin(idx), 'cluster'] = 1000
        hi_markers = cluster_expression(
            df=markers, markers=abx_channels, cluster=1000, num_proteins=3
            )
        print(hi_markers)

        return self.data

    def cellDensities(self, args):

        df = self.data

        sample_metadata = self.sample_metadata

        facet_input = df[df['cluster'] >= 0]

        facet_input = facet_input.groupby(
            ['Condition', 'Sample', 'cluster']).size().reset_index()
        facet_input.rename(columns={0: 'count'}, inplace=True)

        # divide cluster cell counts by total number of cells per sample
        for name, group in facet_input.groupby(['Sample']):

            total_cells = group['count'].sum()

            facet_input.loc[group.index, 'cluster_conc'] = (
                group['count']/total_cells
                )

        # pad counts tables, some clusters may be absent from samples
        pad = pd.DataFrame()
        for cluster in sorted(facet_input['cluster'].unique()):
            to_append = pd.DataFrame(
                {'Sample': natsorted(facet_input['Sample'].unique()),
                 'cluster': [cluster]*len(facet_input['Sample'].unique())})
            pad = pad.append(to_append)
        pad.reset_index(drop=True, inplace=True)
        pad['Condition'] = [sample_metadata[i][0] for i in pad['Sample']]

        facet_input = facet_input.merge(
            pad, how='right', on=['Condition', 'Sample', 'cluster'])
        facet_input.fillna(value=0, inplace=True)

        new_df = pd.DataFrame()
        for cluster in facet_input.groupby(['cluster']):
            for name in cluster[1].groupby(['Condition']):
                cluster[1].loc[
                    name[1]['cluster_conc'].index, 'ave_cluster_conc'
                    ] = name[1]['cluster_conc'].mean()
            conds = cluster[1].drop_duplicates(subset=['Condition']).copy()
            conds.sort_values(
                by='ave_cluster_conc', ascending=False, inplace=True
                )
            conds.reset_index(drop=True, inplace=True)
            conds.reset_index(drop=False, inplace=True)

            conds_dict = dict(zip(conds['Condition'], conds['index']))
            cluster[1]['index'] = [
                conds_dict[i] for i in cluster[1]['Condition']
                ]

            new_df = new_df.append(cluster[1])

        new_df.sort_values(
            by=['cluster', 'ave_cluster_conc'],
            ascending=[True, False], inplace=True)

        # build cmap
        cmap = categorical_cmap(
            numUniqueSamples=len(facet_input['Condition'].unique()),
            numCatagories=10,
            cmap='tab10',
            continuous=False
            )

        sns.set(font_scale=0.4)

        g = sns.FacetGrid(
            data=new_df, col='cluster', col_wrap=10,
            sharex=False, sharey=False, height=1.5, aspect=1.3
            )

        g.map(
            sns.barplot, 'Condition', 'cluster_conc',
            linewidth=0.0, order=None, errwidth=0.2, palette=cmap.colors)

        [
            plt.setp(ax.get_xticklabels(), rotation=90, size=0.5) for
            ax in g.axes.flat
            ]

        new_bar_width = 0.6
        for ax, title in zip(
            g.axes.flatten(), new_df['cluster'].unique()
          ):

            ax.set_title(title, size=6, weight='bold')
            ax.set_xlabel('')

            for patch in ax.patches:
                current_width = patch.get_width()
                diff = current_width - new_bar_width

                # change the bar width
                patch.set_width(new_bar_width)

                # recenter the bar
                patch.set_x(patch.get_x() + diff * 0.5)

        plt.tight_layout()

        plt.savefig(
            os.path.join(self.outDir, 'facetGrid.png'), dpi=800
            )

        plt.show(block=True)

        return self.data

    def frequencyStats(self, args):

        df = read_dataframe(outDir=self.outDir)

        stats_input = df[df['cluster'] >= 0]

        frequency_dir = os.path.join(self.outDir, 'frequency_stats')
        if not os.path.exists(frequency_dir):
            os.mkdir(frequency_dir)

        conditions = sorted(list(set(self.samples.values())))

        # create single-column dataFrame of all sample names
        # to pad counts tables with zeros if a celltype is not in a tissue
        pad = pd.DataFrame(
            sorted(stats_input['sample'].unique())).rename(
                columns={0: 'sample'}
                )

        cluster_list = []
        ratio_list = []
        dif_list = []
        pval_list = []

        # intialize a dataframe to collect catplot data
        catplot_input = pd.DataFrame()

        # loop over clusters
        for w, group in stats_input.groupby('cluster'):

            print(
                f'Calculating log2({conditions[1]}/{conditions[0]})'
                f'of mean cell density for the {str(w)} cluster.'
                )

            group = group.groupby(
                ['sample', 'condition', 'replicate', 'cluster']
                ).size().reset_index(
                drop=False).rename(columns={0: 'count'}).sort_values(
                by='count', ascending=False)

            group = group.merge(pad, how='right', on='sample')

            # guard against NaNs induced by the absence
            # of a given cluster in one or
            # more of the tissue samples
            group['count'] = [
                0 if np.isnan(i) else i for i in group['count']]
            group['condition'] = [
                'cd' if 'cd' in i else 'hfd' for i in group['sample']
                ]
            group['replicate'] = [
                re.sub("\\D", "", i) for i in group['sample']
                ]
            group['cluster'] = w

            # get denominator cell count of each sample
            if self.denominator_cluster is None:
                group['tissue_count'] = [
                    len(stats_input[stats_input['sample'] == i]) for
                    i in group['sample']]
            else:
                group['tissue_count'] = [
                    len(stats_input[(stats_input['sample'] == i) &
                        (stats_input[
                            'cluster'] == self.denominator_cluster)]) for
                    i in group['sample']]

            # compute density of cells per sample corresponding to
            # current variable
            group['density'] = group['count']/group['tissue_count']

            # append group data to catplot_input
            catplot_input = catplot_input.append(group)

            cnd1_values = group['density'][group['condition'] == conditions[0]]
            cnd2_values = group['density'][group['condition'] == conditions[1]]

            stat, pval = ttest_ind(
                cnd1_values, cnd2_values,
                axis=0, equal_var=True, nan_policy='propagate')

            cnd1_mean = np.mean(cnd1_values)
            cnd2_mean = np.mean(cnd2_values)

            ratio = np.log2((cnd2_mean + 0.000001)/(cnd1_mean + 0.000001))

            dif = cnd2_mean-cnd1_mean

            cluster_list.append(w)
            ratio_list.append(ratio)
            dif_list.append(dif)
            pval_list.append(pval)

        statistics = pd.DataFrame(
            list(zip(cluster_list, ratio_list, dif_list, pval_list)),
            columns=['cluster', 'ratio', 'dif', 'pval']).sort_values(
                by='dif')

        statistics.to_csv(
            os.path.join(
                frequency_dir, 'stats.csv'), index=False)

        stats = importr('stats')

        p_adjust = stats.p_adjust(
            FloatVector(statistics['pval'].tolist()),
            method='BH')

        statistics['qval'] = p_adjust

        if self.FDRCorrection:
            stat = 'qval'

        else:
            stat = 'pval'

        significant = statistics[
            statistics[stat] <= 0.05].sort_values(by='dif')

        significant.to_csv(
            os.path.join(
                frequency_dir, 'sig_difs.csv'), index=False)

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.scatter(abs(significant['dif']), significant['ratio'])

        for label, qval, x, y in zip(
          significant['cluster'], significant[stat],
          abs(significant['dif']), significant['ratio']):

            plt.annotate(
                (label, f'{stat[0]}=' + str(round(qval, 4))), size=3,
                xy=(x, y), xytext=(10, 10),
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.1', fc='yellow',
                          alpha=0.0))

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(5)

        plt.title(f'cnd1 vs. cnd2 {stat[0]}<0.05)', fontsize=12)
        plt.xlabel(f'abs({conditions[1]} - {conditions[0]})', fontsize=10)
        plt.ylabel(f'log2({conditions[1]} / {conditions[0]})', fontsize=10)
        plt.savefig(os.path.join(frequency_dir, 'plot.pdf'))
        plt.close()

        catplot_input.reset_index(drop=True, inplace=True)

        catplot_input[stat] = [
             '' if i not in significant['cluster'].unique() else
             round(
                significant[stat][significant['cluster'] == i].values[0], 6)
             for i in catplot_input['cluster']]

        # catplot_input['label'] = catplot_input['cluster'].map(str) + \
        #     ', ' + 'q=' + catplot_input['qval'].astype(str)
        #
        # catplot_input['label'] = [
        #     i.split(',')[0] if not i.split(',')[0] in
        #     significant['cluster'].unique().astype(str) else i for
        #     i in catplot_input['label']]

        catplot_input.sort_values(
            ['cluster', 'condition', 'replicate'], inplace=True)

        # catplot_input.drop('cluster', axis=1, inplace=True)

        # catplot_input.rename(columns={'label': 'cluster'}, inplace=True)

        sns.set(font_scale=0.4)
        g = sns.catplot(
            x='condition', y='density',
            hue='replicate', col='cluster', col_wrap=14,
            data=catplot_input, kind='strip', palette='tab20',
            height=2, aspect=0.8, sharey=False, legend=False)

        g.set(ylim=(0.0, None))
        plt.legend(markerscale=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(frequency_dir, 'catplot.pdf'))
        plt.close('all')

    def clusterBoxplots(self, args):

        cmap = categorical_cmap(
            numCatagories=10, numSubcatagories=2,
            cmap='tab10', continuous=False
            )

        df = read_dataframe(outDir=self.outDir)

        boxplot_input = df[df['cluster'] >= 0]

        # get tidy input data
        boxplot_input = (
            boxplot_input[[
                'cluster', 'sample', 'condition'] + self.markers]
            .melt(
                id_vars=['cluster', 'sample', 'condition'],
                var_name='protein', value_name='log10(intensity)'))

        boxplot_input = boxplot_input.sample(frac=1.0)

        boxplot_input['label'] = (
            boxplot_input['protein'] + '_'
            + boxplot_input['cluster'].map(str) + '_'
            + boxplot_input['condition']
            )

        boxplot_input.sort_values(by='label', inplace=True)

        boxplot_input.rename(
            columns={'log10(intensity)': '$log_{10}(intensity)$'},
            inplace=True)

        sns.set(font_scale=0.27)
        sns.set_style('whitegrid')

        g = sns.FacetGrid(
            boxplot_input, row='cluster', col='protein',
            sharex=False, sharey=False, height=1, aspect=1.5)

        hue_dict = dict(
            zip(boxplot_input['cluster'].unique(), cmap.colors)
            )

        g.map(
            lambda x, y, z, color:
                sns.boxplot(
                    data=boxplot_input, x=x, y=y,
                    hue=z, palette=hue_dict,
                    linewidth=0.95, width=0.75,
                    fliersize=0.95),
                'condition', '$log_{10}(intensity)$', 'cluster')

        def statistics(label):

            conditions = sorted(list(set(self.samples.values())))
            cond_dict = {
                conditions[0]: conditions[1], conditions[1]: conditions[0]
                }
            label = label.values[0]
            cond_name = label.split('_')[-1]
            opposite_label = label.replace(cond_name, cond_dict[cond_name])

            means1 = boxplot_input[
                    boxplot_input['label'] == label].groupby(
                        'sample').mean()['$log_{10}(intensity)$']

            means2 = boxplot_input[
                    boxplot_input['label'] == opposite_label].groupby(
                        'sample').mean()['$log_{10}(intensity)$']

            # perform Welch's unequal variances t-test
            stat, pval = ttest_ind(
                means1, means2,
                axis=0, equal_var=False, nan_policy='propagate')

            if self.bonferroniCorrection:

                # perform Bonferroni correction
                p_adj = pval * len(boxplot_input['label'].unique())/2

                if p_adj <= 0.05:
                    ax = plt.gca()
                    ax.text(
                        0.5, 0.85,
                        r'$p_{adj} = $' + '%.1E' % Decimal(str(p_adj)),
                        fontweight='normal', fontsize=11.0,
                        color='k', ha='center', va='center',
                        transform=ax.transAxes
                        )
            else:
                # DO NOT perform Bonferroni correction
                if pval <= 0.05:
                    ax = plt.gca()
                    ax.text(
                        0.5, 0.85, r'$p = $' + '%.1E' % Decimal(str(pval)),
                        fontweight='normal', fontsize=11.0,
                        color='k', ha='center', va='center',
                        transform=ax.transAxes
                        )

        g.map(
            lambda label, color:
                statistics(label=label), 'label')

        plt.savefig(
            os.path.join(self.outDir, 'cluster_boxplots.pdf'), bbox='tight')
        plt.close('all')

    def curateThumbnails(self, args):

        df = self.data

        markers, dna1, dna_prefix = read_markers(
            markers_filepath=self.markers_filepath
            )

        abx_channels = [
            f'{i}_cellMask' for i in markers['marker_name'] if
            dna_prefix not in i
            ]

        zs = loadZarrs(
            df=df, outDir=self.outDir, markers_filepath=self.markers_filepath
            )

        thumbnail_input = df[df['cluster'] >= 0]

        thumbnails_dir = os.path.join(self.outDir, 'thumbnails')
        if not os.path.exists(thumbnails_dir):
            os.mkdir(thumbnails_dir)

        for cluster in sorted(thumbnail_input['cluster'].unique()):

            print(cluster)

            markers_to_show = cluster_expression(
                df=thumbnail_input, markers=abx_channels,
                cluster=cluster, num_proteins=3
                )

            markers_to_show = [
                '_'.join(i.split('_')[0:-1]) if '_' in i else i for i in
                [dna1] + markers_to_show
                ]

            color_dict = {}
            for i, j, k in zip(
              markers_to_show,

              [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0),
               (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
              ['gray', 'green', 'red', 'blue']
              ):
                color_dict[i] = j

            long_table = pd.DataFrame()

            for sample in thumbnail_input['Sample'].unique():

                if '-' in sample:
                    img_name = sample.split('-')[1]
                else:
                    img_name = sample

                # initialize a zeros-array with the dimensions of the image
                overlay = np.zeros(
                    (zs[f'{img_name}_{dna1}'][:].shape[0],
                     zs[f'{img_name}_{dna1}'][:].shape[1]))

                # convert to rgb
                overlay = gray2rgb(overlay)

                # loop over the channels to create a dict of images
                for marker in markers_to_show:
                    marker_img = img_as_float(
                        zs[f'{img_name}_{marker}']
                        )
                    marker_img = gray2rgb(marker_img)
                    marker_img = (
                        marker_img * color_dict[marker]
                        )

                    # loop over channels to add each image to overlay seed
                    overlay += marker_img

                # crop out thumbnail images
                sample_cluster_subset = thumbnail_input[
                    (thumbnail_input['Sample'] == sample)
                    & (thumbnail_input['cluster'] == cluster)
                    ]

                sample_cluster_subset.reset_index(drop=True, inplace=True)

                if self.numFingernails > len(sample_cluster_subset):
                    dif = self.numFingernails - len(sample_cluster_subset)
                    extra_rows = pd.DataFrame(
                        data=0,
                        index=list(range(dif)),
                        columns=sample_cluster_subset.columns
                        )
                    sample_cluster_subset = sample_cluster_subset.append(
                        extra_rows
                        )
                    sample_cluster_subset.reset_index(
                        drop=True, inplace=True
                        )
                else:
                    sample_cluster_subset = sample_cluster_subset.sample(
                        n=self.numFingernails, random_state=3
                        )

                # add centroid mask to image overlay
                centroids = sample_cluster_subset[['X_centroid', 'Y_centroid']]

                centroid_img = np.zeros(
                    (overlay.shape[0],
                     overlay.shape[1]))

                centroid_dist = 1  # in pixels
                for example, centroid in enumerate(centroids.iterrows()):

                    ystart_centroid = int(
                        centroid[1]['Y_centroid'] - centroid_dist
                        )
                    ystop_centroid = int(
                        centroid[1]['Y_centroid'] + centroid_dist
                        )

                    xstart_centroid = int(
                        centroid[1]['X_centroid'] - centroid_dist
                        )
                    xstop_centroid = int(
                        centroid[1]['X_centroid'] + centroid_dist
                        )

                    centroid_img[
                        ystart_centroid:ystop_centroid,
                        xstart_centroid:xstop_centroid
                        ] = 1

                # convert to rgb and colorize
                centroid_img = gray2rgb(centroid_img)
                centroid_img = (centroid_img * (1.0, 1.0, 1.0))

                # add to overlay
                overlay += centroid_img

                # crop thumbnails
                window_dist = 35  # in pixels
                for example, centroid in enumerate(centroids.iterrows()):

                    if (
                        (centroid[1]['X_centroid'] == 0.0) &
                        (centroid[1]['Y_centroid'] == 0.0)
                    ):

                        blank_img = np.ones(
                            (window_dist,
                             window_dist))

                        long_table = long_table.append(
                            {'sample': sample, 'example': int(example),
                             'image': blank_img},
                            ignore_index=True
                            )

                    else:

                        # specify window x, y ranges
                        ystart_window = int(
                            centroid[1]['Y_centroid'] - window_dist
                            )
                        ystop_window = int(
                            centroid[1]['Y_centroid'] + window_dist
                            )

                        xstart_window = int(
                            centroid[1]['X_centroid'] - window_dist
                            )
                        xstop_window = int(
                            centroid[1]['X_centroid'] + window_dist
                            )

                        # crop overlay image to window size
                        thumbnail = overlay[
                            ystart_window:ystop_window,
                            xstart_window:xstop_window
                            ]

                        long_table = long_table.append(
                            {'sample': sample, 'example': int(example),
                             'image': thumbnail},
                            ignore_index=True
                            )
            print()

            # plot facet grid
            fig, ax = plt.subplots()

            g = sns.FacetGrid(
                long_table, row='sample', col='example',
                sharex=False, sharey=False,
                gridspec_kws={'hspace': 0.1, 'wspace': 0.05})

            g.map(
                lambda x, **kwargs: (
                    plt.imshow(x.values[0]), plt.grid(False)), 'image')

            for ax in g.axes.flatten():
                ax.get_xaxis().set_ticks([])
                ax.set_xlabel('')
                ax.get_yaxis().set_ticks([])
                ax.set_ylabel('')

            g.set_titles(
                col_template="{col_name}", row_template="{row_name}",
                fontweight='bold', size=8)

            custom_lines = []
            for k, v in color_dict.items():
                custom_lines.append(
                    Line2D([0], [0], color=v, lw=6))

            ax.legend(
                custom_lines, list(color_dict.keys()), prop={'size': 12},
                loc='upper right', bbox_to_anchor=(1.0, 1.0)
                )

            plt.savefig(
                os.path.join(
                    thumbnails_dir,
                    'cluster' + str(cluster) + '_thumbnails.pdf'),
                bbox_inches='tight')
            plt.close('all')
            del g
            gc.collect()

        return self.data

    def spatialAnalysis(self, args):

        df = read_dataframe(outDir=self.outDir)

        zs = loadZarrs(df=df, outDir=self.outDir)

        spatial_dir = os.path.join(self.outDir, 'spatial_analysis')
        if not os.path.exists(spatial_dir):
            os.makedirs(spatial_dir)

        stats = pd.DataFrame(
            columns=[
                'protein', 'sample', 'celltype', 'ratio_r',
                'centroids', 'ratio_f', 'points']
            )

        stats_row_idx = 0
        for protein, binary_cutoff in self.spatialDict1.items():
            for sample in df['sample'].unique():

                if sample in self.cropDict.keys():
                    section = self.cropDict[sample][0]
                    cut_point = self.cropDict[sample][1]

                    if section == 'bottom':
                        dna = zs[f'{sample}_dna'][cut_point:]
                        img = zs[f'{sample}_{protein}'][cut_point:]

                    elif section == 'top':
                        dna = zs[f'{sample}_dna'][:cut_point]
                        img = zs[f'{sample}_{protein}'][:cut_point]

                else:
                    dna = zs[f'{sample}_dna'][:]
                    img = zs[f'{sample}_{protein}'][:]

                for celltype, cluster in self.spatialDict2.items():

                    print(celltype, protein, sample)

                    total_centroids = df[['x', 'y']][
                        (df['sample'] == sample) &
                        (df['cluster'] == cluster)].astype(int)

                    if sample in self.cropDict.keys():
                        section = self.cropDict[sample][0]
                        cut_point = self.cropDict[sample][1]

                        if section == 'bottom':
                            total_centroids = total_centroids[
                                total_centroids['y'] >= cut_point]

                            total_centroids['y'] = (
                                total_centroids['y']-cut_point
                                )

                        elif section == 'top':
                            total_centroids = total_centroids[
                                total_centroids['y'] < cut_point]

                    if len(total_centroids) > 1:
                        y_min = total_centroids['y'].min()
                        y_max = total_centroids['y'].max()
                        y_range = y_max - y_min

                        print(y_min, y_max)

                        x_min = total_centroids['x'].min()
                        x_max = total_centroids['x'].max()
                        x_range = x_max - x_min

                        dna_blurred = gaussian(dna, sigma=12)
                        dna_mask = np.where(dna_blurred > 0.05, 1, 0)

                        area_mask = dna_mask[y_min:y_max, :]
                        area_mask = area_mask[:, x_min:x_max]

                        inside_tumor = np.argwhere(area_mask == 1)
                        outside_tumor = np.argwhere(area_mask == 0)

                        frac_area_out = outside_tumor.shape[0]/(
                            outside_tumor.shape[0] + inside_tumor.shape[0])

                        img_blurred = gaussian(img, sigma=12)
                        img_mask_real = np.where(
                            img_blurred > binary_cutoff, 1, 0
                            )
                        img_mask_fake = img_mask_real.copy()

                        (img_mask_real[total_centroids['y'],
                         total_centroids['x']]
                         ) += 10

                        outside_centroids = np.argwhere(img_mask_real == 10)
                        outside_centroids = pd.DataFrame(
                            outside_centroids, columns=['y', 'x'])

                        inside_centroids = np.argwhere(img_mask_real == 11)
                        inside_centroids = pd.DataFrame(
                            inside_centroids, columns=['y', 'x'])

                        radii = []
                        num_points = []
                        for radius in range(
                            self.radiusRange[0], self.radiusRange[1]
                          ):
                            total_poisson_points = pd.DataFrame(
                                poisson_disc_samples(
                                    width=x_range, height=y_range, r=radius),
                                columns=['x', 'y'])
                            print(radius, len(total_poisson_points))
                            radii.append(radius)
                            num_points.append(len(total_poisson_points))

                        def closest(lst, K):
                            return lst[
                                min(
                                    range(
                                        len(lst)),
                                    key=lambda i: abs(lst[i]-K))]

                        optimal_points = closest(
                            num_points, len(total_centroids) +
                            (int(len(total_centroids) * frac_area_out)))
                        optimal_radius = radii[
                            num_points.index(optimal_points)-2
                            ]

                        total_poisson_points = pd.DataFrame(
                            poisson_disc_samples(
                                width=x_range, height=y_range,
                                r=optimal_radius),
                            columns=['x', 'y']).astype(int)

                        # ensure simulation contains at least 2 points
                        while len(total_poisson_points) < 2:
                            optimal_radius -= 1
                            total_poisson_points = pd.DataFrame(
                                poisson_disc_samples(
                                    width=x_range, height=y_range,
                                    r=optimal_radius),
                                columns=['x', 'y']).astype(int)

                        total_poisson_points['x'] = (
                            total_poisson_points['x'] + x_min
                            )
                        total_poisson_points['y'] = (
                            total_poisson_points['y'] + y_min
                            )

                        (dna_mask[total_poisson_points['y'],
                                  total_poisson_points['x']]
                         ) += 10

                        total_poisson_points = np.argwhere(dna_mask == 11)
                        total_poisson_points = pd.DataFrame(
                            total_poisson_points, columns=['y', 'x'])

                        (img_mask_fake[total_poisson_points['y'],
                         total_poisson_points['x']]
                         ) += 10

                        outside_poisson_points = np.argwhere(
                            img_mask_fake == 10
                            )
                        outside_poisson_points = pd.DataFrame(
                            outside_poisson_points, columns=['y', 'x'])

                        inside_poisson_points = np.argwhere(
                            img_mask_fake == 11
                            )
                        inside_poisson_points = pd.DataFrame(
                            inside_poisson_points, columns=['y', 'x'])

                        rgb_img = img_as_float(img)
                        rgb_img = gray2rgb(rgb_img)
                        rgb_img = (rgb_img * (0.5, 0.5, 0.5))

                        plt.imshow(rgb_img)
                        plt.scatter(
                            inside_centroids['x'],
                            inside_centroids['y'], s=0.5, ec='none', c='g'
                            )

                        plt.scatter(
                            outside_centroids['x'],
                            outside_centroids['y'], s=0.5, ec='none', c='r')

                        legend_elements = []
                        legend_elements.append(
                            Line2D([0], [0], marker='o', color='none',
                                   label='inside',
                                   markerfacecolor='g',
                                   markeredgecolor='none', lw=0.001,
                                   markersize=6)
                                   )
                        legend_elements.append(
                            Line2D([0], [0], marker='o', color='none',
                                   label='outside',
                                   markerfacecolor='r',
                                   markeredgecolor='none', lw=0.001,
                                   markersize=6)
                                   )

                        plt.legend(
                            handles=legend_elements, prop={'size': 6},
                            bbox_to_anchor=[1.0, 1.0])

                        ratio_real = str(
                            round(
                                len(inside_centroids)/len(total_centroids), 2)
                                )
                        total_cells = str(len(total_centroids))
                        title_statement = (
                            f'{sample}-{celltype}, inside/total: {ratio_real},'
                            + f' total cells={total_cells}'
                            )
                        plt.title(title_statement, fontsize=10)

                        plt.grid(False)
                        plt.savefig(
                            os.path.join(
                                spatial_dir,
                                f'{protein}_{celltype}_{sample}.png'), dpi=800)
                        plt.close('all')

                        plt.imshow(rgb_img)

                        plt.scatter(
                            inside_poisson_points['x'],
                            inside_poisson_points['y'],
                            s=0.5, ec='none', c='g'
                            )

                        plt.scatter(
                            outside_poisson_points['x'],
                            outside_poisson_points['y'],
                            s=0.5, ec='none', c='r'
                            )

                        legend_elements = []
                        legend_elements.append(
                            Line2D([0], [0], marker='o', color='none',
                                   label='inside',
                                   markerfacecolor='g',
                                   markeredgecolor='none', lw=0.001,
                                   markersize=6)
                                   )
                        legend_elements.append(
                            Line2D([0], [0], marker='o', color='none',
                                   label='outside',
                                   markerfacecolor='r',
                                   markeredgecolor='none', lw=0.001,
                                   markersize=6)
                                   )

                        plt.legend(
                            handles=legend_elements, prop={'size': 6},
                            bbox_to_anchor=[1.0, 1.0])

                        if not len(total_poisson_points) == 0:
                            ratio_fake = str(
                                round(
                                    len(inside_poisson_points) /
                                    len(total_poisson_points), 2))
                        else:
                            ratio_fake = 'divide by zero'

                        poisson_points = str(len(total_poisson_points))
                        title_statement = (
                            f'{sample}-{celltype}_Poisson-disc,' +
                            f' inside/total: {ratio_fake},' +
                            f' total points={poisson_points}'
                            )
                        plt.title(title_statement, fontsize=10)

                        plt.grid(False)
                        plt.savefig(
                            os.path.join(
                                spatial_dir,
                                f'{protein}_{celltype}_Poisson-disc_' +
                                f'{sample}.png'),
                            dpi=800
                            )
                        plt.close('all')

                        stats.loc[stats_row_idx] = (
                            protein, sample, celltype, float(ratio_real),
                            len(total_centroids), float(ratio_fake),
                            len(total_poisson_points)
                            )
                        stats_row_idx += 1

        stats.to_csv(os.path.join(spatial_dir, 'stats.csv'))
        # stats = pd.read_csv(
        #     os.path.join(spatial_dir, 'stats.csv'), index_col=0)

        stats = stats[stats['ratio_f'] != 'divide by zero']
        stats['ratio_f'] = [float(i) for i in stats['ratio_f']]

        stats['condition'] = [
            re.findall(r"[^\W\d_]+|\d+", i)[0] for
            i in stats['sample']
            ]

        simulation = stats.copy()

        stats.rename(
            columns={
                'centroids': '# measured cells',
                'points': '# simulated cells'},
            inplace=True
            )

        sns.set_style('whitegrid')
        sns.scatterplot(
            x=stats['# simulated cells'],
            y=stats['# measured cells'],
            color='k',
            data=stats
            )
        plt.savefig(
            os.path.join(
                spatial_dir, 'points_v_centroids.pdf'),
            dpi=600, bbox_inches='tight')
        plt.close('all')

        stats[
            '% cells overlapping immunomarker (normalized to simulation)'] = (
            stats['ratio_r'] - stats['ratio_f']
            )

        stats.sort_values(
            by=['protein', 'celltype', 'condition'],
            ascending=[True, False, True],
            inplace=True
            )
        stats.reset_index(drop=True, inplace=True)

        stats['label'] = (
            stats['protein'] + '_' +
            stats['celltype']
            )

        sns.set_style('whitegrid')
        sns.swarmplot(
            x='label',
            y='% cells overlapping immunomarker (normalized to simulation)',
            hue='condition',
            size=3,
            dodge=True,
            palette=['lightgray', 'firebrick'],
            data=stats
            )
        plt.xticks(rotation=90, size=5)

        plt.savefig(
            os.path.join(
                spatial_dir,
                'percent cells overlapping immunomarker' +
                ' (normalized to simulation).pdf'),
            dpi=600, bbox_inches='tight')
        plt.close('all')

        condition_sig = {}
        for name, group in stats.groupby(['protein', 'celltype']):
            a = group[
                '% cells overlapping immunomarker (normalized to simulation)'][
                    group['condition'] == 'cd']
            b = group[
                '% cells overlapping immunomarker (normalized to simulation)'][
                    group['condition'] == 'hfd']
            t, pval = ttest_ind(
                a, b, axis=0, equal_var=True, nan_policy='propagate')
            condition_sig[f'{name}'] = round(pval, 3)

        condition_sig_df = pd.DataFrame(condition_sig, index=range(0, 27)).T[0]
        condition_sig_df.rename('pval', inplace=True)
        condition_sig_df.to_csv(
            os.path.join(spatial_dir, 'treatment_sig.csv'), header=True)

        stats.to_csv(
            os.path.join(spatial_dir, 'stats_normalized.csv'), index=False
            )

        simulation['dtype'] = 'simulation'

        measurement = simulation.copy()
        measurement['dtype'] = 'measurement'

        cols = [
            col for col in simulation.columns if
            col not in ['ratio_r', 'centroids']
            ]
        simulation = simulation[cols]
        simulation.rename(
            columns={'ratio_f': '% cells overlapping immunomarker'},
            inplace=True)

        cols = [
            col for col in measurement.columns if
            col not in ['ratio_f', 'points']
            ]
        measurement = measurement[cols]

        measurement.rename(
            columns={
                'ratio_r': '% cells overlapping immunomarker',
                'centroids': 'points'}, inplace=True
                )

        stats = simulation.append(measurement, sort=True, ignore_index=True)
        del simulation
        del measurement

        stats['% cells overlapping immunomarker'] = (
            stats['% cells overlapping immunomarker'] * 100
            )

        stats.sort_values(
            by=['protein', 'celltype', 'condition', 'dtype'],
            ascending=[True, False, True, False],
            inplace=True
            )
        stats.reset_index(drop=True, inplace=True)

        stats['label'] = (
            stats['protein'] + '_' +
            stats['celltype'] + '_' +
            stats['condition'] + '_' +
            stats['dtype']
            )

        cols = [
            'protein', 'celltype', 'condition', 'dtype',
            '% cells overlapping immunomarker', 'points',
            'sample', 'label'
            ]

        stats = stats[cols]

        sns.set_style('whitegrid')
        sns.swarmplot(
            x='label',
            y='% cells overlapping immunomarker',
            hue='dtype',
            size=3,
            dodge=False,
            palette=['lightgray', 'firebrick'],
            data=stats
            )
        plt.xticks(rotation=90, size=5)

        plt.savefig(
            os.path.join(
                spatial_dir, 'percent cells overlapping immunomarker.pdf'),
            dpi=600, bbox_inches='tight')
        plt.close('all')

        model_sig = {}
        for name, group in stats.groupby(['protein', 'celltype', 'condition']):
            a = group[
                '% cells overlapping immunomarker'][
                    group['dtype'] == 'simulation']
            b = group[
                '% cells overlapping immunomarker'][
                    group['dtype'] == 'measurement']
            t, pval = ttest_ind(
                a, b, axis=0, equal_var=True, nan_policy='propagate')
            model_sig[f'{name}'] = round(pval, 3)

        model_sig_df = pd.DataFrame(model_sig, index=range(0, 27)).T[0]
        model_sig_df.rename('pval', inplace=True)
        model_sig_df.to_csv(
            os.path.join(spatial_dir, 'model_sig.csv'), header=True
            )

        stats.to_csv(os.path.join(
            spatial_dir, 'stats_unnormalized.csv'), index=False
            )
