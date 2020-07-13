
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io
import pickle
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.widgets import Slider, Button
from matplotlib.widgets import TextBox
from matplotlib.text import Text
import re
import subprocess
from natsort import natsorted, order_by_index, index_natsorted
from sklearn.preprocessing import MinMaxScaler
import yaml
import zarr
from numcodecs import Blosc
from skimage.util.dtype import img_as_float
from skimage.util.dtype import img_as_uint
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize as norm
from datetime import datetime
from sklearn.manifold import TSNE
import hdbscan
from joblib import Memory
from matplotlib.colors import ListedColormap
import math
from scipy.stats import ttest_ind
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from decimal import Decimal
from skimage.color import gray2rgb
import gc
from skimage.filters import gaussian
from bridson import poisson_disc_samples

from utilsQC import *

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
                 markers=None,
                 samples=None,
                 replicates=None,
                 cycleConfig=None,
                 omeroSettings=None,
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
                 denominator_clusters=2,
                 FDRCorrection=False,
                 BonferroniCorrection=False,
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
                    'CD11B_MYC': 2, 'CD8T': 7
                    },
                 radiusRange=[40, 600],
                 dfSaveCount=1
                 ):
        """
        Args:
            input_dir: directory containing single-cell data files
            output_dir: directory save QC output files
            markers: molecular probes to be considered in the analysis
        """

        # assert(SOMETHING)  # placeholder for now

        self.inDir = inDir
        self.outDir = outDir
        self.markers = markers
        self.samples = samples
        self.replicates = replicates
        self.cycleConfig = cycleConfig
        self.omeroSettings = omeroSettings
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
        self.denominator_clusters = denominator_clusters
        self.FDRCorrection = FDRCorrection
        self.BonferroniCorrection = BonferroniCorrection
        self.numFingernails = numFingernails
        self.cropDict = cropDict
        self.spatialDict1 = spatialDict1
        self.spatialDict2 = spatialDict2
        self.radiusRange = radiusRange
        self.dfSaveCount = dfSaveCount

    def getSingleCellData(self, args):

        files = dataset_files(self.inDir)

        sample_to_condition = self.samples
        sample_to_replicate = self.replicates

        df_list = []
        for file in files:
            sample = file.split('.')[0]

            print(f'Importing single-cell data for sample {sample}.')

            data = pd.read_csv(os.path.join(self.inDir, file))

            data['sample'] = sample

            # append dataframe to list
            df_list.append(data)
            del data

        # stack dataframes row-wise
        df = pd.concat(df_list, axis=0)
        del df_list

        # assign global index
        df.reset_index(drop=True, inplace=True)

        # add condition column
        df['condition'] = [
            sample_to_condition[s] for s in df['sample']]

        # add replicate column
        df['replicate'] = [
            int(sample_to_replicate[s]) for s in df['sample']]

        # reindex dna cycles starting at 1, not 0
        dna_dict = {}
        idx_seed = 1
        for cycle in df.columns.sort_values():
            if 'dna_cycle' in cycle:
                dna_dict[cycle] = 'dna_cycle' + str(idx_seed)
                idx_seed += 1

        df.rename(columns=dna_dict, inplace=True)

        # organize columns
        cols = [
            'x', 'y', 'area', 'sample', 'condition', 'replicate', 'mask_id',
            'dna_cycle1', 'dna_cycle2', 'dna_cycle3', 'dna_cycle4',
            'dna_cycle5', 'dna_cycle6', 'dna_cycle7', 'dna_cycle8',
            'dna_cycle9', 'dna_cycle10'] + self.markers
        df = df[cols]

        self.dfSaveCount = save_dataframe(
            df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
            )

    def lassoROIs(self, args):

        df = read_dataframe(outDir=self.outDir)

        if not os.path.exists(os.path.join(self.outDir, 'lasso_dict.pkl')):

            lasso_dict = {}
            for name, group in df.groupby('sample'):
                data = df[['x', 'y']][df['sample'] == name]
                data.reset_index(drop=False, inplace=True)

                subplot_kw = dict(
                    xlim=(data['x'].min(), data['x'].max()),
                    ylim=(data['y'].min(), data['y'].max()),
                    autoscale_on=False)

                subplot_kw = dict()

                fig, ax = plt.subplots(subplot_kw=subplot_kw)

                # read ashlar_output dna channel
                dna = io.imread(f'{self.inDir}/{name}_dna.tif')

                ax.imshow(dna, cmap='gray')

                plt.grid(False)
                pts = ax.scatter(data['x'], data['y'], c='lime', s=0.0005)
                selector = SelectFromCollection(ax, pts)

                def accept(event):
                    if event.key == "enter":
                        print("Selected points:")
                        print(selector.xys[selector.ind])
                        selector.disconnect()
                        ax.set_title("")
                        fig.canvas.draw()

                fig.canvas.mpl_connect("key_press_event", accept)
                ax.set_title(
                    f'Sample {name}. Press enter to accept selected points.')
                ax.set_aspect('equal')
                plt.show(block=True)

                # filter dataframe
                drop_idx = (
                    set(np.array(list(range(selector.xys.data.shape[0]))))
                    - set(selector.ind))

                idx = data['index'].iloc[list(drop_idx)].values
                df.drop(idx, inplace=True)
                lasso_dict[name] = idx

            os.chdir(self.outDir)
            f = open('lasso_dict.pkl', 'wb')
            pickle.dump(lasso_dict, f)
            f.close()

            self.dfSaveCount = save_dataframe(
                df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
                )

        else:

            os.chdir(self.outDir)
            pickle_in = open('lasso_dict.pkl', 'rb')
            lasso_dict = pickle.load(pickle_in)

            for key, value in lasso_dict.items():
                df.drop(lasso_dict[key], inplace=True)

            self.dfSaveCount = save_dataframe(
                df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
                )

    def dnaIntensityCutoff(self, args):

        df = read_dataframe(outDir=self.outDir)

        bins = 100
        rnge = [0, 0.55]
        histtype = 'stepfilled'

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        plt.hist(
            df['dna_cycle1'], bins=bins,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype=histtype,
            range=rnge, label='before'
            )

        plt.title('median DNA intensity')
        plt.ylabel('count')

        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

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
                (df['dna_cycle1'] > lowerCutoff) &
                (df['dna_cycle1'] < upperCutoff)
                ]

            if text in df['sample'].unique():
                dna = io.imread(
                    f'{self.inDir}/{text}_dna.tif')

                fig, ax = plt.subplots()
                ax.imshow(dna, cmap='gray')
                ax.grid(False)
                coords = df_test[
                    ['x', 'y', 'dna_cycle1']][df_test['sample'] == text]
                sp = ax.scatter(
                    coords['x'], coords['y'], s=coords['dna_cycle1']*40,
                    c=coords['dna_cycle1'], cmap='viridis')
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
        text_box.on_text_change(lambda val: submit(val, df))
        plt.show(block=True)

        lowerCutoff, upperCutoff = update(val=None)

        fig, ax = plt.subplots()
        # plot cell area histogram BEFORE filtering
        plt.hist(
            df['dna_cycle1'], bins=bins,
            density=False, color='b', ec='none',
            alpha=0.5, histtype=histtype,
            range=rnge, label='before'
            )

        # apply lower and upper cutoffs
        df = df[
            (df['dna_cycle1'] > lowerCutoff) &
            (df['dna_cycle1'] < upperCutoff)
            ]

        # plot cell area histogram AFTER filtering
        plt.hist(
            df['dna_cycle1'], bins=bins, color='r', ec='none',
            alpha=0.5, histtype=histtype, range=rnge, label='after')
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

        self.dfSaveCount = save_dataframe(
            df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
            )

    def nuclearAreaCutoff(self, args):

        df = read_dataframe(outDir=self.outDir)

        bins = 100
        rnge = [0, 1000]
        histtype = 'stepfilled'

        sns.set_style('whitegrid')
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.25, bottom=0.25)

        plt.hist(
            df['area'], bins=bins,
            density=False, color='grey', ec='none',
            alpha=0.75, histtype=histtype,
            range=rnge, label='before'
            )

        plt.title('nuclear area')
        plt.ylabel('count')

        axcolor = 'lightgoldenrodyellow'
        axLowerCutoff = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
        axUpperCutoff = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

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
                (df['area'] > lowerCutoff) &
                (df['area'] < upperCutoff)
                ]

            if text in df['sample'].unique():
                dna = io.imread(
                    f'{self.inDir}/{text}_dna.tif')

                fig, ax = plt.subplots()
                ax.imshow(dna, cmap='gray')
                ax.grid(False)
                coords = df_test[
                    ['x', 'y', 'area']][df_test['sample'] == text]
                sp = ax.scatter(
                    coords['x'], coords['y'], s=coords['area']/100,
                    c=coords['area'], cmap='viridis')
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
        text_box.on_text_change(lambda val: submit(val, df))
        plt.show(block=True)

        lowerCutoff, upperCutoff = update(val=None)

        fig, ax = plt.subplots()
        # plot cell area histogram BEFORE filtering
        plt.hist(
            df['area'], bins=bins,
            density=False, color='b', ec='none',
            alpha=0.5, histtype=histtype,
            range=rnge, label='before'
            )

        # apply lower and upper cutoffs
        df = df[
            (df['area'] > lowerCutoff) &
            (df['area'] < upperCutoff)
            ]

        # plot cell area histogram AFTER filtering
        plt.hist(
            df['area'], bins=bins, color='r', ec='none',
            alpha=0.5, histtype=histtype, range=rnge, label='after')
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

        for sample in df['sample'].unique():

            dna = io.imread(
                f'{self.inDir}/{sample}_dna.tif')

            fig, ax = plt.subplots()
            ax.imshow(dna, cmap='gray')
            ax.grid(False)
            coords = df[['x', 'y', 'area']][df['sample'] == sample]
            sp = ax.scatter(
                coords['x'], coords['y'], s=coords['area']/100,
                c=coords['area'], cmap='viridis')
            plt.title(
                f'Sample {sample}. Lasso colored by nuclear area')
            plt.colorbar(sp)
            plt.savefig(
                os.path.join(
                    lasso_dir, f'{sample}.png'), dpi=1000)
            plt.close('all')

        self.dfSaveCount = save_dataframe(
            df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
            )

    def crossCyleCorrelation(self, args):

        df = read_dataframe(outDir=self.outDir)

        ratios = pd.DataFrame(
            [np.log10(
                (df['dna_cycle1'] + 0.00001) /
                (df[i] + 0.00001)) for i in
                df.columns[df.columns.str.contains('dna_cycle')]]).T

        list1 = [i for i in ratios.columns if i.startswith('Unnamed')]
        list2 = [f'1/{i+1}' for i in range(1, len(list1)+1)]
        ratio_columns = dict(zip(list1, list2))
        ratios.rename(columns=ratio_columns, inplace=True)
        ratios.drop('dna_cycle1', axis=1, inplace=True)
        ratios['sample'] = df['sample']

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

        # get dataframe index
        total_indices = df.index

        # grab dna and sample columns of dataframe
        facet_input = df.loc[:, df.columns.str.contains('dna_|sample')]

        # initialize a set to append indices to drop
        indices_to_drop = set()

        # loop over samples
        for sample in df['sample'].unique():

            # slice sample-specific data
            sample_df = df[df['sample'] == sample]

            # loop over dna cycle columns
            for col_name in sample_df.columns:
                if 'dna_cycle' in col_name:

                    # get cycle number
                    cycle_num = str(re.search(r'\d+', col_name).group())

                    if cycle_num != '1':

                        # get ratios
                        ratios = np.log10(
                            (sample_df['dna_cycle1'] + 0.00001) /
                            (sample_df['dna_cycle' + cycle_num] + 0.00001))

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
            :, df.columns.str.contains('dna_|sample')].copy()

        # plot cell dropout facet (per cycle)
        facet_per_cycle_melt = (
            facet_input.drop(['sample'], axis=1)
            .melt(var_name='cycle', value_name='signal'))

        fig, ax = plt.subplots(figsize=(5, 5))
        g = sns.FacetGrid(facet_per_cycle_melt, col='cycle', col_wrap=4)

        g.map(
            lambda y, color: plt.scatter(
                y, facet_input['dna_cycle1'],
                s=0.25, alpha=0.1, linewidth=None,
                marker='o', c='r'), 'signal')
        plt.savefig(
            os.path.join(
                self.outDir, 'cycle_correlation(perCycle).png'), dpi=1000)
        plt.close('all')

        # plot cell dropout facet (per sample per cycle)
        # take a fraction of the total dataframe for plotting
        facet_per_sample_per_cycle_melt = (
            facet_input.sample(frac=0.1)
            .melt(id_vars=['sample'],
                  var_name='cycle', value_name='signal')
            )

        facet_per_sample_per_cycle_melt['sample'] = pd.Categorical(
            facet_per_sample_per_cycle_melt['sample'], ordered=True,
            categories=natsorted(
                facet_per_sample_per_cycle_melt['sample'].unique()))

        facet_per_sample_per_cycle_melt['cycle'] = pd.Categorical(
            facet_per_sample_per_cycle_melt['cycle'], ordered=True,
            categories=natsorted(
                facet_per_sample_per_cycle_melt['cycle'].unique()))

        facet_per_sample_per_cycle_melt = (
            facet_per_sample_per_cycle_melt.sort_values(['sample', 'cycle'])
            )

        fig, ax = plt.subplots(figsize=(5, 5))

        cmap = categorical_cmap(
            numCatagories=10, numSubcatagories=3,
            cmap='tab10', continuous=False
            )

        # trim qualitative cmap to the number of unique samples
        trim = len(cmap.colors) - len(
            facet_per_sample_per_cycle_melt['sample'].unique())

        sample_color_dict = dict(
            zip(
                natsorted(facet_per_sample_per_cycle_melt['sample'].unique()),
                cmap.colors[:-trim])
                )

        g = sns.FacetGrid(
            facet_per_sample_per_cycle_melt, col='cycle', hue='sample'
            )

        g.map(
            lambda sam, y, color, **kwargs: plt.scatter(
                y, facet_per_sample_per_cycle_melt.loc[
                    (facet_per_sample_per_cycle_melt['sample'] ==
                     sam.unique()[0])
                    & (facet_per_sample_per_cycle_melt['cycle'] ==
                       'dna_cycle1'), 'signal'],
                c=np.reshape(sample_color_dict[sam.unique()[0]], (-1, 3)),
                s=0.25, linewidth=None, marker='o', **kwargs),
            'sample', 'signal'
            )

        plt.legend(markerscale=10, bbox_to_anchor=(1.1, 1.05))

        plt.savefig(
            os.path.join(
                self.outDir, 'cycle_correlation(perSample).png'), dpi=600,
            bbox_inches='tight')
        plt.close('all')

        self.dfSaveCount = save_dataframe(
            df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
            )

    def log10transform(self, args):

        df = read_dataframe(outDir=self.outDir)

        df[self.markers] += 0.00000000001
        df[self.markers] = np.log10(df[self.markers])

        self.dfSaveCount = save_dataframe(
            df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
            )

    def pruneOutliers(self, args):

        df = read_dataframe(outDir=self.outDir)

        # plot raw signal intensity histrograms for inspection
        hist_facet = (
            df[['sample', 'area'] + self.markers]
            .sample(frac=1.0)
            .melt(
                id_vars=['sample', 'area'],
                var_name='channel', value_name='signal')
            )

        # sort naturally by sample column
        hist_facet = hist_facet.reindex(
            index=order_by_index(
                hist_facet.index, index_natsorted(hist_facet['sample'])))

        g = sns.FacetGrid(
            hist_facet, row='sample', col='channel',
            height=0.5, aspect=1.0, sharex=True, sharey=False
            )

        g.map(
            lambda x, y, color: plt.scatter(
                x, y, ec='none', s=0.001, c='r'), 'signal', 'area')

        g.set_titles(
            col_template="{col_name}", row_template="{row_name}",
            fontweight='bold', size=1.5, pad=0.75)

        for ax in g.axes.flatten():
            # ax.get_xaxis().set_ticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.tick_params(axis='both', which='major', labelsize=0.2, pad=-5)
            ax.xaxis.label.set_size(2.0)
            ax.yaxis.label.set_size(2.0)
            ax.grid(linewidth=0.075)

        plt.subplots_adjust(
            left=0.01, bottom=0.01, right=0.99,
            top=0.99, hspace=0.4, wspace=0.4
            )

        plt.savefig(
            os.path.join(
                self.outDir, 'scatter_plots(raw).png'), dpi=1000)
        plt.close()

        subprocess.call(
            ['open', '-a', 'Preview', os.path.join(
                self.outDir, 'scatter_plots(raw).png')])

        def submit(text, df):

            lowerPercentileCutoff = float(text.split(',')[0])
            upperPercentileCutoff = float(text.split(',')[1][1:])

            indices_to_drop = []
            for s in df['sample'].unique():
                for col in df[self.markers]:

                    data = df[col][df['sample'] == s]

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

            # rescale signal intensities
            for s in df['sample'].unique():
                data = df[df['sample'] == s]
                scaler = MinMaxScaler(
                        feature_range=(0, 1), copy=True).fit(
                            data[self.markers].values)
                rescaled_channel_data = scaler.transform(
                    data[self.markers].values)

                rescaled_channel_data = pd.DataFrame(
                    data=rescaled_channel_data,
                    index=data[self.markers].index,
                    columns=data[self.markers].columns
                    )

                df.update(rescaled_channel_data)

            # plot pruned and rescaled signal intensity histrograms
            hist_facet = (
                df[['sample', 'area'] + self.markers]
                .sample(frac=1.0)
                .melt(
                    id_vars=['sample', 'area'],
                    var_name='channel', value_name='signal')
                )

            # sort naturally by sample column
            hist_facet = hist_facet.reindex(
                index=order_by_index(
                    hist_facet.index,
                    index_natsorted(hist_facet['sample']))
                    )

            g = sns.FacetGrid(
                hist_facet, row='sample', col='channel',
                height=0.5, aspect=1.0, sharex=True, sharey=False
                )

            g.map(
                lambda x, y, color: plt.scatter(
                    x, y, ec='none', s=0.001, c='r'), 'signal', 'area')

            g.set_titles(
                col_template='{col_name}', row_template='{row_name}',
                fontweight='bold', size=1.5, pad=0.75)

            for ax in g.axes.flatten():
                # ax.get_xaxis().set_ticks([0, 0.25, 0.5, 0.75, 1.0])
                ax.tick_params(
                    axis='both', which='major', labelsize=0.2, pad=-5)
                ax.xaxis.label.set_size(2.0)
                ax.yaxis.label.set_size(2.0)
                ax.grid(linewidth=0.075)

            plt.subplots_adjust(
                left=0.01, bottom=0.01, right=0.99,
                top=0.99, hspace=0.4, wspace=0.4
                )

            plt.savefig(
                os.path.join(
                    self.outDir,
                    'scatterplots(pruned_rescaled).png'), dpi=1000)
            plt.close()

            subprocess.call(
                ['open', '-a', 'Preview',
                 os.path.join(
                    self.outDir, 'scatterplots(pruned_rescaled).png')]
                 )

            self.dfSaveCount = save_dataframe(
                df=df, outDir=self.outDir, dfSaveCount=self.dfSaveCount
                )

        plt.rcParams['figure.figsize'] = (7, 3)
        axbox = plt.axes([0.4, 0.525, 0.35, 0.15])
        text_box = TextBox(
            axbox,
            'lowerPercentileCutoff, upperPercentileCutoff',
            initial='',
            color='0.95',
            hovercolor='1.0',
            label_pad=0.05
            )
        text_box.label.set_size(10)

        text_box.on_submit(lambda text: submit(text, df))

        plt.show(block=True)

    def getOmeroImages(self, args):

        # get config_template.yml from cluster
        subprocess.call(['rsync', '-avP', self.cycleConfig, self.outDir])

        # load configuration file
        config = yaml.load(
            open(f'{self.outDir}/config_template.yml'))

        # set cycle map from configuration file to a variable
        cycle_map = config['cycle_map']

        df = read_dataframe(outDir=self.outDir)

        # make directory to store zarr arrays
        zarrs_dir = os.path.join(self.outDir, 'ashlar_zarrs')
        if not os.path.exists(zarrs_dir):
            os.mkdir(zarrs_dir)

            # initialize a dictionary to index and access zarr arrays
            zs = {}

            # select data compression algorithm for zarr arrays
            compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE)

            # loop over tissue samples
            for s in df['sample'].unique():

                # read segmentation mask for sample to crop
                # ashlar_output images to the same size
                segmentation = io.imread(
                    os.path.join(self.inDir, f'{s}_segmentation.tif'))

                # read DNA ashlar_output image
                img = io.imread(
                    os.path.join(
                        self.inDir, f'{s}_dna.tif'))

                # crop DNA image to size of U-Net segmentation mask
                img = img[0:segmentation.shape[0], 0:segmentation.shape[1]]

                # initialize zarr array for DNA image, append to a dictionary
                zs[f'{s}_dna'] = zarr.open(
                    f'{zarrs_dir}/{s}_dna.zarr', mode='w',
                    shape=(
                        img.shape[0], img.shape[1]),
                    chunks=(1000, 1000), dtype='uint16', compressor=compressor
                    )

                # update zarr array with DNA image data
                zs[f'{s}_dna'][:] = img

                # update zarr array with all antibody images
                for k, v in cycle_map.items():

                    print(f'Storing sample {s} {v}.tif')

                    # split key into cycle_number and channel_number components
                    cycle_number = str(int(k.split('_')[0])-1)
                    channel_number = k.split('_')[1]

                    # read ashlar_output image
                    img = io.imread(
                        os.path.join(
                            self.inDir, f'{s}_cycle_{cycle_number}_channel_'
                            f'{channel_number}.tif'))

                    # crop image to size of U-Net segmentation mask
                    img = img[0:segmentation.shape[0], 0:segmentation.shape[1]]

                    # initialize zarr array for image, append to a dictionary
                    zs[f'{s}_{v}'] = zarr.open(
                        f'{zarrs_dir}/{s}_{v}.zarr', mode='w',
                        shape=(
                            img.shape[0], img.shape[1]),
                        chunks=(1000, 1000), dtype='uint16',
                        compressor=compressor
                        )

                    # update zarr array with image data
                    zs[f'{s}_{v}'][:] = img
            print()

            # apply omero bottom and top signal intensity cutoffs
            subprocess.call(['rsync', '-avP', self.omeroSettings, self.outDir])

            channel_settings = yaml.load(
                open(f'{self.outDir}/sample_2_settings.yml'))

            # set channel start to 0.0 to see all background,
            # or use: channel_settings['channels'][channel]['start']
            settings_dict = {}
            for channel in channel_settings['channels']:
                settings_dict[
                    channel_settings['channels'][channel]['label']] = ((
                        channel_settings['channels'][channel]['start'],
                        channel_settings['channels'][channel]['end']),
                        channel_settings['channels'][channel]['color'])

            for k, v in zs.items():

                print(f'Applying OMERO signal intensity cutoffs for {k}')

                param_map_key = k.split('_')[1]

                bottom_omero_cutoff = settings_dict[param_map_key][0][0]
                top_omero_cutoff = settings_dict[param_map_key][0][1]

                temp = img_as_float(zs[k])
                temp -= (bottom_omero_cutoff/65535)
                temp /= ((top_omero_cutoff/65535)-(bottom_omero_cutoff/65535))
                temp = np.clip(temp, 0, 1)

                zs[k][:] = img_as_uint(temp)

        zs = {}
        for s in df['sample'].unique():
            zs[f'{s}_dna'] = zarr.open(f'{zarrs_dir}/{s}_dna.zarr', mode='r')
            for k, v in cycle_map.items():
                zs[f'{s}_{v}'] = zarr.open(
                    f'{zarrs_dir}/{s}_{v}.zarr', mode='r'
                    )

    def performPCA(self, args):

        df = read_dataframe(outDir=self.outDir)

        medians = df.groupby(
            ['sample']).median()[self.markers]

        # set grid background style
        sns.set_style('whitegrid')

        # specify PCA parameters
        pca = PCA(self.numPCAComponents)

        idx = medians.index

        # normalize signal intensities across samples (axis=0)
        if self.normalize is True:
            df = norm(medians, norm='l2', axis=0, copy=True, return_norm=False)
        else:
            medians = medians.values

        # apply PCA parameters to data
        projected = pca.fit_transform(medians)

        # generate dataframe for plot input
        scatter_input = pd.DataFrame(data=projected, index=idx)
        scatter_input.rename(columns={0: 'PC1', 1: 'PC2'}, inplace=True)

        # plot scores plot for first 2 PCs
        sns.scatterplot(
            data=scatter_input, x='PC1', y='PC2',
            hue=scatter_input.index,
            palette=[
                self.condHueDict['cd'] if 'cd' in i else
                self.condHueDict['hfd'] for i in medians.index],
            edgecolor='k', s=125, alpha=1.0, legend=False)

        # annotate data points
        if self.labelPoints is True:
            for label, x, y in zip(
              scatter_input.index,
              scatter_input['PC1'], scatter_input['PC2']):

                plt.annotate(
                    label, xy=(x, y), xytext=(3, 5), size=7.0,
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

    def performTSNE(self, args):

        df = read_dataframe(outDir=self.outDir)

        if not os.path.join(self.outDir, 'embedding.npy'):
            startTime = datetime.now()
            embedded = TSNE(
                n_components=self.numTSNEComponents,
                init='pca',
                perplexity=self.perplexity,
                early_exaggeration=self.earlyExaggeration,
                learning_rate=self.learningRate,
                metric=self.metric,
                random_state=5,
                n_jobs=-1).fit_transform(df[self.markers])
            print('Embedding completed in ' + str(datetime.now() - startTime))

            np.save(os.path.join(self.outDir, 'embedding'), embedded)
            df['emb1'] = embedded[:, 0]
            df['emb2'] = embedded[:, 1]

        else:
            embedded = np.load(os.path.join(self.outDir, 'embedding.npy'))
            df['emb1'] = embedded[:, 0]
            df['emb2'] = embedded[:, 1]

        sns.set_style('white')

        def submit(text, df):

            numerical_input = text.split('.')[0].strip()
            tup = tuple(map(int, numerical_input.split('-')))

            if len(tup) == 1:

                mylist = [tup[0]]

                scatter_point_size = 0.045

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
                    for color_by in ['cluster', 'sample']:

                        highlight = 'none'

                        if color_by == 'cluster':

                            num_colors_required = len(
                                df[color_by].unique()) - 1
                            numSubcatagories = math.ceil(
                                num_colors_required/10
                                )

                            # get cmap
                            cmap = categorical_cmap(
                                numCatagories=10,
                                numSubcatagories=numSubcatagories,
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
                                s=scatter_point_size,
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
                                    df=df, markers=self.markers,
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

                        elif color_by == 'sample':

                            num_colors_required = len(df[color_by].unique())
                            numSubcatagories = math.ceil(
                                num_colors_required/10
                                )

                            cmap = categorical_cmap(
                                numCatagories=10,
                                numSubcatagories=numSubcatagories,
                                cmap='tab10',
                                continuous=False
                                )

                            # trim qualitative cmap to number of unique samples
                            trim = (
                                len(cmap.colors) - len(df['sample'].unique())
                                )
                            cmap = ListedColormap(
                                cmap.colors[:-trim]
                                )

                            sample_dict = dict(
                                zip(
                                    natsorted(df['sample'].unique()),
                                    list(range(len(df['sample'].unique()))))
                                    )

                            c = [sample_dict[i] for i in df[color_by]]

                            ax2.scatter(
                                df['emb1'],
                                df['emb2'],
                                c=c,
                                cmap=cmap,
                                s=scatter_point_size,
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

                    if '.save' in text:

                        plt.savefig(
                            os.path.join(
                                self.outDir,
                                f'tsne_{color_by}.png'),
                            dpi=300
                            )

                        self.dfSaveCount = save_dataframe(
                            df=df, outDir=self.outDir,
                            dfSaveCount=self.dfSaveCount
                            )

                        plt.close('all')

                    plt.show(block=False)

            else:

                # df = df.sample(frac=0.01, random_state=22)

                mylist = list(range(tup[0], tup[1] + 1, 1))
                mylist.reverse()  # run higher sizes first for plot order

                scatter_point_size = 0.045

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

    def getClustermap(self, args):

        df = read_dataframe(outDir=self.outDir)

        clustermap_input = df[df['cluster'] != -1]

        cluster_heatmap_input = clustermap_input[
            self.markers + ['cluster']].groupby('cluster').mean()

        sns.set(font_scale=1.1)
        g = sns.clustermap(
            cluster_heatmap_input, cmap='viridis', standard_scale=1,
            square=False, yticklabels=1, linewidth=0.75, cbar=True
            )

        plt.savefig(
            os.path.join(
                self.outDir, f'clustermap.pdf'), bbox_inches='tight')

        plt.show(block=True)

    def lassoClusters(self, args):

        df = read_dataframe(outDir=self.outDir)

        lasso_dict = {}

        subplot_kw = dict(
            xlim=(df['emb1'].min(), df['emb1'].max()),
            ylim=(df['emb2'].min(), df['emb2'].max()),
            autoscale_on=False)

        subplot_kw = dict()

        fig, ax = plt.subplots(subplot_kw=subplot_kw)

        cmap = categorical_cmap(
            numCatagories=10, numSubcatagories=3,
            cmap='tab10', continuous=False
            )

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
        markers['cluster'].loc[~df.index.isin(idx)] = 100
        hi_markers = cluster_expression(
            df=markers, markers=self.markers, cluster=100, num_proteins=3
            )
        print(hi_markers)

    def cellDensities(self, args):

        df = read_dataframe(outDir=self.outDir)

        facet_input = df[df['cluster'] >= 0]

        facet_input = facet_input.groupby(
            ['condition', 'sample', 'cluster']).size().reset_index()
        facet_input.rename(columns={0: 'count'}, inplace=True)

        # divide cluster cell counts by total number of cells per sample
        for name, group in facet_input.groupby(['sample']):

            total_cells = group['count'].sum()

            facet_input.loc[group.index, 'cluster_conc'] = (
                group['count']/total_cells
                )

        # pad counts tables, some clusters may be absent from samples
        pad = pd.DataFrame()
        for cluster in sorted(facet_input['cluster'].unique()):
            to_append = pd.DataFrame(
                {'sample': natsorted(facet_input['sample'].unique()),
                 'cluster': [cluster]*len(facet_input['sample'].unique())})
            pad = pad.append(to_append)
        pad.reset_index(drop=True, inplace=True)
        pad['condition'] = [
            'cd' if 'cd' in i else 'hfd' for i in pad['sample']
            ]

        facet_input = facet_input.merge(
            pad, how='right', on=['condition', 'sample', 'cluster'])
        facet_input.fillna(value=0, inplace=True)

        means = {}
        for cluster in facet_input.groupby('cluster'):
            means[cluster[0]] = cluster[1]['cluster_conc'].median()
        facet_input['ave_cluster_conc'] = [
            means[i] for i in facet_input['cluster']
            ]

        facet_input.sort_values(
            by=['condition', 'cluster', 'cluster_conc'],
            ascending=[True, True, True], inplace=True)
        facet_input.reset_index(drop=True, inplace=True)

        sns.set(font_scale=0.4)

        g = sns.FacetGrid(
            data=facet_input, col='cluster', col_wrap=5,
            sharex=False, sharey=False, height=1.5, aspect=1.3
            )

        g.map(
            sns.barplot, 'sample', 'cluster_conc', 'condition',
            order=None, hue_order=None, ec='k', lw=0.25)

        plt.legend()

        [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]

        new_bar_width = 0.6
        for ax, title in zip(
            g.axes.flatten(), facet_input['cluster'].unique()
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
        plt.savefig(os.path.join(self.outDir, f'facetGrid.pdf'))
        plt.show(block=True)

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
                re.sub("\D", "", i) for i in group['sample']
                ]
            group['cluster'] = w

            # get denominator cell count of each sample
            if self.denominator_clusters is None:
                group['tissue_count'] = [
                    len(stats_input[stats_input['sample'] == i]) for
                    i in group['sample']]
            else:
                group['tissue_count'] = [
                    len(stats_input[(stats_input['sample'] == i) &
                        (stats_input[
                            'cluster'] == self.denominator_clusters)]) for
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

            if self.BonferroniCorrection:

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

    def curateFingernails(self, args):

        df = read_dataframe(outDir=self.outDir)

        zs = getZarrs(df=df, cycleConfig=self.cycleConfig, outDir=self.outDir)

        thumbnail_input = df[df['cluster'] >= 0]

        thumbnails_dir = os.path.join(self.outDir, 'thumbnails')
        if not os.path.exists(thumbnails_dir):
            os.mkdir(thumbnails_dir)

        for cluster in sorted(thumbnail_input['cluster'].unique()):

            markers_to_show = cluster_expression(
                df=thumbnail_input, markers=self.markers,
                cluster=cluster, num_proteins=3
                )

            markers_to_show = ['dna'] + markers_to_show

            color_dict = {}
            for i, j, k in zip(
              markers_to_show,

              [(0.5, 0.5, 0.5), (0.0, 1.0, 0.0),
               (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
              ['gray', 'green', 'red', 'blue']
              ):
                color_dict[i] = j

            print(cluster, markers_to_show)

            long_table = pd.DataFrame()

            for sample in thumbnail_input['sample'].unique():

                # initialize a zeros-array with the dimensions of the image
                overlay = np.zeros(
                    (zs[f'{sample}_dna'][:].shape[0],
                     zs[f'{sample}_dna'][:].shape[1]))

                # convert to rgb
                overlay = gray2rgb(overlay)

                # loop over the channels to create a dict of images
                for marker in markers_to_show:
                    print(f'{sample}_{marker}')
                    marker_img = img_as_float(
                        zs[f'{sample}_{marker}']
                        )
                    marker_img = gray2rgb(marker_img)
                    marker_img = (
                        marker_img * color_dict[marker]
                        )

                    # loop over channels to add each image to overlay seed
                    overlay += marker_img

                # crop out thumbnail images
                sample_cluster_subset = thumbnail_input[
                    (thumbnail_input['sample'] == sample)
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
                centroids = sample_cluster_subset[['x', 'y']]

                centroid_img = np.zeros(
                    (overlay.shape[0],
                     overlay.shape[1]))

                centroid_dist = 1  # in pixels
                for example, centroid in enumerate(centroids.iterrows()):

                    ystart_centroid = int(
                        centroid[1]['y'] - centroid_dist
                        )
                    ystop_centroid = int(
                        centroid[1]['y'] + centroid_dist
                        )

                    xstart_centroid = int(
                        centroid[1]['x'] - centroid_dist
                        )
                    xstop_centroid = int(
                        centroid[1]['x'] + centroid_dist
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
                        (centroid[1]['x'] == 0.0) &
                        (centroid[1]['y'] == 0.0)
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
                            centroid[1]['y'] - window_dist
                            )
                        ystop_window = int(
                            centroid[1]['y'] + window_dist
                            )

                        xstart_window = int(
                            centroid[1]['x'] - window_dist
                            )
                        xstop_window = int(
                            centroid[1]['x'] + window_dist
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
                bbox_to_anchor=(2.0, 30.0)
                )

            plt.savefig(
                os.path.join(
                    thumbnails_dir,
                    'cluster' + str(cluster) + '_thumbnails.pdf'),
                bbox_inches='tight')
            plt.close('all')
            del g
            gc.collect()

    def spatialAnalysis(self, args):

        df = read_dataframe(outDir=self.outDir)

        zs = getZarrs(df=df, cycleConfig=self.cycleConfig, outDir=self.outDir)

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
