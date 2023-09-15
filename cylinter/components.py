import logging
import functools

import matplotlib.pyplot as plt
import seaborn as sns

from cylinter.modules.aggregateData import aggregateData
from cylinter.modules.selectROIs import selectROIs
from cylinter.modules.intensityFilter import intensityFilter
from cylinter.modules.areaFilter import areaFilter
from cylinter.modules.cycleCorrelation import cycleCorrelation
from cylinter.modules.logTransform import logTransform
from cylinter.modules.pruneOutliers import pruneOutliers
from cylinter.modules.metaQC import metaQC
from cylinter.modules.PCA import PCA
from cylinter.modules.clustering import clustering
from cylinter.modules.clustermap import clustermap
from cylinter.modules.gating import gating
from cylinter.modules.setContrast import setContrast
from cylinter.modules.frequencyStats import frequencyStats
from cylinter.modules.curateThumbnails import curateThumbnails

logger = logging.getLogger(__name__)

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

# Pipeline module order, to be filled in by the @module decorator.
pipeline_modules = []
pipeline_module_names = []


def module(func):
    """
    Annotation for pipeline module functions.

    This function adds the given function to the registry list. It also wraps
    the given function to log a pre/post-call banner.

    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logger.info("=" * 70)
        logger.info("RUNNING MODULE: %s", func.__name__)
        result = func(*args, **kwargs)
        logger.info("=" * 70)
        logger.info("")
        return result
    pipeline_modules.append(wrapper)
    pipeline_module_names.append(wrapper.__name__)
    return wrapper


class QC(object):
    def __init__(self,

                 # config.yaml —
                 inDir=None,
                 outDir=None,
                 startModule=None,
                 sampleNames=None,
                 sampleConditions=None,
                 sampleConditionAbbrs=None,
                 sampleStatuses=None,
                 sampleReplicates=None,
                 samplesToExclude=None,
                 markersToExclude=None,

                 # selectROIs -
                 delintMode=None,
                 showAbChannels=None,
                 samplesForROISelection=None,

                 # intensityFilter -
                 numBinsIntensity=None,

                 # intensityArea -
                 numBinsArea=None,

                 # cycleCorrelation -
                 numBinsCorrelation=None,

                 # pruneOutliers -
                 hexbins=None,
                 hexbinGridSize=None,

                 # metaQC -
                 metaQC=None,
                 default_mcs=200,
                 default_reclass_tuple='0.75, 0.75',
                 embeddingAlgorithmQC=None,
                 channelExclusionsClusteringQC=None,
                 samplesToRemoveClusteringQC=None,
                 fracForEmbeddingQC=None,
                 dimensionEmbeddingQC=None,
                 topMarkersQC=None,
                 colormapAnnotationQC=None,
                 metricQC=None,
                 perplexityQC=None,
                 earlyExaggerationQC=None,
                 learningRateTSNEQC=None,

                 randomStateQC=None,
                 nNeighborsQC=None,
                 learningRateUMAPQC=None,
                 minDistQC=None,
                 repulsionStrengthQC=None,

                 # PCA module —
                 channelExclusionsPCA=None,
                 samplesToRemovePCA=None,
                 dimensionPCA=None,
                 pointSize=None,
                 labelPoints=None,
                 distanceCutoff=None,
                 conditionsToSilhouette=None,

                 # gating module —
                 gating=None,
                 channelExclusionsGating=None,
                 samplesToRemoveGating=None,
                 vectorThreshold=None,
                 classes=None,

                 # clustering module —
                 embeddingAlgorithm=None,
                 channelExclusionsClustering=None,
                 samplesToRemoveClustering=None,
                 normalizeTissueCounts=None,
                 fracForEmbedding=None,
                 dimensionEmbedding=None,
                 topMarkers=None,
                 colormapAnnotationClustering=None,
                 colormapAnnotation=None,
                 perplexity=None,
                 earlyExaggeration=None,
                 learningRateTSNE=None,
                 metric=None,
                 randomStateTSNE=None,
                 nNeighbors=None,
                 learningRateUMAP=None,
                 minDist=None,
                 repulsionStrength=None,
                 randomStateUMAP=None,

                 # frequencyStats —
                 controlGroups=None,
                 denominatorCluster=None,
                 FDRCorrection=None,

                 # curateThumbnails —
                 numThumbnails=None,
                 topMarkersThumbnails=None,
                 windowSize=None,
                 segOutlines=None,
                 ):

        assert topMarkers in ['channels', 'clusters'], \
            'Invalid input for topMarkers configuration parameter.'

        self.inDir = inDir
        self.outDir = outDir
        self.startModule = startModule
        self.sampleNames = sampleNames
        self.sampleConditions = sampleConditions
        self.sampleConditionAbbrs = sampleConditionAbbrs
        self.sampleStatuses = sampleStatuses
        self.sampleReplicates = sampleReplicates
        self.samplesToExclude = samplesToExclude
        self.markersToExclude = markersToExclude

        self.delintMode = delintMode
        self.showAbChannels = showAbChannels
        self.samplesForROISelection = samplesForROISelection

        self.numBinsIntensity = numBinsIntensity

        self.numBinsArea = numBinsArea

        self.numBinsCorrelation = numBinsCorrelation

        self.hexbins = hexbins
        self.hexbinGridSize = hexbinGridSize

        self.metaQC = metaQC
        self.default_mcsQC = default_mcs
        self.default_reclass_tuple = default_reclass_tuple
        self.embeddingAlgorithmQC = embeddingAlgorithmQC
        self.channelExclusionsClusteringQC = channelExclusionsClusteringQC
        self.samplesToRemoveClusteringQC = samplesToRemoveClusteringQC
        self.fracForEmbeddingQC = fracForEmbeddingQC
        self.dimensionEmbeddingQC = dimensionEmbeddingQC
        self.topMarkersQC = topMarkersQC
        self.colormapAnnotationQC = colormapAnnotationQC
        self.metricQC = metricQC
        self.perplexityQC = perplexityQC
        self.earlyExaggerationQC = earlyExaggerationQC
        self.learningRateTSNEQC = learningRateTSNEQC
        self.randomStateQC = randomStateQC
        self.nNeighborsQC = nNeighborsQC
        self.learningRateUMAPQC = learningRateUMAPQC
        self.minDistQC = minDistQC
        self.repulsionStrengthQC = repulsionStrengthQC

        self.channelExclusionsPCA = channelExclusionsPCA
        self.samplesToRemovePCA = samplesToRemovePCA
        self.dimensionPCA = dimensionPCA
        self.pointSize = pointSize
        self.labelPoints = labelPoints
        self.distanceCutoff = distanceCutoff
        self.conditionsToSilhouette = conditionsToSilhouette

        self.gating = gating
        self.channelExclusionsGating = channelExclusionsGating
        self.samplesToRemoveGating = samplesToRemoveGating
        self.vectorThreshold = vectorThreshold
        self.classes = classes

        self.embeddingAlgorithm = embeddingAlgorithm
        self.channelExclusionsClustering = channelExclusionsClustering
        self.samplesToRemoveClustering = samplesToRemoveClustering
        self.normalizeTissueCounts = normalizeTissueCounts
        self.fracForEmbedding = fracForEmbedding
        self.dimensionEmbedding = dimensionEmbedding
        self.topMarkers = topMarkers
        self.colormapAnnotationClustering = colormapAnnotationClustering
        self.perplexity = perplexity
        self.earlyExaggeration = earlyExaggeration
        self.learningRateTSNE = learningRateTSNE
        self.metric = metric
        self.randomStateTSNE = randomStateTSNE
        self.nNeighbors = nNeighbors
        self.learningRateUMAP = learningRateUMAP
        self.minDist = minDist
        self.repulsionStrength = repulsionStrength
        self.randomStateUMAP = randomStateUMAP

        self.controlGroups = controlGroups
        self.denominatorCluster = denominatorCluster
        self.FDRCorrection = FDRCorrection

        self.numThumbnails = numThumbnails
        self.topMarkersThumbnails = topMarkersThumbnails
        self.windowSize = windowSize
        self.segOutlines = segOutlines

    module(aggregateData)
    module(selectROIs)
    module(intensityFilter)
    module(areaFilter)
    module(cycleCorrelation)
    module(logTransform)
    module(pruneOutliers)
    module(metaQC)
    module(PCA)
    module(setContrast)
    module(gating)
    module(clustering)
    module(clustermap)
    module(frequencyStats)
    module(curateThumbnails)
