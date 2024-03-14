import pathlib
import yaml
from dataclasses import dataclass


@dataclass(frozen=True)
class BooleanTerm:
    name: str
    negated: bool

    @classmethod
    def parse_str(cls, s):
        if s.startswith('+'):
            negated = False
            name = s[1:]
        elif s.startswith('-'):
            negated = True
            name = s[1:]
        else:
            negated = None
            name = s
        return cls(name, negated)

    def __repr__(self):
        s = self.name
        if self.negated:
            s = '~' + self.name
        return s

    def __invert__(self):
        return BooleanTerm(self.name, ~self.negated)


class Config:

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @classmethod
    def from_path(cls, path):
        config = cls()
        with open(path) as f:
            data = yaml.safe_load(f)
        config.inDir = pathlib.Path(data['inDir']).resolve()
        config.outDir = pathlib.Path(data['outDir']).resolve()
        config._parse_sample_metadata(data['sampleMetadata'])
        config.samplesToExclude = list(data['samplesToExclude'])
        config.counterstainChannel = str(data['counterstainChannel'])
        config.markersToExclude = list(data['markersToExclude'])

        # CLASS MODULE CONFIGURATIONS
        
        config.delintMode = bool(data['delintMode'])
        config.showAbChannels = bool(data['showAbChannels'])
        config.samplesForROISelection = list(data['samplesForROISelection'])
        config.autoArtifactDetection = bool(data['autoArtifactDetection'])
        config.artifactDetectionMethod = str(data['artifactDetectionMethod'])

        config.numBinsIntensity = int(data['numBinsIntensity'])

        config.numBinsArea = int(data['numBinsArea'])

        config.numBinsCorrelation = int(data['numBinsCorrelation'])

        config.hexbins = bool(data['hexbins'])
        config.hexbinGridSize = int(data['hexbinGridSize'])

        config.metaQC = bool(data['metaQC'])

        config.channelExclusionsPCA = list(data['channelExclusionsPCA'])
        config.samplesToRemovePCA = list(data['samplesToRemovePCA'])
        config.dimensionPCA = int(data['dimensionPCA'])
        config.pointSize = float(data['pointSize'])
        config.labelPoints = bool(data['labelPoints'])
        config.distanceCutoff = float(data['distanceCutoff'])
        config.conditionsToSilhouette = list(data['conditionsToSilhouette'])

        config.gating = bool(data['gating'])
        config.channelExclusionsGating = list(data['channelExclusionsGating'])
        config.samplesToRemoveGating = list(data['samplesToRemoveGating'])
        config.vectorThreshold = int(data['vectorThreshold'])
        config.vectorThreshold = int(data['vectorThreshold'])
        config._parse_classes(data['classes'])

        config.embeddingAlgorithmQC = str(data['embeddingAlgorithmQC'])
        config.embeddingAlgorithm = str(data['embeddingAlgorithm'])
        config.channelExclusionsClusteringQC = list(
            data['channelExclusionsClusteringQC']
        )
        config.channelExclusionsClustering = list(
            data['channelExclusionsClustering']
        )
        config.samplesToRemoveClusteringQC = list(
            data['samplesToRemoveClusteringQC']
        )
        config.samplesToRemoveClustering = list(
            data['samplesToRemoveClustering']
        )
        config.normalizeTissueCounts = bool(data['normalizeTissueCounts'])
        config.percentDataPerChunk = float(data['percentDataPerChunk'])
        config.fracForEmbedding = float(data['fracForEmbedding'])
        config.dimensionEmbedding = int(data['dimensionEmbedding'])
        config.colormapAnnotationQC = str(
            data['colormapAnnotationQC'])
        config.colormapAnnotationClustering = str(
            data['colormapAnnotationClustering'])

        config.perplexityQC = float(data['perplexityQC'])
        config.perplexity = float(data['perplexity'])
        config.earlyExaggerationQC = float(data['earlyExaggerationQC'])
        config.earlyExaggeration = float(data['earlyExaggeration'])
        config.learningRateTSNEQC = float(data['learningRateTSNEQC'])
        config.learningRateTSNE = float(data['learningRateTSNE'])
        config.metricQC = str(data['metricQC'])
        config.metric = str(data['metric'])
        config.randomStateQC = int(data['randomStateQC'])
        config.randomStateTSNE = int(data['randomStateTSNE'])

        config.nNeighborsQC = int(data['nNeighborsQC'])
        config.nNeighbors = int(data['nNeighbors'])
        config.learningRateUMAPQC = float(data['learningRateUMAPQC'])
        config.learningRateUMAP = float(data['learningRateUMAP'])
        config.minDistQC = float(data['minDistQC'])
        config.minDist = float(data['minDist'])
        config.repulsionStrengthQC = float(data['repulsionStrengthQC'])
        config.repulsionStrength = float(data['repulsionStrength'])
        config.randomStateUMAP = int(data['randomStateUMAP'])

        config.controlGroups = list(data['controlGroups'])
        if (data['denominatorCluster']) is None:
            config.denominatorCluster = (data['denominatorCluster'])
        else:
            config.denominatorCluster = int(data['denominatorCluster'])
        config.FDRCorrection = bool(data['FDRCorrection'])

        config.numThumbnails = int(data['numThumbnails'])
        config.windowSize = int(data['windowSize'])
        config.segOutlines = bool(data['segOutlines'])

        return config

    def _parse_sample_metadata(self, value):
        self.sampleNames = {}
        self.sampleConditions = {}
        self.sampleConditionAbbrs = {}
        self.sampleStatuses = {}
        self.sampleReplicates = {}

        if value is None:
            return

        for file_name, terms in value.items():

            name = str(terms[0])
            condition = str(terms[1])
            abbreviation = str(terms[2])
            status = str(terms[3])
            replicate = int(terms[4])

            self.sampleNames[file_name] = name
            self.sampleConditions[file_name] = condition
            self.sampleConditionAbbrs[file_name] = abbreviation
            self.sampleStatuses[file_name] = status
            self.sampleReplicates[file_name] = replicate

    def _parse_classes(self, value):

        self.classes = {}

        if value is None:
            return
        
        for outer_key, inner_dict in value.items():
            boo = [BooleanTerm.parse_str(t) for t in inner_dict['definition']]
            inner_dict['definition'] = boo
            self.classes[str(outer_key)] = inner_dict

    @property
    def checkpoint_path(self):
        return self.outDir / 'checkpoints'

    def __repr__(self):
        kwargs_str = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"Config({kwargs_str})"
