import pathlib
import yaml


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
        config.randomSampleSize = float(data['randomSampleSize'])
        config._parse_sample_metadata(data['sampleMetadata'])
        config.samplesToExclude = (data['samplesToExclude'])
        config.markersToExclude = (data['markersToExclude'])

        # CLASS MODULE CONFIGURATIONS
        config.viewSample = str(data['viewSample'])

        config.delintMode = bool(data['delintMode'])
        config.showAbChannels = bool(data['showAbChannels'])
        config.samplesForROISelection = list(data['samplesForROISelection'])

        config.yAxisGating = bool(data['yAxisGating'])

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
        config.fracForEmbeddingQC = float(data['fracForEmbeddingQC'])
        config.fracForEmbedding = float(data['fracForEmbedding'])
        config.dimensionEmbeddingQC = int(data['dimensionEmbeddingQC'])
        config.dimensionEmbedding = int(data['dimensionEmbedding'])

        config.perplexityQC = float(data['perplexityQC'])
        config.perplexity = float(data['perplexity'])
        config.earlyExaggerationQC = float(data['earlyExaggerationQC'])
        config.earlyExaggeration = float(data['earlyExaggeration'])
        config.learningRateTSNEQC = float(data['learningRateTSNEQC'])
        config.learningRateTSNE = float(data['learningRateTSNE'])
        config.metricQC = str(data['metricQC'])
        config.metric = str(data['metric'])
        config.randomStateQC = int(data['randomStateQC'])
        config.randomState = int(data['randomState'])

        config.nNeighborsQC = int(data['nNeighborsQC'])
        config.nNeighbors = int(data['nNeighbors'])
        config.learningRateUMAPQC = float(data['learningRateUMAPQC'])
        config.learningRateUMAP = float(data['learningRateUMAP'])
        config.minDistQC = float(data['minDistQC'])
        config.minDist = float(data['minDist'])
        config.repulsionStrengthQC = float(data['repulsionStrengthQC'])
        config.repulsionStrength = float(data['repulsionStrength'])

        config.controlGroups = list(data['controlGroups'])
        if (data['denominatorCluster']) is None:
            config.denominatorCluster = (data['denominatorCluster'])
        else:
            config.denominatorCluster = int(data['denominatorCluster'])
        config.FDRCorrection = bool(data['FDRCorrection'])

        config.numThumbnails = int(data['numThumbnails'])
        config.squareWindowDimension = int(data['squareWindowDimension'])
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

    @property
    def checkpoint_path(self):
        return self.outDir / 'checkpoints'

    def __repr__(self):
        kwargs_str = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"Config({kwargs_str})"
