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
        config.in_dir = pathlib.Path(data['in_dir']).resolve()
        config.out_dir = pathlib.Path(data['out_dir']).resolve()
        config.random_sample_size = float(data['random_sample_size'])
        config.mask_object = str(data['mask_object'])
        config._parse_sample_metadata(data['sample_metadata'])
        config.samples_to_exclude = (data['samples_to_exclude'])
        config.markers_to_exclude = (data['markers_to_exclude'])

        # CLASS MODULE CONFIGURATIONS
        config.view_sample = str(data['view_sample'])

        config.delint_mode = bool(data['delint_mode'])
        config.show_ab_channels = bool(data['show_ab_channels'])

        config.cutoffAxis = str(data['cutoffAxis'])
        if (data['log_ratio_rnge']) is None:
            config.log_ratio_rnge = (data['log_ratio_rnge'])
        else:
            config.log_ratio_rnge = tuple(data['log_ratio_rnge'])

        config.hexbins = bool(data['hexbins'])
        config.hexbin_grid_size = int(data['hexbin_grid_size'])

        config.channelExclusionsPCA = list(data['channelExclusionsPCA'])
        config.samplesToRemovePCA = list(data['samplesToRemovePCA'])
        config.dimensionPCA = int(data['dimensionPCA'])
        config.pointSize = float(data['pointSize'])
        config.normalize = bool(data['normalize'])
        config.labelPoints = bool(data['labelPoints'])
        config.distanceCutoff = float(data['distanceCutoff'])
        config.samplesToSilhouette = list(data['samplesToSilhouette'])

        config.embeddingAlgorithm1 = str(data['embeddingAlgorithm1'])
        config.embeddingAlgorithm2 = str(data['embeddingAlgorithm2'])
        config.channelExclusionsClustering1 = list(
            data['channelExclusionsClustering1']
            )
        config.channelExclusionsClustering2 = list(
            data['channelExclusionsClustering2']
            )
        config.samplesToRemoveClustering1 = list(
            data['samplesToRemoveClustering1']
            )
        config.samplesToRemoveClustering2 = list(
            data['samplesToRemoveClustering2']
            )
        config.normalizeTissueCounts1 = bool(data['normalizeTissueCounts1'])
        config.normalizeTissueCounts2 = bool(data['normalizeTissueCounts2'])
        config.fracForEmbedding1 = float(data['fracForEmbedding1'])
        config.fracForEmbedding2 = float(data['fracForEmbedding2'])
        config.dimensionEmbedding1 = int(data['dimensionEmbedding1'])
        config.dimensionEmbedding2 = int(data['dimensionEmbedding2'])

        config.perplexity1 = float(data['perplexity1'])
        config.perplexity2 = float(data['perplexity2'])
        config.earlyExaggeration1 = float(data['earlyExaggeration1'])
        config.earlyExaggeration2 = float(data['earlyExaggeration2'])
        config.learningRateTSNE1 = float(data['learningRateTSNE1'])
        config.learningRateTSNE2 = float(data['learningRateTSNE2'])
        config.metric1 = str(data['metric1'])
        config.metric2 = str(data['metric2'])
        config.random_state1 = int(data['random_state1'])
        config.random_state2 = int(data['random_state2'])

        config.nNeighbors1 = int(data['nNeighbors1'])
        config.nNeighbors2 = int(data['nNeighbors2'])
        config.learningRateUMAP1 = float(data['learningRateUMAP1'])
        config.learningRateUMAP2 = float(data['learningRateUMAP2'])
        config.minDist1 = float(data['minDist1'])
        config.minDist2 = float(data['minDist2'])
        config.repulsionStrength1 = float(data['repulsionStrength1'])
        config.repulsionStrength2 = float(data['repulsionStrength2'])

        config.controlGroups = list(data['controlGroups'])
        if (data['denominatorCluster1']) is None:
            config.denominatorCluster1 = (data['denominatorCluster1'])
        else:
            config.denominatorCluster1 = int(data['denominatorCluster1'])
        if (data['denominatorCluster2']) is None:
            config.denominatorCluster2 = (data['denominatorCluster2'])
        else:
            config.denominatorCluster2 = int(data['denominatorCluster2'])
        config.FDRCorrection = bool(data['FDRCorrection'])

        config.numThumbnails = int(data['numThumbnails'])
        config.squareWindowDimension = int(data['squareWindowDimension'])

        config.clustersToDrop = list(data['clustersToDrop'])

        config.bonferroniCorrection = bool(data['bonferroniCorrection'])

        config.cropDict = dict(data['cropDict'])
        config.spatialDict1 = dict(data['spatialDict1'])
        config.spatialDict2 = dict(data['spatialDict2'])
        config.radiusRange = list(data['radiusRange'])

        return config

    def _parse_sample_metadata(self, value):
        self.sample_conditions = {}
        self.sample_abbreviations = {}
        self.sample_statuses = {}
        self.sample_replicates = {}

        if value is None:
            return

        for name, terms in value.items():

            condition = str(terms[0])
            abbreviation = str(terms[1])
            status = str(terms[2])
            replicate = int(terms[3])

            self.sample_conditions[name] = condition
            self.sample_abbreviations[name] = abbreviation
            self.sample_statuses[name] = status
            self.sample_replicates[name] = replicate

    @property
    def checkpoint_path(self):
        return self.out_dir / 'checkpoints'

    def __repr__(self):
        kwargs_str = ', '.join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"Config({kwargs_str})"
