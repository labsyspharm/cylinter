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
        if (data['log_ratio_rnge']) is None:
            config.log_ratio_rnge = (data['log_ratio_rnge'])
        else:
            config.log_ratio_rnge = tuple(data['log_ratio_rnge'])
        config.hexbins = bool(data['hexbins'])
        config.hexbin_grid_size = int(data['hexbin_grid_size'])
        config.channelExclusionsPCA = list(data['channelExclusionsPCA'])
        config.numPCAComponents = int(data['numPCAComponents'])
        config.pointSize = float(data['pointSize'])
        config.normalize = bool(data['normalize'])
        config.labelPoints = bool(data['labelPoints'])
        config.distanceCutoff = float(data['distanceCutoff'])
        config.samplesToSilhouette = list(data['samplesToSilhouette'])
        config.channelExclusionsTSNE = list(data['channelExclusionsTSNE'])
        config.fracForEmbedding = float(data['fracForEmbedding'])
        config.numTSNEComponents = int(data['numTSNEComponents'])
        config.perplexity = float(data['perplexity'])
        config.earlyExaggeration = float(data['earlyExaggeration'])
        config.learningRate = float(data['learningRate'])
        config.metric = str(data['metric'])
        config.random_state = int(data['random_state'])
        config.numThumbnails = int(data['numThumbnails'])
        config.squareWindowDimension = int(data['squareWindowDimension'])
        config.controlGroups = list(data['controlGroups'])
        if (data['denominatorCluster']) is None:
            config.denominatorCluster = (data['denominatorCluster'])
        else:
            config.denominatorCluster = int(data['denominatorCluster'])
        config.FDRCorrection = bool(data['FDRCorrection'])
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
