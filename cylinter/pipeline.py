import pyarrow
import pyarrow.parquet
import pandas as pd
from . import components


def save_checkpoint(data, config, module):
    module_name = module.__name__
    path = config.checkpoint_path / f"{module_name}.parquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    # Ideally we would have just used pandas' to_parquet instead of calling
    # pyarrow directly, but to_parquet has as an over-zealous validity check on
    # the input dataframe that errors with a column MultiIndex. If that bug is
    # resolved we can switch to use just the following commented line.
    # data.to_parquet(path, index=True)
    table = pyarrow.Table.from_pandas(data)
    pyarrow.parquet.write_table(table, path)


def run_pipeline(config, start_module_name):
    if (
        start_module_name is None
        or start_module_name == components.pipeline_module_names[0]
    ):
        start_index = 0
        data = None
    else:
        start_index = components.pipeline_module_names.index(start_module_name)
        previous_module_name = components.pipeline_module_names[start_index - 1]
        checkpoint_file_path = (
            config.checkpoint_path / f"{previous_module_name}.parquet"
        )
        if not checkpoint_file_path.exists():
            raise Exception(
                f"Checkpoint file for module {previous_module_name} not found"
            )
        data = pd.read_parquet(checkpoint_file_path)

    # module_order = [
    #  'getSingleCellData',
    #  'setContrast',
    #  'selectROIs',
    #  'dnaIntensityCutoff',
    #  'dnaAreaCutoff',
    #  'crossCycleCorrelation',
    #  'log10transform',
    #  'pruneOutliers',
    #  'performPCA',
    #  'performClustering',
    #  'getClustermap',
    #  'lassoClusters',
    #  'curateThumbnails',
    #  ]

    # make instance of the QC class
    qc = components.QC(
        in_dir=config.in_dir,
        out_dir=config.out_dir,
        random_sample_size=config.random_sample_size,
        mask_object=config.mask_object,
        sample_conditions=config.sample_conditions,
        sample_abbreviations=config.sample_abbreviations,
        sample_statuses=config.sample_statuses,
        sample_replicates=config.sample_replicates,
        samples_to_exclude=config.samples_to_exclude,
        markers_to_exclude=config.markers_to_exclude,

        view_sample=config.view_sample,

        delint_mode=config.delint_mode,
        show_ab_channels=config.show_ab_channels,

        cutoffAxis=config.cutoffAxis,
        log_ratio_rnge=config.log_ratio_rnge,

        hexbins=config.hexbins,
        hexbin_grid_size=config.hexbin_grid_size,

        channelExclusionsPCA=config.channelExclusionsPCA,
        samplesToRemovePCA=config.samplesToRemovePCA,
        dimensionPCA=config.dimensionPCA,
        pointSize=config.pointSize,
        normalize=config.normalize,
        labelPoints=config.labelPoints,
        distanceCutoff=config.distanceCutoff,
        samplesToSilhouette=config.samplesToSilhouette,

        embeddingAlgorithm1=config.embeddingAlgorithm1,
        embeddingAlgorithm2=config.embeddingAlgorithm2,
        channelExclusionsClustering1=config.channelExclusionsClustering1,
        channelExclusionsClustering2=config.channelExclusionsClustering2,
        normalizeTissueCounts1=config.normalizeTissueCounts1,
        normalizeTissueCounts2=config.normalizeTissueCounts2,
        samplesToRemoveClustering1=config.samplesToRemoveClustering1,
        samplesToRemoveClustering2=config.samplesToRemoveClustering2,
        fracForEmbedding1=config.fracForEmbedding1,
        fracForEmbedding2=config.fracForEmbedding2,
        dimensionEmbedding1=config.dimensionEmbedding1,
        dimensionEmbedding2=config.dimensionEmbedding2,

        perplexity1=config.perplexity1,
        perplexity2=config.perplexity2,
        earlyExaggeration1=config.earlyExaggeration1,
        earlyExaggeration2=config.earlyExaggeration2,
        learningRateTSNE1=config.learningRateTSNE1,
        learningRateTSNE2=config.learningRateTSNE2,
        metric1=config.metric1,
        metric2=config.metric2,
        random_state1=config.random_state1,
        random_state2=config.random_state2,

        nNeighbors1=config.nNeighbors1,
        nNeighbors2=config.nNeighbors2,
        learningRateUMAP1=config.learningRateUMAP1,
        learningRateUMAP2=config.learningRateUMAP2,
        minDist1=config.minDist1,
        minDist2=config.minDist2,
        repulsionStrength1=config.repulsionStrength1,
        repulsionStrength2=config.repulsionStrength2,

        controlGroups=config.controlGroups,
        denominatorCluster1=config.denominatorCluster1,
        denominatorCluster2=config.denominatorCluster2,
        FDRCorrection=config.FDRCorrection,

        numThumbnails=config.numThumbnails,
        squareWindowDimension=config.squareWindowDimension,

        clustersToDrop=config.clustersToDrop,

        bonferroniCorrection=config.bonferroniCorrection,

        cropDict=config.cropDict,
        spatialDict1=config.spatialDict1,
        spatialDict2=config.spatialDict2,
        radiusRange=config.radiusRange,
        )

    # start_idx = module_order[start_index:]
    for module in components.pipeline_modules[start_index:]:
        print(f'Running: {module}')
        data = module(data, qc, config)  # getattr(qc, module)
        # data(config)
        save_checkpoint(data, config, module)
