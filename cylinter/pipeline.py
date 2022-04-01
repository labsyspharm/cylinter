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
        inDir=config.inDir,
        outDir=config.outDir,
        sampleNames=config.sampleNames,
        sampleConditions=config.sampleConditions,
        sampleConditionAbbrs=config.sampleConditionAbbrs,
        sampleStatuses=config.sampleStatuses,
        sampleReplicates=config.sampleReplicates,
        samplesToExclude=config.samplesToExclude,
        markersToExclude=config.markersToExclude,

        delintMode=config.delintMode,
        showAbChannels=config.showAbChannels,
        samplesForROISelection=config.samplesForROISelection,

        yAxisGating=config.yAxisGating,

        hexbins=config.hexbins,
        hexbinGridSize=config.hexbinGridSize,

        metaQC=config.metaQC,

        channelExclusionsPCA=config.channelExclusionsPCA,
        samplesToRemovePCA=config.samplesToRemovePCA,
        dimensionPCA=config.dimensionPCA,
        pointSize=config.pointSize,
        labelPoints=config.labelPoints,
        distanceCutoff=config.distanceCutoff,
        conditionsToSilhouette=config.conditionsToSilhouette,

        embeddingAlgorithmQC=config.embeddingAlgorithmQC,
        embeddingAlgorithm=config.embeddingAlgorithm,
        channelExclusionsClusteringQC=config.channelExclusionsClusteringQC,
        channelExclusionsClustering=config.channelExclusionsClustering,
        normalizeTissueCounts=config.normalizeTissueCounts,
        samplesToRemoveClusteringQC=config.samplesToRemoveClusteringQC,
        samplesToRemoveClustering=config.samplesToRemoveClustering,
        fracForEmbeddingQC=config.fracForEmbeddingQC,
        fracForEmbedding=config.fracForEmbedding,
        dimensionEmbeddingQC=config.dimensionEmbeddingQC,
        dimensionEmbedding=config.dimensionEmbedding,
        topMarkersQC=config.topMarkersQC,
        topMarkers=config.topMarkers,
        colormapChannel=config.colormapChannel,

        perplexityQC=config.perplexityQC,
        perplexity=config.perplexity,
        earlyExaggerationQC=config.earlyExaggerationQC,
        earlyExaggeration=config.earlyExaggeration,
        learningRateTSNEQC=config.learningRateTSNEQC,
        learningRateTSNE=config.learningRateTSNE,
        metricQC=config.metricQC,
        metric=config.metric,
        randomStateQC=config.randomStateQC,
        randomStateTSNE=config.randomStateTSNE,

        nNeighborsQC=config.nNeighborsQC,
        nNeighbors=config.nNeighbors,
        learningRateUMAPQC=config.learningRateUMAPQC,
        learningRateUMAP=config.learningRateUMAP,
        minDistQC=config.minDistQC,
        minDist=config.minDist,
        repulsionStrengthQC=config.repulsionStrengthQC,
        repulsionStrength=config.repulsionStrength,
        randomStateUMAP=config.randomStateUMAP,

        controlGroups=config.controlGroups,
        denominatorCluster=config.denominatorCluster,
        FDRCorrection=config.FDRCorrection,

        viewSample=config.viewSample,

        numThumbnails=config.numThumbnails,
        topMarkersThumbnails=config.topMarkersThumbnails,
        windowSize=config.windowSize,
        segOutlines=config.segOutlines,
        )

    # start_idx = module_order[start_index:]
    for module in components.pipeline_modules[start_index:]:
        print(f'Running: {module}')
        data = module(data, qc, config)  # getattr(qc, module)
        # data(config)
        save_checkpoint(data, config, module)
