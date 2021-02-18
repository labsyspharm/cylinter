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
    #  'performTSNE',
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
        log_ratio_rnge=config.log_ratio_rnge,
        hexbins=config.hexbins,
        hexbin_grid_size=config.hexbin_grid_size,
        channelExclusionsPCA=config.channelExclusionsPCA,
        numPCAComponents=config.numPCAComponents,
        pointSize=config.pointSize,
        normalize=config.normalize,
        labelPoints=config.labelPoints,
        distanceCutoff=config.distanceCutoff,
        samplesToSilhouette=config.samplesToSilhouette,
        channelExclusionsTSNE=config.channelExclusionsTSNE,
        fracForEmbedding=config.fracForEmbedding,
        numTSNEComponents=config.numTSNEComponents,
        perplexity=config.perplexity,
        earlyExaggeration=config.earlyExaggeration,
        learningRate=config.learningRate,
        metric=config.metric,
        random_state=config.random_state,
        controlGroups=config.controlGroups,
        denominatorCluster=config.denominatorCluster,
        FDRCorrection=config.FDRCorrection,
        bonferroniCorrection=config.bonferroniCorrection,
        numThumbnails=config.numThumbnails,
        squareWindowDimension=config.squareWindowDimension,
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
