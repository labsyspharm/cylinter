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
        sample_conditions=config.conditions,
        sample_abbreviations=config.abbreviations,
        sample_replicates=config.replicates,
        samples_to_exclude=config.samples_to_exclude,
        markers_to_exclude=config.markers_to_exclude,
        modules=components.pipeline_module_names,
        )

    # start_idx = module_order[start_index:]
    for module in components.pipeline_modules[start_index:]:
        print(f'Running: {module}')
        data = module(data, qc, config)  # getattr(qc, module)
        # data(config)
        save_checkpoint(data, config, module)
