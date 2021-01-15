import os
import yaml
import json

import fire

from qc_class import QC


def run(config_filepath: str):
    # load config file
    with open(config_filepath, 'r') as f:
        config = yaml.load(f, Loader=yaml.Loader)

    # parse param files
    markers_filepath = os.path.join(config['in_dir'], 'markers.csv')
    with open(markers_filepath, 'r') as f:
        markers_filepath = str(f.name)

    # create directories
    if not os.path.exists(config['out_dir']):
        os.makedirs(config['out_dir'])

    df_dir = os.path.join(config['out_dir'], 'dataframe_archive')
    if not os.path.exists(df_dir):
        os.makedirs(df_dir)

    modules = [
     'getSingleCellData',
     'setContrast',
     'selectROIs',
     'dnaIntensityCutoff',
     'dnaAreaCutoff',
     'crossCycleCorrelation',
     'log10transform',
     'pruneOutliers',
     'performPCA',
     'performTSNE',
     'getClustermap',
     'lassoClusters',
     'curateThumbnails',
     # 'makeZarrs',
     # 'cellDensities',
     # 'frequencyStats',
     # 'clusterBoxplots',
     # 'spatialAnalysis',
     ]

    # make instance of the QC class
    qc = QC(
        in_dir=config['in_dir'],
        out_dir=config['out_dir'],
        random_sample_size=config['random_sample_size'],
        mask_object=config['mask_object'],
        start_module=config['start_module'],
        sample_metadata=config['sample_metadata'],
        samples_to_exclude=config['samples_to_exclude'],
        markers_to_exclude=config['markers_to_exclude'],
        modules=modules,
        )

    start_idx = modules.index(config['start_module'])
    for i in modules[start_idx:]:
        print(f'Running: {i}')
        result = getattr(qc, i)
        result(config)


if __name__ == '__main__':
    fire.Fire(run)
