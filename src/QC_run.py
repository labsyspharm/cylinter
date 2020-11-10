import os
import yaml
import json

import fire

from QC_class import QC


def run(config_filepath: str):
    # load config file
    with open(config_filepath, 'r') as f:
        config = yaml.load(f, Loader=yaml.Loader)

    # parse param files
    with open(config['markers_filepath'], 'r') as f:
        config['markers_filepath'] = str(f.name)

    # create directories
    if not os.path.exists(config['outDir']):
        os.makedirs(config['outDir'])

    df_dir = os.path.join(config['outDir'], 'dataframe_archive')
    if not os.path.exists(df_dir):
        os.makedirs(df_dir)

    modules = [
     'getSingleCellData',
     'setContrast',
     'selectROIs',
     'dnaIntensityCutoff',
     'nuclearAreaCutoff',
     'crossCycleCorrelation',
     'log10transform',
     'pruneOutliers',
     'performPCA',
     'performTSNE',
     'getClustermap',
     'lassoClusters',
     'curateThumbnails',
     'makeZarrs',
     # 'cellDensities',
     # 'frequencyStats',
     # 'clusterBoxplots',
     # 'spatialAnalysis',
     ]

    # make instance of the QC class
    qc = QC(
        inDir=config['inDir'],
        outDir=config['outDir'],
        configDir=config['configDir'],
        markers_filepath=config['markers_filepath'],
        mask_object=config['mask_object'],
        start_module=config['start_module'],
        sample_metadata=config['sample_metadata'],
        randomSampleSize=config['randomSampleSize'],
        modules=modules,
        )

    start_idx = modules.index(config['start_module'])
    for i in modules[start_idx:]:
        print(f'Running: {i}')
        result = getattr(qc, i)
        result(config)


if __name__ == '__main__':
    fire.Fire(run)
