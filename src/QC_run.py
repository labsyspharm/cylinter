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

    # make instance of the QC class
    qc = QC(
        inDir=config['inDir'],
        outDir=config['outDir'],
        markers_filepath=config['markers_filepath'],
        mask_object=config['mask_object'],
        sample_metadata=config['sample_metadata'],
        randomSampleSize=config['randomSampleSize']
        )

    qc.getSingleCellData(config)
    qc.lassoROIs(config)
    qc.dnaIntensityCutoff(config)
    qc.nuclearAreaCutoff(config)
    qc.crossCyleCorrelation(config)
    qc.log10transform(config)
    qc.pruneOutliers(config)
    qc.makeZarrs(config)
    qc.performPCA(config)
    qc.performTSNE(config)
    qc.getClustermap(config)
    qc.lassoClusters(config)
    qc.cellDensities(config)
    # qc.frequencyStats(config)
    # qc.clusterBoxplots(config)
    qc.curateThumbnails(config)
    # qc.spatialAnalysis(config)


if __name__ == '__main__':
    fire.Fire(run)
