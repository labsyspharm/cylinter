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
        config['markers'] = [line.strip() for line in f]

    with open(config['samples_filepath'], 'r') as f:
        config['samples'] = yaml.load(f, Loader=yaml.Loader)

    with open(config['replicates_filepath'], 'r') as f:
        config['replicates'] = yaml.load(f, Loader=yaml.Loader)

    # create directory
    if not os.path.exists(config['outDir']):
        os.makedirs(config['outDir'])

    df_dir = os.path.join(config['outDir'], 'dataframe_archive')
    if not os.path.exists(df_dir):
        os.makedirs(df_dir)

    # make instance of the QC class
    qc = QC(
        inDir=config['inDir'],
        outDir=config['outDir'],
        markers=config['markers'],
        samples=config['samples'],
        replicates=config['replicates'],
        cycleConfig=config['cycleConfig'],
        omeroSettings=config['omeroSettings'],
        )

    qc.getSingleCellData(config)
    qc.lassoROIs(config)
    qc.dnaIntensityCutoff(config)
    qc.nuclearAreaCutoff(config)
    qc.crossCyleCorrelation(config)
    qc.log10transform(config)
    qc.pruneOutliers(config)
    qc.getOmeroImages(config)
    qc.performPCA(config)
    qc.performTSNE(config)
    qc.getClustermap(config)
    qc.lassoClusters(config)
    qc.cellDensities(config)
    qc.frequencyStats(config)
    qc.clusterBoxplots(config)
    qc.curateFingernails(config)
    qc.spatialAnalysis(config)


if __name__ == '__main__':
    fire.Fire(run)
