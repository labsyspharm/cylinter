import argparse
import os
import json

from classQC import QC

parser = argparse.ArgumentParser()
parser.add_argument('--inDir', type=str, default=None)
parser.add_argument('--outDir', type=str, default='QC_output')
parser.add_argument('--markers', nargs='*', type=str, default=[])
parser.add_argument('--samples', type=json.loads, default={})
parser.add_argument('--replicates', type=json.loads, default={})
parser.add_argument('--cycleConfig', type=str, default=None)
parser.add_argument('--omeroSettings', type=str, default=None)

args = parser.parse_args()

# assert(SOMETHING)  # placeholder for now

if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)

df_dir = os.path.join(args.outDir, 'dataframe_archive')
if not os.path.exists(df_dir):
    os.makedirs(df_dir)

# make instance of the QC class
qc = QC(
    inDir=args.inDir,
    outDir=args.outDir,
    markers=args.markers,
    samples=args.samples,
    replicates=args.replicates,
    cycleConfig=args.cycleConfig,
    omeroSettings=args.omeroSettings,
    )
qc.getSingleCellData(args)
qc.lassoROIs(args)
qc.dnaIntensityCutoff(args)
qc.nuclearAreaCutoff(args)
qc.crossCyleCorrelation(args)
qc.log10transform(args)
qc.pruneOutliers(args)
qc.getOmeroImages(args)
qc.performPCA(args)
qc.performTSNE(args)
qc.getClustermap(args)
qc.lassoClusters(args)
qc.cellDensities(args)
qc.frequencyStats(args)
qc.clusterBoxplots(args)
qc.curateFingernails(args)
qc.spatialAnalysis(args)
