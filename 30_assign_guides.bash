source env_vars.sh

mkdir -p $GUIDEBENDER_RES_DIR

# we'll run each nextflow job from a separate subfolder so we can
# manually run them in parallel if we want to
mkdir -p nextflow_rundirs

GUIDEBENDER_SIM_DIR=$(pwd -P $GUIDEBENDER_SIM_DIR)
GUIDEBENDER_RES_DIR=$(pwd -P $GUIDEBENDER_RES_DIR)

# test run
cd nextflow_rundirs
mkdir -p test
cd test
nextflow ../../main.nf \
         -resume \
         -profile slurm \
         --sample_sheet $GUIDEBENDER_SIM_DIR/test/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/test
cd ../..

# 2k cells, high-gRNA regime, varying number of guides
cd nextflow_rundirs
mkdir -p vary_nguides_2k_high_expr
cd vary_nguides_2k_high_expr
nextflow ../../main.nf \
         -resume \
         -profile slurm \
         --sample_sheet $GUIDEBENDER_SIM_DIR/numCells2000_medUmi100_snr4_endo75_varyNumGuides/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/numCells2000_medUmi100_snr4_endo75_varyNumGuides
cd ../..

# 2k cells, low-gRNA regime, varying number of guides
# NOTE: do NOT use -resume, as we use this run for timing.
cd nextflow_rundirs
mkdir -p vary_nguides_2k_low_expr
cd vary_nguides_2k_low_expr
nextflow ../../main.nf \
         -profile slurm \
         --sample_sheet $GUIDEBENDER_SIM_DIR/numCells2000_medUmi20_snr1_endo25_varyNumGuides/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/numCells2000_medUmi20_snr1_endo25_varyNumGuides
cd ../..

# 20k cells, high-gRNA regime, varying number of guides
# note demuxEM requires a lot of memory in this scenario
cd nextflow_rundirs
mkdir -p vary_nguides_20k_high_expr
cd vary_nguides_20k_high_expr
nextflow ../../main.nf \
         -resume \
         -profile slurm \
         --demuxemMemFactor 4 \
         --maxMemFactor 2 \
         --skipCrispatBig \
         --skipCleanser \
         --sample_sheet $GUIDEBENDER_SIM_DIR/numCells20k_medUmi100_snr4_endo75_varyNumGuides/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/numCells20k_medUmi100_snr4_endo75_varyNumGuides
cd ../..

# 20k cells, low-gRNA regime, varying number of guides
# NOTE: do NOT use -resume, as we use this run for timing.
cd nextflow_rundirs
mkdir -p vary_nguides_20k_low_expr
cd vary_nguides_20k_low_expr
nextflow ../../main.nf \
         -profile slurm \
         --maxMemFactor 2 \
         --skipCrispatBig \
         --skipCleanser \
         --sample_sheet $GUIDEBENDER_SIM_DIR/numCells20k_medUmi20_snr1_endo25_varyNumGuides/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/numCells20k_medUmi20_snr1_endo25_varyNumGuides
cd ../..

# varying MOI
cd nextflow_rundirs
mkdir -p vary_moi_2k
cd vary_moi_2k
nextflow ../../main.nf \
         -resume \
         -profile slurm \
         --sample_sheet $GUIDEBENDER_SIM_DIR/numCells2000_numGuides100_varyMOI/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/numCells2000_numGuides100_varyMOI
cd ../..

# barnyard run
#
# NOTE: A fake ground-truth was created so that the simulation
# pipeline would run to completion (as the nextflow pipeline computes
# the accuracy at the end). In particular, all entries with >0 counts
# were labeled as being truly present.These accuracy metrics should be
# ignored as it is just a hack to run the pipeline!  Instead, extract
# the full results and manually compare them to the true species after
# running this (this is done in 40_analysis.R).
#
# (TODO: refactor the pipeline so it can be run without requiring
# ground truth and computing accuracy metrics)
cd nextflow_rundirs
mkdir -p barnyard
cd barnyard
nextflow ../../main.nf \
         -resume \
         -profile slurm \
         --sample_sheet $GUIDEBENDER_SIM_DIR/cleanser_barnyard_data/sample_sheet.csv \
         --outdir $GUIDEBENDER_RES_DIR/cleanser_barnyard_data
cd ../..
