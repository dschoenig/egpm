#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --array=1-200
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.2

Rscript 2_egp.R tweedie_imbalance_low egp_som50 tweedie \
	0.01 50 1000 500 250 TRUE \
	$SLURM_CPUS_PER_TASK $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
