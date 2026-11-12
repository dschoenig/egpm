#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --array=1-200
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2026 gcc/15.3.0 gdal/3.13.3 geos/3.14.1 proj/9.9.0 udunits/2.2.28 python/3.15.2 r/4.6.1

Rscript 2_egp.R imbalance_high egp_som50 normal \
	0.01 50 1000 500 250 TRUE \
	$SLURM_CPUS_PER_TASK $SLURM_ARRAY_TASK_ID ${N_TASKS:-$SLURM_ARRAY_TASK_COUNT}
