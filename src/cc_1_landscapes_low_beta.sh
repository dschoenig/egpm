#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-10
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023 gcc/12.3 gdal/3.7.2 geos/3.12.0 python/3.11.5 udunits/2.2.28 r/4.3.1

Rscript 1_landscapes_beta.R imbalance_low beta_imbalance_low $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
