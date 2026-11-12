#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-10
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2026 gcc/15.3.0 gdal/3.13.3 geos/3.14.1 proj/9.9.0 udunits/2.2.28 python/3.15.2 r/4.6.1

Rscript 1_landscapes_binary.R noeff_imbalance_high noeff_binary_imbalance_high $SLURM_ARRAY_TASK_ID ${N_TASKS:-$SLURM_ARRAY_TASK_COUNT}
