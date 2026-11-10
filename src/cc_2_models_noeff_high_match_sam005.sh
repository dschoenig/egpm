#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=14G
#SBATCH --time=3:00:00
#SBATCH --array=1-100
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2026 gcc/15.3.0 gdal/3.13.3 geos/3.14.1 proj/9.9.0 udunits/2.2.28 python/3.15.2 r/4.6.1

Rscript 2_match.R noeff_imbalance_high match_sam005 normal 0.005 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
