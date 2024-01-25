#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --time=18:00:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=brm19

module load StdEnv/2023 gcc/12.3 gdal/3.7.2 geos/3.12.0 python/3.11.5 udunits/2.2.28 r/4.3.1

Rscript 5_compare_brm.R 19
