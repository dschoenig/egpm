#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=0:30:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023 gcc/12.3 gdal/3.4.3 geos/3.12.0 python/3.11.5 udunits/2.2.28 r/4.2.1

Rscript 1_landscapes_1000_4cov_nl.R 1 100
