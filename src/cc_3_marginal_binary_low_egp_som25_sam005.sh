#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 gdal/3.5.1 geos/3.10.2 python/3.10 udunits/2.2.28 r/4.2.2

Rscript 3_marginal_egp.R binary_imbalance_low egp_som25_sam005