#!/bin/bash
#SBATCH --account=def-cricrime 
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-user=schonig.daniel@courrier.uqam.ca
#SBATCH --mail-type=ALL

cd ../landscapes

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_imbalance_high -g imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_imbalance_low -g imbalance_low

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_tweedie_imbalance_high -g tweedie_imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_tweedie_imbalance_low -g tweedie_imbalance_low

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_binary_imbalance_high -g binary_imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_binary_imbalance_low -g binary_imbalance_low

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_noeff_imbalance_high -g noeff_imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_noeff_imbalance_low -g noeff_imbalance_low

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_noeff_tweedie_imbalance_high -g noeff_tweedie_imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_noeff_tweedie_imbalance_low -g noeff_tweedie_imbalance_low

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_noeff_binary_imbalance_high -g noeff_binary_imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_noeff_binary_imbalance_low -g noeff_binary_imbalance_low

dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_beta_imbalance_high -g beta_imbalance_high
dar -w -c /home/schoed/nearline/def-cricrime/schoed/egp_landscapes_beta_imbalance_low -g beta_imbalance_low
