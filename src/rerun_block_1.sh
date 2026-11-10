#!/bin/bash
# Rerun block 1 (see rerun.txt): Tweedie/binary landscapes and Gaussian models.
# Needs only the existing Gaussian landscapes (imbalance_* and noeff_imbalance_*).
set -euo pipefail
cd "$(dirname "$0")"

scripts=(
  cc_1_landscapes_low_tweedie.sh
  cc_1_landscapes_high_tweedie.sh
  cc_1_landscapes_low_binary.sh
  cc_1_landscapes_high_binary.sh
  cc_1_landscapes_noeff_low_tweedie.sh
  cc_1_landscapes_noeff_high_tweedie.sh
  cc_1_landscapes_noeff_low_binary.sh
  cc_1_landscapes_noeff_high_binary.sh
  cc_2_models_low_egp_som25.sh
  cc_2_models_low_egp_som10.sh
  cc_2_models_low_egp_som50.sh
  cc_2_models_low_egp_som25_unequal.sh
  cc_2_models_high_egp_som25.sh
  cc_2_models_high_egp_som10.sh
  cc_2_models_high_egp_som50.sh
  cc_2_models_high_egp_som25_unequal.sh
  cc_2_models_noeff_low_egp_som25.sh
  cc_2_models_noeff_high_egp_som25.sh
)

for s in "${scripts[@]}"; do
  [ -f "$s" ] || { echo "Missing script: $s" >&2; exit 1; }
done

for s in "${scripts[@]}"; do
  echo "Submitting $s"
  sbatch "$s"
done
