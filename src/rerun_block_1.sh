#!/bin/bash
# Rerun block 1 of 6 (980 jobs): Tweedie/binary landscapes, first Gaussian models.
# Needs only the existing Gaussian landscapes (imbalance_*). The four no-effect
# Tweedie/binary landscape scripts fail here (no noeff_imbalance_* yet) and are
# resubmitted in block 3.
set -euo pipefail
cd "$(dirname "$0")"
source rerun_lib.sh

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
)
prev_scripts=()
inter_files=()

submit_block
