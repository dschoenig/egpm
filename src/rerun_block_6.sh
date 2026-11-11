#!/bin/bash
# Rerun block 6 of 6 (4 jobs): marginals for the binary high models of block 5.
set -euo pipefail
cd "$(dirname "$0")"
source rerun_lib.sh

scripts=(
  cc_3_marginal_binary_high_egp_som25.sh
  cc_3_marginal_binary_high_egp_som10.sh
  cc_3_marginal_binary_high_egp_som50.sh
  cc_3_marginal_binary_high_egp_som25_unequal.sh
)
prev_scripts=(
  cc_2_models_binary_high_egp_som25.sh
  cc_2_models_binary_high_egp_som10.sh
  cc_2_models_binary_high_egp_som50.sh
  cc_2_models_binary_high_egp_som25_unequal.sh
  cc_3_marginal_binary_low_egp_som25.sh
  cc_3_marginal_binary_low_egp_som10.sh
  cc_3_marginal_binary_low_egp_som50.sh
  cc_3_marginal_binary_low_egp_som25_unequal.sh
  cc_3_marginal_noeff_tweedie_low_egp_som25.sh
  cc_3_marginal_noeff_tweedie_high_egp_som25.sh
  cc_3_marginal_noeff_binary_low_egp_som25.sh
  cc_3_marginal_noeff_binary_high_egp_som25.sh
)
inter_files=(
  ../results/binary_imbalance_high/egp_som25.inter.rds
  ../results/binary_imbalance_high/egp_som10.inter.rds
  ../results/binary_imbalance_high/egp_som50.inter.rds
  ../results/binary_imbalance_high/egp_som25_unequal.inter.rds
)

submit_block
