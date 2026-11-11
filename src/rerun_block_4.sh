#!/bin/bash
# Rerun block 4 of 6 (910 jobs): binary low models, no-effect Tweedie/binary models
# and summaries, marginals for the models of block 3.
set -euo pipefail
cd "$(dirname "$0")"
source rerun_lib.sh

scripts=(
  cc_2_models_binary_low_egp_som25.sh
  cc_2_models_binary_low_egp_som10.sh
  cc_2_models_binary_low_egp_som50.sh
  cc_2_models_binary_low_egp_som25_unequal.sh
  cc_2_models_noeff_tweedie_low_egp_som25.sh
  cc_2_models_noeff_tweedie_high_egp_som25.sh
  cc_2_models_noeff_binary_low_egp_som25.sh
  cc_2_models_noeff_binary_high_egp_som25.sh
  cc_3_summary_noeff_tweedie_low.sh
  cc_3_summary_noeff_tweedie_high.sh
  cc_3_summary_noeff_binary_low.sh
  cc_3_summary_noeff_binary_high.sh
  cc_3_marginal_noeff_low_egp_som25.sh
  cc_3_marginal_noeff_high_egp_som25.sh
  cc_3_marginal_tweedie_high_egp_som25.sh
  cc_3_marginal_tweedie_high_egp_som10.sh
  cc_3_marginal_tweedie_high_egp_som50.sh
  cc_3_marginal_tweedie_high_egp_som25_unequal.sh
)
prev_scripts=(
  cc_1_landscapes_noeff_low_tweedie.sh
  cc_1_landscapes_noeff_high_tweedie.sh
  cc_1_landscapes_noeff_low_binary.sh
  cc_1_landscapes_noeff_high_binary.sh
  cc_3_summary_noeff_low.sh
  cc_3_summary_noeff_high.sh
  cc_2_models_noeff_low_egp_som25.sh
  cc_2_models_noeff_high_egp_som25.sh
  cc_2_models_tweedie_high_egp_som25.sh
  cc_2_models_tweedie_high_egp_som10.sh
  cc_2_models_tweedie_high_egp_som50.sh
  cc_2_models_tweedie_high_egp_som25_unequal.sh
  cc_3_marginal_high_egp_som25_unequal.sh
  cc_3_marginal_tweedie_low_egp_som25.sh
  cc_3_marginal_tweedie_low_egp_som10.sh
  cc_3_marginal_tweedie_low_egp_som50.sh
  cc_3_marginal_tweedie_low_egp_som25_unequal.sh
)
inter_files=(
  ../results/noeff_imbalance_low/egp_som25.inter.rds
  ../results/noeff_imbalance_high/egp_som25.inter.rds
  ../results/tweedie_imbalance_high/egp_som25.inter.rds
  ../results/tweedie_imbalance_high/egp_som10.inter.rds
  ../results/tweedie_imbalance_high/egp_som50.inter.rds
  ../results/tweedie_imbalance_high/egp_som25_unequal.inter.rds
)

submit_block
