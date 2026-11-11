#!/bin/bash
# Rerun block 3 of 6 (747 jobs): no-effect Tweedie/binary landscapes (resubmitted;
# their block 1 runs lacked the no-effect Gaussian input), no-effect Gaussian models
# and summaries, Tweedie high models, marginals for the models of block 2.
set -euo pipefail
cd "$(dirname "$0")"
source rerun_lib.sh

scripts=(
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
prev_scripts=(
  cc_1_landscapes_noeff_low.sh
  cc_1_landscapes_noeff_high.sh
  cc_2_models_high_egp_som25_unequal.sh
  cc_2_models_tweedie_low_egp_som25.sh
  cc_2_models_tweedie_low_egp_som10.sh
  cc_2_models_tweedie_low_egp_som50.sh
  cc_2_models_tweedie_low_egp_som25_unequal.sh
  cc_3_summary_tweedie_low.sh
  cc_3_summary_tweedie_high.sh
  cc_3_summary_binary_low.sh
  cc_3_summary_binary_high.sh
  cc_3_marginal_low_egp_som25.sh
  cc_3_marginal_low_egp_som10.sh
  cc_3_marginal_low_egp_som50.sh
  cc_3_marginal_low_egp_som25_unequal.sh
  cc_3_marginal_high_egp_som25.sh
  cc_3_marginal_high_egp_som10.sh
  cc_3_marginal_high_egp_som50.sh
)
inter_files=(
  ../results/imbalance_high/egp_som25_unequal.inter.rds
  ../results/tweedie_imbalance_low/egp_som25.inter.rds
  ../results/tweedie_imbalance_low/egp_som10.inter.rds
  ../results/tweedie_imbalance_low/egp_som50.inter.rds
  ../results/tweedie_imbalance_low/egp_som25_unequal.inter.rds
)

submit_block
