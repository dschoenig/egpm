# EGP 10 all landscapes
## Normal 
cp cc_2_models_high_egp_som10.sh cc_2_models_low_egp_som10.sh
sed -i "s/imbalance_high/imbalance_low/" cc_2_models_low_egp_som10.sh 
sed -i "s/normal/normal/" cc_2_models_low_egp_som10.sh 


## Binary
cp cc_2_models_high_egp_som10.sh cc_2_models_binary_low_egp_som10.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_2_models_binary_low_egp_som10.sh 
sed -i "s/normal/binary/" cc_2_models_binary_low_egp_som10.sh 
cp cc_2_models_high_egp_som10.sh cc_2_models_binary_high_egp_som10.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_2_models_binary_high_egp_som10.sh 
sed -i "s/normal/binary/" cc_2_models_binary_high_egp_som10.sh 

## Tweedie
cp cc_2_models_high_egp_som10.sh cc_2_models_tweedie_low_egp_som10.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_2_models_tweedie_low_egp_som10.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_low_egp_som10.sh 
cp cc_2_models_high_egp_som10.sh cc_2_models_tweedie_high_egp_som10.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_2_models_tweedie_high_egp_som10.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_high_egp_som10.sh 


sbatch cc_2_models_low_egp_som10.sh
sbatch cc_2_models_high_egp_som10.sh
sbatch cc_2_models_binary_low_egp_som10.sh
sbatch cc_2_models_binary_high_egp_som10.sh
sbatch cc_2_models_tweedie_low_egp_som10.sh
sbatch cc_2_models_tweedie_high_egp_som10.sh


# EGP 50 all landscapes
## Normal 
cp cc_2_models_high_egp_som50.sh cc_2_models_low_egp_som50.sh
sed -i "s/imbalance_high/imbalance_low/" cc_2_models_low_egp_som50.sh 
sed -i "s/normal/normal/" cc_2_models_low_egp_som50.sh 

## Binary
cp cc_2_models_high_egp_som50.sh cc_2_models_binary_low_egp_som50.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_2_models_binary_low_egp_som50.sh 
sed -i "s/normal/binary/" cc_2_models_binary_low_egp_som50.sh 
cp cc_2_models_high_egp_som50.sh cc_2_models_binary_high_egp_som50.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_2_models_binary_high_egp_som50.sh 
sed -i "s/normal/binary/" cc_2_models_binary_high_egp_som50.sh 

## Tweedie
cp cc_2_models_high_egp_som50.sh cc_2_models_tweedie_low_egp_som50.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_2_models_tweedie_low_egp_som50.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_low_egp_som50.sh 
cp cc_2_models_high_egp_som50.sh cc_2_models_tweedie_high_egp_som50.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_2_models_tweedie_high_egp_som50.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_high_egp_som50.sh 


sbatch cc_2_models_low_egp_som50.sh
sbatch cc_2_models_high_egp_som50.sh
sbatch cc_2_models_binary_low_egp_som50.sh
sbatch cc_2_models_binary_high_egp_som50.sh
sbatch cc_2_models_tweedie_low_egp_som50.sh
sbatch cc_2_models_tweedie_high_egp_som50.sh



# EGP 25 unequal all landscapes
## Normal 
cp cc_2_models_high_egp_som25_unequal.sh cc_2_models_low_egp_som25_unequal.sh
sed -i "s/imbalance_high/imbalance_low/" cc_2_models_low_egp_som25_unequal.sh 
sed -i "s/normal/normal/" cc_2_models_low_egp_som25_unequal.sh 

## Binary
cp cc_2_models_high_egp_som25_unequal.sh cc_2_models_binary_low_egp_som25_unequal.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_2_models_binary_low_egp_som25_unequal.sh 
sed -i "s/normal/binary/" cc_2_models_binary_low_egp_som25_unequal.sh 
cp cc_2_models_high_egp_som25_unequal.sh cc_2_models_binary_high_egp_som25_unequal.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_2_models_binary_high_egp_som25_unequal.sh 
sed -i "s/normal/binary/" cc_2_models_binary_high_egp_som25_unequal.sh 

## Tweedie
cp cc_2_models_high_egp_som25_unequal.sh cc_2_models_tweedie_low_egp_som25_unequal.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_2_models_tweedie_low_egp_som25_unequal.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_low_egp_som25_unequal.sh 
cp cc_2_models_high_egp_som25_unequal.sh cc_2_models_tweedie_high_egp_som25_unequal.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_2_models_tweedie_high_egp_som25_unequal.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_high_egp_som25_unequal.sh 


sbatch cc_2_models_low_egp_som25_unequal.sh
sbatch cc_2_models_high_egp_som25_unequal.sh
sbatch cc_2_models_binary_low_egp_som25_unequal.sh
sbatch cc_2_models_binary_high_egp_som25_unequal.sh
sbatch cc_2_models_tweedie_low_egp_som25_unequal.sh
sbatch cc_2_models_tweedie_high_egp_som25_unequal.sh




# EGP large sample all landscapes
## Normal 
cp cc_2_models_high_egp_som25_sam02.sh cc_2_models_low_egp_som25_sam02.sh
sed -i "s/imbalance_high/imbalance_low/" cc_2_models_low_egp_som25_sam02.sh 
sed -i "s/normal/normal/" cc_2_models_low_egp_som25_sam02.sh 

## Binary
cp cc_2_models_high_egp_som25_sam02.sh cc_2_models_binary_low_egp_som25_sam02.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_2_models_binary_low_egp_som25_sam02.sh 
sed -i "s/normal/binary/" cc_2_models_binary_low_egp_som25_sam02.sh 
cp cc_2_models_high_egp_som25_sam02.sh cc_2_models_binary_high_egp_som25_sam02.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_2_models_binary_high_egp_som25_sam02.sh 
sed -i "s/normal/binary/" cc_2_models_binary_high_egp_som25_sam02.sh 

## Tweedie
cp cc_2_models_high_egp_som25_sam02.sh cc_2_models_tweedie_low_egp_som25_sam02.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_2_models_tweedie_low_egp_som25_sam02.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_low_egp_som25_sam02.sh 
cp cc_2_models_high_egp_som25_sam02.sh cc_2_models_tweedie_high_egp_som25_sam02.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_2_models_tweedie_high_egp_som25_sam02.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_high_egp_som25_sam02.sh


sbatch cc_2_models_low_egp_som25_sam02.sh
sbatch cc_2_models_high_egp_som25_sam02.sh
sbatch cc_2_models_binary_low_egp_som25_sam02.sh
sbatch cc_2_models_binary_high_egp_som25_sam02.sh
sbatch cc_2_models_tweedie_low_egp_som25_sam02.sh
sbatch cc_2_models_tweedie_high_egp_som25_sam02.sh


# EGP small sample all landscapes
## Normal 
cp cc_2_models_high_egp_som25_sam005.sh cc_2_models_low_egp_som25_sam005.sh
sed -i "s/imbalance_high/imbalance_low/" cc_2_models_low_egp_som25_sam005.sh 
sed -i "s/normal/normal/" cc_2_models_low_egp_som25_sam005.sh 

## Binary
cp cc_2_models_high_egp_som25_sam005.sh cc_2_models_binary_low_egp_som25_sam005.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_2_models_binary_low_egp_som25_sam005.sh 
sed -i "s/normal/binary/" cc_2_models_binary_low_egp_som25_sam005.sh 
cp cc_2_models_high_egp_som25_sam005.sh cc_2_models_binary_high_egp_som25_sam005.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_2_models_binary_high_egp_som25_sam005.sh 
sed -i "s/normal/binary/" cc_2_models_binary_high_egp_som25_sam005.sh 

## Tweedie
cp cc_2_models_high_egp_som25_sam005.sh cc_2_models_tweedie_low_egp_som25_sam005.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_2_models_tweedie_low_egp_som25_sam005.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_low_egp_som25_sam005.sh 
cp cc_2_models_high_egp_som25_sam005.sh cc_2_models_tweedie_high_egp_som25_sam005.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_2_models_tweedie_high_egp_som25_sam005.sh 
sed -i "s/normal/tweedie/" cc_2_models_tweedie_high_egp_som25_sam005.sh 


sbatch cc_2_models_low_egp_som25_sam005.sh
sbatch cc_2_models_high_egp_som25_sam005.sh
sbatch cc_2_models_binary_low_egp_som25_sam005.sh
sbatch cc_2_models_binary_high_egp_som25_sam005.sh
sbatch cc_2_models_tweedie_low_egp_som25_sam005.sh
sbatch cc_2_models_tweedie_high_egp_som25_sam005.sh


# ----------

# EGP SOM 10
## Normal
cp cc_3_marginal_high_egp_som10.sh cc_3_marginal_low_egp_som10.sh
sed -i "s/imbalance_high/imbalance_low/" cc_3_marginal_low_egp_som10.sh 

## Binary
cp cc_3_marginal_high_egp_som10.sh cc_3_marginal_binary_low_egp_som10.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_3_marginal_binary_low_egp_som10.sh 
cp cc_3_marginal_high_egp_som10.sh cc_3_marginal_binary_high_egp_som10.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_3_marginal_binary_high_egp_som10.sh 

## Tweedie
cp cc_3_marginal_high_egp_som10.sh cc_3_marginal_tweedie_low_egp_som10.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_3_marginal_tweedie_low_egp_som10.sh 
cp cc_3_marginal_high_egp_som10.sh cc_3_marginal_tweedie_high_egp_som10.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_3_marginal_tweedie_high_egp_som10.sh 

sbatch cc_3_marginal_low_egp_som10.sh
sbatch cc_3_marginal_high_egp_som10.sh
sbatch cc_3_marginal_binary_low_egp_som10.sh
sbatch cc_3_marginal_binary_high_egp_som10.sh
sbatch cc_3_marginal_tweedie_low_egp_som10.sh
sbatch cc_3_marginal_tweedie_high_egp_som10.sh

# rm -rf ../models/imbalance_low/egp_som10*
# rm -rf ../models/binary_imbalance_low/egp_som10*
# rm -rf ../models/binary_imbalance_high/egp_som10*
# rm -rf ../models/tweedie_imbalance_low/egp_som10*
# rm -rf ../models/tweedie_imbalance_high/egp_som10*


# EGP SOM 50
## Normal
cp cc_3_marginal_high_egp_som50.sh cc_3_marginal_low_egp_som50.sh
sed -i "s/imbalance_high/imbalance_low/" cc_3_marginal_low_egp_som50.sh 

## Binary
cp cc_3_marginal_high_egp_som50.sh cc_3_marginal_binary_low_egp_som50.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_3_marginal_binary_low_egp_som50.sh 
cp cc_3_marginal_high_egp_som50.sh cc_3_marginal_binary_high_egp_som50.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_3_marginal_binary_high_egp_som50.sh 

## Tweedie
cp cc_3_marginal_high_egp_som50.sh cc_3_marginal_tweedie_low_egp_som50.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_3_marginal_tweedie_low_egp_som50.sh 
cp cc_3_marginal_high_egp_som50.sh cc_3_marginal_tweedie_high_egp_som50.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_3_marginal_tweedie_high_egp_som50.sh 

sbatch cc_3_marginal_low_egp_som50.sh
sbatch cc_3_marginal_high_egp_som50.sh
sbatch cc_3_marginal_binary_low_egp_som50.sh
sbatch cc_3_marginal_binary_high_egp_som50.sh
sbatch cc_3_marginal_tweedie_low_egp_som50.sh
sbatch cc_3_marginal_tweedie_high_egp_som50.sh

# rm -rf ../models/imbalance_low/egp_som50*
# rm -rf ../models/binary_imbalance_low/egp_som50*
# rm -rf ../models/binary_imbalance_high/egp_som50*
# rm -rf ../models/tweedie_imbalance_low/egp_som50*
# rm -rf ../models/tweedie_imbalance_high/egp_som50*


# EGP 25 unequal all landscapes
## Normal
cp cc_3_marginal_high_egp_som25_unequal.sh cc_3_marginal_low_egp_som25_unequal.sh
sed -i "s/imbalance_high/imbalance_low/" cc_3_marginal_low_egp_som25_unequal.sh 

## Binary
cp cc_3_marginal_high_egp_som25_unequal.sh cc_3_marginal_binary_low_egp_som25_unequal.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_3_marginal_binary_low_egp_som25_unequal.sh 
cp cc_3_marginal_high_egp_som25_unequal.sh cc_3_marginal_binary_high_egp_som25_unequal.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_3_marginal_binary_high_egp_som25_unequal.sh 

## Tweedie
cp cc_3_marginal_high_egp_som25_unequal.sh cc_3_marginal_tweedie_low_egp_som25_unequal.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_3_marginal_tweedie_low_egp_som25_unequal.sh 
cp cc_3_marginal_high_egp_som25_unequal.sh cc_3_marginal_tweedie_high_egp_som25_unequal.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_3_marginal_tweedie_high_egp_som25_unequal.sh 

sbatch cc_3_marginal_low_egp_som25_unequal.sh
sbatch cc_3_marginal_high_egp_som25_unequal.sh
sbatch cc_3_marginal_binary_low_egp_som25_unequal.sh
sbatch cc_3_marginal_binary_high_egp_som25_unequal.sh
sbatch cc_3_marginal_tweedie_low_egp_som25_unequal.sh
sbatch cc_3_marginal_tweedie_high_egp_som25_unequal.sh

# rm -rf ../models/imbalance_low/egp_som25_unequal*
# rm -rf ../models/binary_imbalance_low/egp_som25_unequal*
# rm -rf ../models/binary_imbalance_high/egp_som25_unequal*
# rm -rf ../models/tweedie_imbalance_low/egp_som25_unequal*
# rm -rf ../models/tweedie_imbalance_high/egp_som25_unequal*


# EGP large sample all landscapes
## Normal
cp cc_3_marginal_high_egp_som25_sam02.sh cc_3_marginal_low_egp_som25_sam02.sh
sed -i "s/imbalance_high/imbalance_low/" cc_3_marginal_low_egp_som25_sam02.sh 

## Binary
cp cc_3_marginal_high_egp_som25_sam02.sh cc_3_marginal_binary_low_egp_som25_sam02.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_3_marginal_binary_low_egp_som25_sam02.sh 
cp cc_3_marginal_high_egp_som25_sam02.sh cc_3_marginal_binary_high_egp_som25_sam02.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_3_marginal_binary_high_egp_som25_sam02.sh 

## Tweedie
cp cc_3_marginal_high_egp_som25_sam02.sh cc_3_marginal_tweedie_low_egp_som25_sam02.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_3_marginal_tweedie_low_egp_som25_sam02.sh 
cp cc_3_marginal_high_egp_som25_sam02.sh cc_3_marginal_tweedie_high_egp_som25_sam02.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_3_marginal_tweedie_high_egp_som25_sam02.sh 

sbatch cc_3_marginal_low_egp_som25_sam02.sh
sbatch cc_3_marginal_high_egp_som25_sam02.sh
sbatch cc_3_marginal_binary_low_egp_som25_sam02.sh
sbatch cc_3_marginal_binary_high_egp_som25_sam02.sh
sbatch cc_3_marginal_tweedie_low_egp_som25_sam02.sh
sbatch cc_3_marginal_tweedie_high_egp_som25_sam02.sh

# rm -rf ../models/imbalance_low/egp_som25_sam02*
# rm -rf ../models/binary_imbalance_low/egp_som25_sam02*
# rm -rf ../models/binary_imbalance_high/egp_som25_sam02*
# rm -rf ../models/tweedie_imbalance_low/egp_som25_sam02*
# rm -rf ../models/tweedie_imbalance_high/egp_som25_sam02*


# EGP large sample all landscapes
## Normal
cp cc_3_marginal_high_egp_som25_sam005.sh cc_3_marginal_low_egp_som25_sam005.sh
sed -i "s/imbalance_high/imbalance_low/" cc_3_marginal_low_egp_som25_sam005.sh 

## Binary
cp cc_3_marginal_high_egp_som25_sam005.sh cc_3_marginal_binary_low_egp_som25_sam005.sh
sed -i "s/imbalance_high/binary_imbalance_low/" cc_3_marginal_binary_low_egp_som25_sam005.sh 
cp cc_3_marginal_high_egp_som25_sam005.sh cc_3_marginal_binary_high_egp_som25_sam005.sh
sed -i "s/imbalance_high/binary_imbalance_high/" cc_3_marginal_binary_high_egp_som25_sam005.sh 

## Tweedie
cp cc_3_marginal_high_egp_som25_sam005.sh cc_3_marginal_tweedie_low_egp_som25_sam005.sh
sed -i "s/imbalance_high/tweedie_imbalance_low/" cc_3_marginal_tweedie_low_egp_som25_sam005.sh 
cp cc_3_marginal_high_egp_som25_sam005.sh cc_3_marginal_tweedie_high_egp_som25_sam005.sh
sed -i "s/imbalance_high/tweedie_imbalance_high/" cc_3_marginal_tweedie_high_egp_som25_sam005.sh 

sbatch cc_3_marginal_low_egp_som25_sam005.sh
sbatch cc_3_marginal_high_egp_som25_sam005.sh
sbatch cc_3_marginal_binary_low_egp_som25_sam005.sh
sbatch cc_3_marginal_binary_high_egp_som25_sam005.sh
sbatch cc_3_marginal_tweedie_low_egp_som25_sam005.sh
sbatch cc_3_marginal_tweedie_high_egp_som25_sam005.sh

# rm -rf ../models/imbalance_low/egp_som25_sam005*
# rm -rf ../models/binary_imbalance_low/egp_som25_sam005*
# rm -rf ../models/binary_imbalance_high/egp_som25_sam005*
# rm -rf ../models/tweedie_imbalance_low/egp_som25_sam005*
# rm -rf ../models/tweedie_imbalance_high/egp_som25_sam005*



# Default model
sbatch cc_2_models_low_egp_som25.sh
sbatch cc_2_models_high_egp_som25.sh


sbatch cc_3_marginal_low_egp_som25.sh
sbatch cc_3_marginal_high_egp_som25.sh
sbatch cc_3_marginal_binary_low_egp_som25.sh
sbatch cc_3_marginal_binary_high_egp_som25.sh
sbatch cc_3_marginal_tweedie_low_egp_som25.sh
sbatch cc_3_marginal_tweedie_high_egp_som25.sh



# rm -rf ../models/imbalance_low/egp_som25.*
# rm -rf ../models/binary_imbalance_low/egp_som25.*
# rm -rf ../models/binary_imbalance_high/egp_som25.*
# rm -rf ../models/tweedie_imbalance_low/egp_som25.*
# rm -rf ../models/tweedie_imbalance_high/egp_som25.*
