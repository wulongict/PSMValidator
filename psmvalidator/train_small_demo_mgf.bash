# run the script to train a model
#

set -v

# delete some results file if exist!
basename=human_hcd_selected_new_small
rm ../train/${basename}__2.000000_frag.model
rm ../train/${basename}.mgffeatures_.feat
rm ../train/${basename}.mgffeatures_trN_90000__PNcombined.feat
rm ../train/decoy_${basename}.pep.xmlfeatures_.feat
rm ../train/decoy_${basename}.pep.xml

current_path=`pwd`
output_file_conf=psmvalidator_train_small_human_demo_model.conf
basename=human_hcd_selected_new_small

ntree=900
mtry=8

cat <<-EOF > $output_file_conf
task = all

# ----------------------------------------------------
# DATA: NIST human data(small demo)
#----------------------------------------------------
inputfile=${current_path}/../train/${basename}.mgf

# -----------------------------------------------
# DATA: Comet search parameter comet2016
#       NIST human spectral library can be regarded as LOW mass accuracy.
#-------------------------------------------------
cometparam=${current_path}/../train/comet16low.param
cometbinary=${current_path}/../bin/comet

#--------------------------------------------------
# PeakPair Model:
#------------------------------------------------------
ghostpeak=true
use_output=true

# -------------------------------------------
# DATA: Huamn Uniprot sequence database
#-----------------------------------------------
targetdb=${current_path}/../train/uniprot-human-2020-12.fasta
decoydb=${current_path}/../train/uniprot-human-2020-12TD_only_decoy.fasta

#----------------------------
# Training with ranger.
#----------------------------
trainValidator=true
rangerbinary=${current_path}/../bin/ranger
ntree=${ntree}
mtry=${mtry}

featurelistfile=${current_path}/../train/features.txt

writeValidatorModel=${current_path}/../train/${basename}_mtry_${mtry}_ntree_${ntree}_features.forest

EOF

PSMValidator=${current_path}/../bin/psmvalidator

$PSMValidator  --config $output_file_conf

