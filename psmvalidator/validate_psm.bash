set -x
current_path=`pwd`

model=../model/nist_human_hcd_selected_features.txtmtry_8_ntree_900_trN_20000.forest # todo
featurelist=../train/features.txt
truthfile=../predict/TUM_first_pool_12.pep.list # todo
pepxml=../predict/interact-01625b_GD2-TUM_first_pool_12_01_01-3xHCD-1h-R1.ipro.pep.xml # todo

output_file_conf=psmvalidator_validate_psm.conf
ntree=900
mtry=8

cat <<-EOF > $output_file_conf
task = validator
inputfile=${pepxml}
trainValidator=true
rangerbinary=${current_path}/../bin/ranger
ntree=${ntree}
mtry=${mtry}
featurelistfile=${current_path}/${featurelist}
validatorModel=${current_path}/${model}
truthfile=${truthfile}
EOF

PSMValidator=${current_path}/../bin/psmvalidator
$PSMValidator  --config $output_file_conf