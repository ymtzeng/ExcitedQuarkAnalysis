#!/bin/csh


set outputFolder_=/castor/cern.ch/user/y/ymtzeng/ExcitedQuark_8j1bor2b_v6
set location_=`pwd`


if ( !(-e ${location_}/submitSamples) ) then
mkdir ${location_}/submitSamples
endif

set sample_=Summer11_MC_JES_V9_bpkit1009
#set sample_=CMSSW42X_data_DCSONLY_PFLeptonAddLooseLep1012
set path_=/castor/cern.ch/user/y/ymtzeng/${sample_}
set list_=`rfdir ${path_} | awk '{print $9}'`
#set list_=`rfdir ${path_} | awk '{print $9}' | grep SingleMu`

rfmkdir ${outputFolder_}
rfchmod 777 ${outputFolder_}

foreach lt_($list_)

rfmkdir ${outputFolder_}/${lt_}
rfchmod 777 ${outputFolder_}/${lt_}
mkdir ${location_}/submitSamples/${lt_}

set filelist_=`rfdir ${path_}/${lt_} | awk '{print $9}'`

foreach fl_(${filelist_})

set fl__=`echo $fl_ | awk -F"." '{print $1}'`


set check_file_in_castor=`rfdir ${outputFolder_}/${lt_}/${lt_}_${fl__}.root | grep ${fl__} | wc | awk '{print $1}'`
if ( ${check_file_in_castor} == 0 ) then

echo "${lt_}_${fl__}.root"

if ( -e ${location_}/submitSamples/${lt_}/tprimeanalysis_cfg_${lt_}_${fl__}.py ) then
rm ${location_}/submitSamples/${lt_}/tprimeanalysis_cfg_${lt_}_${fl__}.py
endif

sed "s/INPUT_FOLDER/${lt_}/g" tprimeanalysis_cfg.py |sed "s/CASTOR_FOLDER/${sample_}/g" | sed "s/ROOTFILE/${fl_}/g" | sed "s/OUTPUT_FILE/${lt_}_${fl_}/g" >& ${location_}/submitSamples/${lt_}/tprimeanalysis_cfg_${lt_}_${fl__}.py

if ( -e ${location_}/submitSamples/${lt_}/submit_${lt_}_${fl__}.csh ) then 
rm ${location_}/submitSamples/${lt_}/submit_${lt_}_${fl__}.csh
endif

echo "#\!/bin/csh \
cd ${location_}/submitSamples/${lt_} \
cmsenv \
cmsRun ${location_}/submitSamples/${lt_}/tprimeanalysis_cfg_${lt_}_${fl__}.py \
rfcp /tmp/ymtzeng/${lt_}_${fl__}.root ${outputFolder_}/${lt_}/${lt_}_${fl__}.root" >& ${location_}/submitSamples/${lt_}/submit_${lt_}_${fl__}.csh

chmod +x ${location_}/submitSamples/${lt_}/submit_${lt_}_${fl__}.csh

bsub -q 8nh -J ${lt_}_${fl__} < ${location_}/submitSamples/${lt_}/submit_${lt_}_${fl__}.csh

endif

end
end
