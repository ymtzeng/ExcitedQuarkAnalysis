#!/bin/csh

set workspace=`pwd`
set sourceFiles=/nasdata2/ymtzeng/ExcitedQuark_v4
set IsAll_=0

if ( "$1" == "" ) then
echo "Please use ./link.csh SingleMu"
echo "or ./link.csh all"
exit
else if ( "$1" == "all" ) then
set IsAll_=1
endif

set file_folder=`ls ${sourceFiles} `

if ( ${IsAll_} == "0" ) then
set file_folder=`ls ${sourceFiles} | grep $1`
endif

mkdir ${workspace}/REDUCE_DATA
foreach ff_(${file_folder})
if ( -e ${workspace}/REDUCE_DATA/${ff_}_results_.root ) then
	rm ${workspace}/REDUCE_DATA/${ff_}_results_.root
endif

set Nfiles_=`ls ${sourceFiles}/${ff_} | wc | awk '{print $1}' `

if ( ${Nfiles_} > 999 ) then
	
hadd ${workspace}/REDUCE_DATA/${ff_}_results_1.root ${sourceFiles}/${ff_}/*results_[1-9]_*root
hadd ${workspace}/REDUCE_DATA/${ff_}_results_2.root ${sourceFiles}/${ff_}/*results_[1-9]?_*root
hadd ${workspace}/REDUCE_DATA/${ff_}_results_3.root ${sourceFiles}/${ff_}/*results_[1-9]??_*root
hadd ${workspace}/REDUCE_DATA/${ff_}_results_4.root ${sourceFiles}/${ff_}/*results_1???_*root
hadd ${workspace}/REDUCE_DATA/${ff_}_results_.root ${workspace}/REDUCE_DATA/${ff_}_results_?.root  
rm ${workspace}/REDUCE_DATA/${ff_}_results_?.root

else 
hadd ${workspace}/REDUCE_DATA/${ff_}_results_.root ${sourceFiles}/${ff_}/*root
endif

end

#hadd ${workspace}/REDUCE_DATA/SingleMu_results.root ${workspace}/REDUCE_DATA/SingleMu*root
