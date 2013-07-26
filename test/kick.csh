#!/bin/tcsh

set path_=/castor/cern.ch/user/y/ymtzeng/Summer11_MC_JES_V9_bpkit1009/
set list_=`rfdir ${path_} | awk '{print $9}'`

foreach lt_($list_)

set file_=`rfdir ${path_}${lt_} | awk '{print $9}'`

foreach fl_($file_)

rfcp ${path_}${lt_}/${fl_} /tmp/ymtzeng/${fl_} &
#sleep 1
ps aux | grep ymtzeng | grep ".root" | awk '{print "kill -9 "$2}' |csh
echo "${path_}${lt_}/${fl_}"

if ( -e /tmp/ymtzeng/${fl_} ) then
rm /tmp/ymtzeng/${fl_}
endif

end
end
