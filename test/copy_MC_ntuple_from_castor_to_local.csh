#!/bin/tcsh

if ( $1 == "" ) then
echo "PLEASE INPUT FOLDER name like : "
echo "copy_MC_ntuple_from_castor_to_local.csh DYJets 3"
exit
endif

#set folder_name=HitFit/20120228_test/Data
set folder_name=HitFit/20120228_test/MC

set PATH_=/castor/cern.ch/user/t/twang/$folder_name/
#set folder_name=vorobiev
#set PATH_=/castor/cern.ch/user/v/$folder_name/
if ( $1 != "all" ) then
set datasets=`rfdir $PATH_ | grep $1 | awk '{print $9}'`
endif
if ( $1 == "all" ) then
set datasets=`rfdir $PATH_ | awk '{print $9}'`
endif

set date_=`date | sed 's/ /_/g' | sed 's/:/_/g'`
touch error_log_$date_
touch problematic_files_$date_
set isOkay=1

foreach ds($datasets)

   if (!(-e $folder_name/$ds)) then
      mkdir -p $folder_name/$ds
   endif

   set file_list=`rfdir $PATH_/$ds | awk '{print $9}'| grep results`
   foreach flt($file_list)
	set check=`echo "0"`
	@ retries=0
	while ($check == "0")	
		set copy_size=`ls -l $folder_name/$ds/$flt | awk '{print $5}'`
		set original_size=`rfdir $PATH_/$ds/$flt | awk '{print $5}'`

		if ($original_size < 1000 ) then
			rfdir $PATH_/$ds/$flt >> problematic_files_$date_
		endif		

 		echo "original_size = $original_size ; copy_size =  $copy_size"
		if ( "$original_size" != "$copy_size" ) then
      			rfcp  $PATH_/$ds/$flt $folder_name/$ds/$flt
		else if ( "$original_size" == "$copy_size" ) then
			set check=`echo "1"`
		endif
		@ retries++
		if ( $retries == "10" ) then
            echo "Not ready for copying "$PATH_/$ds/$flt >> error_log_$date_
            set isOkay=0
			break
		endif
	end	
   end

end

if ( "$isOkay" == "1" ) then
   rm error_log_$date_
endif

set problem=`ll problematic_files_$date_ | awk '{print $5}'`
if ($problem == 0) then
rm problematic_files_$date_
endif


