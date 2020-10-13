#!/bin/bash

# Script to process merra-2 met fields for 200404
# MERRA2_100.*.19800101.nc4
# MERRA2_200.*.19920101.nc4
# MERRA2_300.*.20010101.nc4
# MERRA2_400.*.20110101.nc4
# leap years: 2000 2004 2008 2012 2016


leap_yrs=(1980 1984 1988 1992 1996 2000 2004 2008 2012 2016 2020)
daysinmonth=(31 28 31 30 31 30 31 31 30 31 30 31)

rawpath=/data10/jmeng/MERRA2_Download/MERRA2_Glob/merra2_raw     #directory where the raw MERRA-2 data is stored
tempdir=/data10/jmeng/MERRA2_Download/MERRA2_Glob/merra2_output/           #directory of outputs               
processdir=/data10/jmeng/MERRA2_Download/MERRA2_Glob/merra2_code            #directory containing the processing code
ctmdir=/data10/ctm
ctmdir_scratch=/stetson-home/jmeng/scratch_geoschem/ctm_scratch

year=1983

GRIDNAME=("GEOS_0.5x0.625_AS.d" "GEOS_0.5x0.625_NA.d" "GEOS_0.5x0.625_EU.d" "GEOS_2x2.5.d" "GEOS_4x5.d")
GRIDNAME_scratch=("GEOS_0.5x0.625.d")


sed -i "95s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "96s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "102s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "103s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "109s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "110s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "116s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "117s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "123s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "124s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "129s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "130s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "135s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
sed -i "136s\.*\/$tempdir\g" $processdir/bin/Merra2_Driver.input
    

if [ $year -lt 1992 ]; then

		sed -i "40s\.*\MERRA2_100.tavg1_2d_flx_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "45s\.*\MERRA2_100.tavg1_2d_lnd_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "50s\.*\MERRA2_100.tavg1_2d_rad_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "55s\.*\MERRA2_100.tavg1_2d_slv_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "60s\.*\MERRA2_100.tavg3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "65s\.*\MERRA2_100.tavg3_3d_cld_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "71s\.*\MERRA2_100.tavg3_3d_mst_Ne.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "76s\.*\MERRA2_100.tavg3_3d_mst_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "81s\.*\MERRA2_100.tavg3_3d_rad_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "86s\.*\MERRA2_100.inst3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input

elif [ $year -ge 1992 ] && [ $year -lt 2001 ]; then

		sed -i "40s\.*\MERRA2_200.tavg1_2d_flx_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "45s\.*\MERRA2_200.tavg1_2d_lnd_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "50s\.*\MERRA2_200.tavg1_2d_rad_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "55s\.*\MERRA2_200.tavg1_2d_slv_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "60s\.*\MERRA2_200.tavg3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "65s\.*\MERRA2_200.tavg3_3d_cld_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "71s\.*\MERRA2_200.tavg3_3d_mst_Ne.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "76s\.*\MERRA2_200.tavg3_3d_mst_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "81s\.*\MERRA2_200.tavg3_3d_rad_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "86s\.*\MERRA2_200.inst3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input

elif [ $year -ge 2001 ] && [ $year -lt 2011 ]; then

		sed -i "40s\.*\MERRA2_300.tavg1_2d_flx_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "45s\.*\MERRA2_300.tavg1_2d_lnd_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "50s\.*\MERRA2_300.tavg1_2d_rad_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "55s\.*\MERRA2_300.tavg1_2d_slv_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "60s\.*\MERRA2_300.tavg3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "65s\.*\MERRA2_300.tavg3_3d_cld_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "71s\.*\MERRA2_300.tavg3_3d_mst_Ne.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "76s\.*\MERRA2_300.tavg3_3d_mst_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "81s\.*\MERRA2_300.tavg3_3d_rad_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "86s\.*\MERRA2_300.inst3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input

elif [ $year -ge 2011 ]; then

		sed -i "40s\.*\MERRA2_400.tavg1_2d_flx_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "45s\.*\MERRA2_400.tavg1_2d_lnd_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "50s\.*\MERRA2_400.tavg1_2d_rad_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "55s\.*\MERRA2_400.tavg1_2d_slv_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "60s\.*\MERRA2_400.tavg3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "65s\.*\MERRA2_400.tavg3_3d_cld_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "71s\.*\MERRA2_400.tavg3_3d_mst_Ne.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "76s\.*\MERRA2_400.tavg3_3d_mst_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "81s\.*\MERRA2_400.tavg3_3d_rad_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "86s\.*\MERRA2_400.inst3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input

fi


for month in $(seq  11  11); do
	   
	 # Change num of days in Feb to 29 if is a leap year
     if [ $month -eq 2 ]; then
        if [[ " ${leap_yrs[*]} " == *" $year "* ]]; then
           daysinmonth[1]=29
         else daysinmonth[1]=28
        fi
     fi
	 
	 if [ ! ${#month} -eq 2 ]; then # ${#month} returns the length of a string
		 strmonth="0$month"
                 else strmonth=$month
	 fi	 
      
      ## make directory in ctm directory 
      for grid in "${GRIDNAME[@]}" ; do
        dir=$ctmdir/$grid/MERRA2/$year/$strmonth/
        if [ ! -d "$dir" ]; then # check if directory exist!
          mkdir -p $dir
        fi
      done

      ## make directory in ctm_scratch directory
      for grid in "${GRIDNAME_scratch[@]}" ; do
           dir=$ctmdir_scratch/$grid/MERRA2/$year/$strmonth/
           if [ ! -d "$dir" ]; then # check if directory exist!
                  mkdir -p $dir
           fi
      done      

     
     # Change the raw met fields directory in input files		  
     sed -i "90s\.*\/$rawpath/$year$strmonth/\g" $processdir/bin/Merra2_Driver.input 
     sed -i "31s\.*\/$rawpath/$year$strmonth/\g" $processdir/perl/doMerra2.input  	 
     
     rawdir=$rawpath/$year$strmonth
     	 
     # Loop over all days in a month to download and process data
	 
	   for day in $(seq -w  22 22); do #${daysinmonth[10#$month-1]}); do     #${daysinmonth[10#$month-1]}
		   
  	  	 if [ ! ${#day} -eq 2 ]; then # ${#month} returns the length of a string
  	  		 strday="0$day"
	     else strday=$day
  	  	 fi
 		   
           date=$year$strmonth$strday
		   cd $processdir/perl
		   
	     # Regrid data of the day
		   printf " Regridding $date \n"	
           	./doMerra2 $date
       done
      
	   # Compress nc files to level 5
	   printf " Compressing $year$strmonth \n"
	   
	   (cd $tempdir; chmod 755 *;

	   for f in *$year$strmonth*.nc4; do
	       
	       mv $f $f.big
	       line=$(ncdump -h $f.big | grep "lat = ")  #retrieve number of lats
	       subline=${line##*= }
	       nlat=${subline% ;*}
	       
	       line=$(ncdump -h $f.big | grep "lon = ")   #retrieve number of lons
	       subline=${line##*= }
	       nlon=${subline% ;*}
	       
	       line=$(ncdump -h $f.big | grep -n "dimensions")  #retrieve number of dimensions
	       line2=$(ncdump -h $f.big | grep -n "variables")
	       ndim=$((${line2%:v*}-${line%:d*}-1))
	       
	       if [ $ndim == 3 ]; then   #2-D fields
	                  
		   nccopy -w -k3 -d5 -c time/1,lat/$nlat,lon/$nlon $f.big $f
	       fi
	       
	       if [ $ndim == 4 ]; then   #3-D fields
	           
		   nccopy -w -k3 -d5 -c time/1,lev/1,lat/$nlat,lon/$nlon $f.big $f
	       
	       fi
	       
	   done)
           
	   
	  
	  (cd $tempdir; chmod 755 *)
		   #mv $tempdir/*$year$strmonth*05x0625.CH*.nc4 $ctmdir/GEOS_0.5x0.625_CH.d/MERRA2/$year/$strmonth/  #No CH created. soft link from AS to CH
		    mv $tempdir/*$year$strmonth*05x0625.NA*.nc4 $ctmdir/GEOS_0.5x0.625_NA.d/MERRA2/$year/$strmonth/
		    mv $tempdir/*$year$strmonth*05x0625.EU*.nc4 $ctmdir/GEOS_0.5x0.625_EU.d/MERRA2/$year/$strmonth/
		    mv $tempdir/*$year$strmonth*05x0625.AS*.nc4 $ctmdir/GEOS_0.5x0.625_AS.d/MERRA2/$year/$strmonth/
		    
                   ## save native resolution global files to scratch space on stetson 
                   #mv $tempdir/*$year$strmonth*05x0625.nc4 $ctmdir/GEOS_0.5x0.625.d/MERRA2/$year/$strmonth/
		    mv $tempdir/*$year$strmonth*05x0625.nc4 $ctmdir_scratch/GEOS_0.5x0.625.d/MERRA2/$year/$strmonth/
		    
	            mv $tempdir/*$year$strmonth*2x25*.nc4 $ctmdir/GEOS_2x2.5.d/MERRA2/$year/$strmonth/
		    mv $tempdir/*$year$strmonth*4x5*.nc4 $ctmdir/GEOS_4x5.d/MERRA2/$year/$strmonth/
		   
                    rm $tempdir/*$year$strmonth*nc4.big
		   
		   
	  
done

#remove raw files to save space
#rm -rf /data10/jmeng/MERRA2_Download/MERRA2_Glob/merra2_raw/* 

printf " End of script!\n"


