#!/bin/bash

GEOS_FP_DOWNLOAD_DIR=/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/MetFieldProcessing/GEOS_FP-raw

ProDir=/data10/jmeng/GEOS_FP_Download/regrid_codes/
year='2020'

ctmdir_scratch=/stetson-home/jmeng/scratch_geoschem/ctm_scratch

daysinmonth=(31 29 31 30 31 30 31 31 30 31 30 31)  #make sure you made modifications in leap years


for month in $(seq -w 8 8);

do
  if [ ! ${#month} -eq 2 ]; then # ${#month} returns the length of a string
		 strmonth="0$month"
                 else strmonth=$month
  fi	
# month=0$month
# Change the raw met fields directory in input files		  
 sed -i "96s\.*\/stetson-home/jmeng/scratch_geoschem/GEOS-FP_raw/$year$strmonth/\g" $ProDir/bin/GeosFpDriver.input
 sed -i "39s\.*\/stetson-home/jmeng/scratch_geoschem/GEOS-FP_raw/$year$strmonth/\g" $ProDir/perl/doGeosFp.input
 
 sed -i "185s\.*\/data10/jmeng/GEOS_FP_Download/regrid_codes/temp/\g" $ProDir/bin/GeosFpDriver.input
 sed -i "186s\.*\/data10/jmeng/GEOS_FP_Download/regrid_codes/temp/\g" $ProDir/bin/GeosFpDriver.input


GRIDNAME=("GEOS_0.25x0.3125_CH.d" "GEOS_0.25x0.3125_NA.d" "GEOS_0.25x0.3125_EU.d"  "GEOS_0.25x0.3125_AS.d"\
   "GEOS_0.5x0.625_CH.d" "GEOS_0.5x0.625_NA.d" "GEOS_0.5x0.625_EU.d" "GEOS_0.5x0.625_AS.d"\
   "GEOS_2x2.5.d" "GEOS_4x5.d")

GRIDNAME_scratch=("GEOS_0.25x0.3125.d")

# Create processed data directory
for grid in "${GRIDNAME[@]}" ; do
 dir=/data10/ctm/$grid/GEOS_FP/$year/$strmonth/
 if [ ! -d "$dir" ]; then # check if directory exist!
    mkdir -p $dir
 fi
done


# Create directory in ctm_scratch directory for global native resolution met field
for grid in "${GRIDNAME_scratch[@]}" ; do
	dir=$ctmdir_scratch/$grid/GEOS_FP/$year/$strmonth/
	if [ ! -d "$dir" ]; then # check if directory exist!
		mkdir -p $dir 
	fi
done

# Process data
echo "---------started-------"
cd $ProDir
export PATH=$ProDir/bin:$PATH
cd perl



 for day in $(seq -w 1 ${daysinmonth[10#$month-1]}); #${daysinmonth[10#$month-1]}

 #for day in '09' '10' '28';
 do
    if [ ! ${#day} -eq 2 ]; then # ${#month} returns the length of a string
		 strday="0$day"
                 else strday=$day
   fi	
   date=$year$strmonth$strday

  ./doGeosFpMulti $date

 done
 

printf "\n nc3 is finished!\n"

# Convert nc3 to nc4

directory=$ProDir/temp
ctmdir=/data10/ctm       
nc4dir=$ProDir/temp/nc4    

if [ ! -d "$nc4dir" ]; then # check if directory exist!
   mkdir -p $nc4dir
fi

(cd $directory; chmod 755 *;

     for f in *$year$strmonth*.nc*; do
	       
	       
	       line=$(ncdump -h $f | grep "lat = ")  #retrieve number of lats
	       subline=${line##*= }
	       nlat=${subline% ;*}
	       
	       line=$(ncdump -h $f | grep "lon = ")   #retrieve number of lons
	       subline=${line##*= }
	       nlon=${subline% ;*}
	       
	       line=$(ncdump -h $f | grep -n "dimensions")  #retrieve number of dimensions
	       line2=$(ncdump -h $f | grep -n "variables")
	       ndim=$((${line2%:v*}-${line%:d*}-1))
	       
	       if [ $ndim == 3 ]; then   #2-D fields
	                  
		   nccopy -w -k3 -d5 -c time/1,lat/$nlat,lon/$nlon $f $nc4dir/$f
	       fi
	       
	       if [ $ndim == 4 ]; then   #3-D fields
	           
		   nccopy -w -k3 -d5 -c time/1,lev/1,lat/$nlat,lon/$nlon $f $nc4dir/$f
	       
	       fi
	       
     done)

if [ ! $(ls $directory/*$year$strmonth*.nc |wc -l) -eq $(ls $nc4dir/*$year$strmonth*.nc |wc -l) ]; then
   printf "\n $directory/*$year$strmonth*.nc3 not same number with $nc4dir/*$year$strmonth*.nc4!\n"
else
   rm $directory/*$year$strmonth*
fi
printf "\nconvert $grid to nc4 is finished!\n"

# Move proessed met fields from my directory to /acenet/shared/ctm
(cd $nc4dir; chmod 755 *)
  
  mv $nc4dir/*$year$strmonth*025x03125.CH* $ctmdir/GEOS_0.25x0.3125_CH.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*025x03125.NA* $ctmdir/GEOS_0.25x0.3125_NA.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*025x03125.EU* $ctmdir/GEOS_0.25x0.3125_EU.d/GEOS_FP/$year/$strmonth/
#  mv $nc4dir/*$year$strmonth*025x03125.SE* $ctmdir/GEOS_0.25x0.3125_SE.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*025x03125.AS* $ctmdir/GEOS_0.25x0.3125_AS.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*05x0625.CH* $ctmdir/GEOS_0.5x0.625_CH.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*05x0625.NA* $ctmdir/GEOS_0.5x0.625_NA.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*05x0625.EU* $ctmdir/GEOS_0.5x0.625_EU.d/GEOS_FP/$year/$strmonth/
#  mv $nc4dir/*$year$strmonth*05x0625.SE* $ctmdir/GEOS_0.5x0.625_SE.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*05x0625.AS* $ctmdir/GEOS_0.5x0.625_AS.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*2x25* $ctmdir/GEOS_2x2.5.d/GEOS_FP/$year/$strmonth/
  mv $nc4dir/*$year$strmonth*4x5* $ctmdir/GEOS_4x5.d/GEOS_FP/$year/$strmonth/
  # move global native met field into scratch space 
  mv $nc4dir/*$year$strmonth*025x03125.nc $ctmdir_scratch/GEOS_0.25x0.3125.d/GEOS_FP/$year/$strmonth/

done

printf "\nEnd of scriptf!\n"
