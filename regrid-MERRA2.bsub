#!/bin/bash
#BSUB -n 12
#BSUB -R "rusage[mem=150GB] span[hosts=1] order[-slots]"
#BSUB -q rvmartin
#BSUB -W 168:00
#BSUB -w "DownloadMERRA2"
#BSUB -a 'docker(registry.gsc.wustl.edu/sleong/base-engineering)'
#BSUB -g /liam.bindle/downloads
#BSUB -J "RegridMERRA2"
#BSUB -Ne
#BSUB -u liam.bindle@wustl.edu
#BSUB -o /storage1/fs1/rvmartin/Active/GEOS-Chem-shared/MetFieldProcessing/logs/regrid-MERRA2_${YEAR}${MONTH}-%J.txt

# Script to process merra-2 met fields for 200404
# MERRA2_100.*.19800101.nc4
# MERRA2_200.*.19920101.nc4
# MERRA2_300.*.20010101.nc4
# MERRA2_400.*.20110101.nc4
# leap years: 2000 2004 2008 2012 2016
. /etc/bashrc

set -x
set -u
: "$YEAR" "$MONTH"

ulimit -s unlimited


rawpath=/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/MetFieldProcessing/MERRA2-raw/$YEAR/$MONTH/
tempdir=/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/MetFieldProcessing/scratch/MERRA2-raw/$YEAR/$MONTH/
processdir=/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/MetFieldProcessing/MERRA2

mkdir -p $tempdir

GRIDNAME=("GEOS_0.5x0.625_AS.d" "GEOS_0.5x0.625_NA.d" "GEOS_0.5x0.625_EU.d" "GEOS_2x2.5.d" "GEOS_4x5.d")
GRIDNAME_scratch=("GEOS_0.5x0.625.d")


sed -i "95s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "96s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "102s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "103s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "109s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "110s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "116s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "117s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "123s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "124s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "129s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "130s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "135s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
sed -i "136s#.*#$tempdir#" $processdir/bin/Merra2_Driver.input
    

if [ $YEAR -lt 1992 ]; then

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

elif [ $YEAR -ge 1992 ] && [ $YEAR -lt 2001 ]; then

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

elif [ $YEAR -ge 2001 ] && [ $YEAR -lt 2011 ]; then

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

elif [ "$YEAR" = "2020" ] && [ "$MONTH" = "09" ]; then

		sed -i "40s\.*\MERRA2_401.tavg1_2d_flx_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "45s\.*\MERRA2_401.tavg1_2d_lnd_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "50s\.*\MERRA2_401.tavg1_2d_rad_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "55s\.*\MERRA2_401.tavg1_2d_slv_Nx.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "60s\.*\MERRA2_401.tavg3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "65s\.*\MERRA2_401.tavg3_3d_cld_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "71s\.*\MERRA2_401.tavg3_3d_mst_Ne.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "76s\.*\MERRA2_401.tavg3_3d_mst_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "81s\.*\MERRA2_401.tavg3_3d_rad_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input
		sed -i "86s\.*\MERRA2_401.inst3_3d_asm_Nv.YYYYMMDD.nc4\g" $processdir/bin/Merra2_Driver.input

elif [ $YEAR -ge 2011 ]; then

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

# Change the raw met fields directory in input files		  
sed -i "90s#.*#$rawpath#" $processdir/bin/Merra2_Driver.input 
sed -i "31s#.*#$rawpath#" $processdir/perl/doMerra2.input  	 

rawdir=$rawpath
   
# Loop over all days in a MONTH to download and process data

DAYS_IN_MONTH=$(cal $MONTH $YEAR | awk 'NF {DAYS = $NF}; END {print DAYS}')

for day in $(seq -w  1 ${DAYS_IN_MONTH}); do
   date=$YEAR$MONTH$(printf '%02i' "${day#0}")
   cd $processdir/perl
   ./doMerra2 $date
done

# Compress nc files to level 5
printf " Compressing $YEAR$MONTH \n"

cd $tempdir

for f in *$YEAR$MONTH*.nc4; do
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
   elif [ $ndim == 4 ]; then   #3-D fields
      nccopy -w -k3 -d5 -c time/1,lev/1,lat/$nlat,lon/$nlon $f.big $f
   fi
   
done

printf " End of script!\n"


