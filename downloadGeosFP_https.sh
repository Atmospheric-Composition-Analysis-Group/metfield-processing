#!/bin/bash

year='2020'
month='08'
days='31'


# Make raw and processed data directoreis if not already exist
# rawdir=/data10/jmeng/GEOS_FP_Download/GEOS-FP_raw/$year$month/
# Now download raw GEOS-FP raw files to scratch space
rawdir=/stetson-home/jmeng/scratch_geoschem/GEOS-FP_raw/$year$month/
if [ ! -d "$rawdir" ]; then # check if directory exist!
    mkdir -p $rawdir
fi

echo "---------Started -------"

for day in $(seq -w '01' $days) ; 

do
  date=$year$month$day

   echo "Downloading $date"
   
   cd /stetson-home/jmeng/scratch_geoschem/GEOS-FP_raw/$year$month/
    wget  -r  --no-parent   --cut-dirs=9 -A "*tavg1_2d_rad_Nx*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget  -r  --no-parent   --cut-dirs=9 -A "*tavg1_2d_lnd_Nx*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent   --cut-dirs=9 -A  "*tavg1_2d_flx_Nx*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent  --cut-dirs=9 -A  "*tavg1_2d_slv_Nx*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent  --cut-dirs=9 -A  "*inst3_3d_asm_Nv*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent   --cut-dirs=9 -A  "*tavg3_3d_asm_Nv*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent  --cut-dirs=9 -A  "*tavg3_3d_cld_Nv*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent   --cut-dirs=9 -A  "*tavg3_3d_mst_Ne*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent  --cut-dirs=9 -A  "*tavg3_3d_rad_Nv*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/
   wget -r  --no-parent  --cut-dirs=9 -A  "*tavg3_3d_mst_Nv*nc4" https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y$year/M$month/D$day/ 
     
done
mv portal.nccs.nasa.gov/* .
rm -rf portal.nccs.nasa.gov
echo "End of code"
