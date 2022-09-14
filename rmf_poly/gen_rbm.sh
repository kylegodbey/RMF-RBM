#!/bin/bash

ulimit -s unlimited

outfile=rmf_rbm_hybrid

rm $outfile.pyx
rm $outfile.f90

cd HybridPolynomials

hybrid=0

nbasis=$hybrid
nfields=4

echo "#!python
#cython: language_level=3
#cython: extra_compile_args=[\"-Ofast\"]

import numpy as np
import cython

cimport numpy as np

np.import_array()

DTYPE = np.float64

ctypedef np.float64_t DTYPE_t
" > ../$outfile.pyx

for nuc in 16O 40Ca 48Ca 68Ni 90Zr 100Sn 116Sn 132Sn 144Sm 208Pb; do
#for nuc in 48Ca 208Pb; do
#3	3		6	6		6	7		7	10		10 11		11	11		11	14		11	16		13	16		16	22
if [ $nuc = "16O" ]; then
nstates=3
pstates=3
elif [ $nuc = "40Ca" ]; then
nstates=6
pstates=6
elif [ $nuc = "48Ca" ]; then
nstates=7
pstates=6
elif [ $nuc  = "68Ni" ]; then
nstates=10
pstates=7
elif [ $nuc  = "90Zr" ]; then
nstates=11
pstates=10
elif [ $nuc  = "100Sn" ]; then
nstates=11
pstates=11
elif [ $nuc  = "116Sn" ]; then
nstates=14
pstates=11
elif [ $nuc  = "132Sn" ]; then
nstates=16
pstates=11
elif [ $nuc  = "144Sm" ]; then
nstates=16
pstates=13
elif [ $nuc  = "208Pb" ]; then
nstates=22
pstates=16
fi

dir="Results_"$nuc*
cd $dir

basisfile="Nucleus_"$nuc"BasisNumbers.txt"

### POLYNOMIAL EQUATION ###

# org_file="Nucleus_"$nuc"_BasisNumber3_CoeffEquations.txt"
# file="Nucleus_"$nuc"_BasisNumber3_CoeffEquations.tmp"
org_file="Nucleus_"$nuc"_CoeffEquations.txt"
file="Nucleus_"$nuc"_CoeffEquations.tmp"
cp $org_file $file

echo "
@cython.boundscheck(False)
@cython.wraparound(False)
def rmf_poly_"$nuc"_"$hybrid"(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params):

	cdef np.ndarray[DTYPE_t] result = np.ones([x0.shape[0]], dtype=DTYPE)
" > rmf_poly-$nuc.tmp

cp ../../parse_math.sh .
./parse_math.sh $nbasis $nfields $nstates $pstates $file $basisfile

echo "" >> $file
i=0
while read l; do
  echo -e "\tresult[$i]=$l" >> py.tmp
  i=$((i+1))
done < $file

cat rmf_poly-$nuc.tmp py.tmp > rmf_poly_$nuc-$hybrid.tmp

echo ""  >> rmf_poly_$nuc-$hybrid.tmp
echo "	return result" >> rmf_poly_$nuc-$hybrid.tmp

cat rmf_poly_$nuc-$hybrid.tmp >> ../../$outfile.pyx

### INIT ###

#cat x0.tmp >> ../../$outfile.pyx


### JACOBIAN ###

org_file="Nucleus_"$nuc"_JacobianKyle.txt"
file="Nucleus_"$nuc"_Jacobian.tmp"
cp $org_file $file

echo "
@cython.boundscheck(False)
@cython.wraparound(False)
def rmf_poly_"$nuc"_"$hybrid"_jac(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params):

	cdef np.ndarray[DTYPE_t,ndim=2] result = np.zeros([x0.shape[0],x0.shape[0]], dtype=DTYPE)
" > rmf_poly-$nuc.tmp

cp ../../parse_math.sh .
./parse_math.sh $nbasis $nfields $nstates $pstates $file $basisfile

echo "" >> $file
length=$(wc -l < $file)
totbasis=$(awk -v len=$length 'BEGIN{print sqrt(len)}')
echo $nuc" basis: $totbasis"
rm py.tmp

i=0
while read l; do
  x=$((i / totbasis))
  y=$((i % totbasis))
  echo -e "\tresult[$y,$x]=$l" >> py.tmp
  i=$((i+1))
done < $file

cat rmf_poly-$nuc.tmp py.tmp > rmf_poly_$nuc-$hybrid.tmp

sed -ri '/^.{,23}$/d' rmf_poly_$nuc-$hybrid.tmp


echo ""  >> rmf_poly_$nuc-$hybrid.tmp
echo "	return result" >> rmf_poly_$nuc-$hybrid.tmp

cat rmf_poly_$nuc-$hybrid.tmp >> ../../$outfile.pyx

### PROTON RADIUS FUNCTION ###

org_file="Nucleus_"$nuc"_ProtonRadius.txt"
file="Nucleus_"$nuc"_ProtonRadius.tmp"
cp $org_file $file

echo "
@cython.boundscheck(False)
@cython.wraparound(False)
def proton_radius_"$nuc"_"$hybrid"(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params):

	cdef double radius
" > rmf_poly-$nuc.tmp

cp ../../parse_math.sh .
./parse_math.sh $nbasis $nfields $nstates $pstates $file $basisfile

rm py.tmp

echo "" >> $file
i=0
while read l; do
  echo -e "\tradius=$l" >> py.tmp
  i=$((i+1))
done < $file

cat rmf_poly-$nuc.tmp py.tmp > rmf_poly_$nuc-$hybrid.tmp

echo ""  >> rmf_poly_$nuc-$hybrid.tmp
echo "	return radius" >> rmf_poly_$nuc-$hybrid.tmp

cat rmf_poly_$nuc-$hybrid.tmp >> ../../$outfile.pyx

### TOTAL ENERGY FUNCTION ###

org_file="Nucleus_"$nuc"_TotalEnergyNucleons.txt"
file="Nucleus_"$nuc"_TotalEnergyNucleons.tmp"
cp $org_file $file

echo "
@cython.boundscheck(False)
@cython.wraparound(False)
def nucleon_energy_"$nuc"_"$hybrid"(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params):

	cdef double energy
" > rmf_poly-$nuc.tmp

cp ../../parse_math.sh .
./parse_math.sh $nbasis $nfields $nstates $pstates $file $basisfile

rm py.tmp

echo "" >> $file
i=0
while read l; do
  echo -e "\tenergy=$l" >> py.tmp
  i=$((i+1))
done < $file

cat rmf_poly-$nuc.tmp py.tmp > rmf_poly_$nuc-$hybrid.tmp

echo ""  >> rmf_poly_$nuc-$hybrid.tmp
echo "	return energy" >> rmf_poly_$nuc-$hybrid.tmp

cat rmf_poly_$nuc-$hybrid.tmp >> ../../$outfile.pyx

### TOTAL FIELDS ENERGY FUNCTION ###

org_file="Nucleus_"$nuc"_TotalEnergyFields.txt"
file="Nucleus_"$nuc"_TotalEnergyFields.tmp"
cp $org_file $file

echo "
@cython.boundscheck(False)
@cython.wraparound(False)
def field_energy_"$nuc"_"$hybrid"(np.ndarray[DTYPE_t] x0, np.ndarray[DTYPE_t] params):

	cdef double energy
" > rmf_poly-$nuc.tmp

cp ../../parse_math.sh .
./parse_math.sh $nbasis $nfields $nstates $pstates $file $basisfile

rm py.tmp

echo "" >> $file
i=0
while read l; do
  echo -e "\tenergy=$l" >> py.tmp
  i=$((i+1))
done < $file

cat rmf_poly-$nuc.tmp py.tmp > rmf_poly_$nuc-$hybrid.tmp

echo ""  >> rmf_poly_$nuc-$hybrid.tmp
echo "	return energy" >> rmf_poly_$nuc-$hybrid.tmp

cat rmf_poly_$nuc-$hybrid.tmp >> ../../$outfile.pyx

### FORTRAN ###

echo -e "function rmf_poly_"$nuc"_$hybrid(x0, params) result(result)
    
    implicit none

    !integer, intent(in) :: n

    real(kind=8), intent(in) :: x0(:)
    real(kind=8), intent(in) :: params(:)
    real(kind=8) :: result(size(x0))


" > rmf_poly_f90.tmp

org_file="Nucleus_"$nuc"_CoeffEquations.txt"
file="Nucleus_"$nuc"_CoeffEquations.tmp"
cp $org_file $file

./parse_math.sh $nbasis $nfields $nstates $pstates $file $basisfile

echo "" >> $file
i=0
while read l; do
  echo -e "\tresult[$i]=$l" >> f90.tmp
  i=$((i+1))
done < $file


for ele in `seq 0 $i`; do

sed -i "s/result\[$ele\]/result($((ele+1)))/g" f90.tmp
sed -i "s/x0\[$ele\]/x0($((ele+1)))/g" f90.tmp
sed -i "s/params\[$ele\]/params($((ele+1)))/g" f90.tmp

done

echo "end function rmf_poly_"$nuc"_"$hybrid >> f90.tmp

cat rmf_poly_f90.tmp f90.tmp >> ../../$outfile.f90

rm *.tmp

cd ..

done

rm *.tmp

cd ..

python3 setup.py build_ext --inplace

mv $outfile.*.so ../$outfile.so

f2py3 -c --f90flags="-march=native -mtune=native -Ofast -ffree-line-length-0" -m $outfile_fort $outfile.f90

mv $outfile_fort.*.so ../$outfile_fort.so