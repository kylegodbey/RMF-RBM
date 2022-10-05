#!/bin/bash

nbasis=$1

nfields=$2

nstates=$3

pstates=$4

file=$5

basisfile=$6

initfile=buffer.tmp

rm $initfile

sed -i 's/\r$//' $basisfile

readarray -t bases < $basisfile

sed -i 's/\*\^/e/g' $file
sed -i 's/\^/**/g' $file 
sed -i 's/0\. //g' $file
sed -i 's/ms/params[0]/g' $file
sed -i 's/mv/3.9654992946/g' $file
sed -i 's/mρ/3.8666785454/g' $file
sed -i 's/gs/params[1]/g' $file
sed -i 's/gv/params[2]/g' $file
sed -i 's/gρ/params[3]/g' $file
sed -i 's/κ/params[4]/g' $file
sed -i 's/λ/params[5]/g' $file
sed -i 's/ζ/params[6]/g' $file
sed -i 's/Λv/params[7]/g' $file
sed -i 's/Λs/0.0/g' $file
sed -i 's/M/4.758526326458217/g' $file


ele=0

for i in `seq 1 $nfields`; do
nbasis=$(echo  ${bases[0]} | awk -v i=$i '{print $i}')
for j in `seq 1 $nbasis`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print (i-1)*nbasis + (j-1)}')
if [ $j -eq 1 ] ; then

sed -i "s/FieldCoeffs0\[$i, $j\]/(x0["${ele}"] + 1)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
else

sed -i "s/FieldCoeffs0\[$i, $j\]/x0["${ele}"]/g" $file

fi

ele=$((ele+1))

done
done

for i in `seq 1 $nstates`; do
nbasis=$(echo  ${bases[2]} | awk -v i=$i '{print $i}')
for j in `seq 1 $nbasis`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + (i-1)*nbasis + (j-1)}')
if [ $j -eq 1 ] ; then

sed -i "s/NeutronsGCoeffs0\[$i, $j\]/(x0["${ele}"] + 1)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
else
sed -i "s/NeutronsGCoeffs0\[$i, $j\]/x0[$ele]/g" $file
fi
ele=$((ele+1))
done
done

for i in `seq 1 $nstates`; do
nbasis=$(echo  ${bases[1]} | awk -v i=$i '{print $i}')
for j in `seq 1 $nbasis`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + nstates*nbasis + (i-1)*nbasis + (j-1)}')
if [ $j -eq 1 ] ; then

sed -i "s/NeutronsFCoeffs0\[$i, $j\]/(x0["${ele}"] + 1)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
else
sed -i "s/NeutronsFCoeffs0\[$i, $j\]/x0[$ele]/g" $file
fi
ele=$((ele+1))
done
done

for i in `seq 1 $pstates`; do
nbasis=$(echo  ${bases[4]} | awk -v i=$i '{print $i}')
for j in `seq 1 $nbasis`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + 2*nstates*nbasis + (i-1)*nbasis + (j-1)}')
if [ $j -eq 1 ] ; then

sed -i "s/ProtonsGCoeffs0\[$i, $j\]/(x0["${ele}"] + 1)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
else

sed -i "s/ProtonsGCoeffs0\[$i, $j\]/x0[$ele]/g" $file
fi
ele=$((ele+1))
done
done

for i in `seq 1 $pstates`; do
nbasis=$(echo  ${bases[3]} | awk -v i=$i '{print $i}')
for j in `seq 1 $nbasis`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + pstates*nbasis + 2*nstates*nbasis + (i-1)*nbasis + (j-1)}')
if [ $j -eq 1 ] ; then

sed -i "s/ProtonsFCoeffs0\[$i, $j\]/(x0["${ele}"] + 1)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
else

sed -i "s/ProtonsFCoeffs0\[$i, $j\]/x0[$ele]/g" $file
fi
ele=$((ele+1))
done
done

for i in `seq 1 $nstates`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + 2*pstates*nbasis + 2*nstates*nbasis + (i-1)}')

sed -i "s/EN\[$i\]/(x0[$ele]+4.5)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
ele=$((ele+1))
done

for i in `seq 1 $pstates`; do

# ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + 2*pstates*nbasis + 2*nstates*nbasis + nstates + (i-1)}')

sed -i "s/EP\[$i\]/(x0[$ele]+4.5)/g" $file
echo -e "\tx0[$ele]=0.0" >> $initfile
ele=$((ele+1))
done

echo "
@cython.boundscheck(False)
@cython.wraparound(False)
def rmf_poly_48Ca_0_x0():

	cdef np.ndarray[DTYPE_t] x0 = np.zeros([$ele], dtype=DTYPE)
" > header.tmp

cat header.tmp $initfile > x0.tmp

echo -e "\treturn x0" >> x0.tmp

sed -i 's/==0/ /g' $file
sed -i 's/{//g' $file
sed -i 's/}//g' $file