#!/bin/bash

nbasis=$1

nfields=$2

nstates=$3

pstates=$4

file=$5

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

for i in `seq 1 $nfields`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print (i-1)*nbasis + (j-1)}')

sed -i "s/FieldCoeffs0\[$i, $j\]/x0["${ele}"]/g" $file

done
done

for i in `seq 1 $nstates`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + (i-1)*nbasis + (j-1)}')

sed -i "s/NeutronsGCoeffs0\[$i, $j\]/x0[$ele]/g" $file

done
done

for i in `seq 1 $nstates`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + nstates*nbasis + (i-1)*nbasis + (j-1)}')

sed -i "s/NeutronsFCoeffs0\[$i, $j\]/x0[$ele]/g" $file

done
done

for i in `seq 1 $pstates`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + 2*nstates*nbasis + (i-1)*nbasis + (j-1)}')

sed -i "s/ProtonsGCoeffs0\[$i, $j\]/x0[$ele]/g" $file

done
done

for i in `seq 1 $pstates`; do
for j in `seq 1 $nbasis`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + pstates*nbasis + 2*nstates*nbasis + (i-1)*nbasis + (j-1)}')

sed -i "s/ProtonsFCoeffs0\[$i, $j\]/x0[$ele]/g" $file

done
done

for i in `seq 1 $nstates`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + 2*pstates*nbasis + 2*nstates*nbasis + (i-1)}')

sed -i "s/EN\[$i\]/x0[$ele]/g" $file

done

for i in `seq 1 $pstates`; do

ele=$(awk -v i=$i -v j=$j -v nbasis=$nbasis -v nfields=$nfields -v nstates=$nstates -v pstates=$pstates 'BEGIN{print nfields*nbasis + 2*pstates*nbasis + 2*nstates*nbasis + nstates + (i-1)}')

sed -i "s/EP\[$i\]/x0[$ele]/g" $file

done


sed -i 's/==0/ /g' $file
sed -i 's/{//g' $file
sed -i 's/}//g' $file