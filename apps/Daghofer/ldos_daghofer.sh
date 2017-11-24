#!/bin/bash
#
#xi0=4.0
#xi1=4.0
#
#xi0=8.0
#xi1=8.0
#
#pairing_d=0.0
#pairing_s=0.025
#pairing_p=0.0
#
#  outfile=ldos_daghofer_p_${xi0}_${xi1}_${pairing_d}_${pairing_s}_${pairing_p}.jld
#  export OMP_NUM_THREADS=1
#  export MKL_NUM_THREADS=1
#
#  if [ ! -e "$outfile" ] ; then
#  julia ldos_daghofer.jl --nx 128 --ny 128 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile &
#  else
#      echo "file $outfile exists"
#  fi
#
#
#xi0=2.0
#xi1=2.0
#
#pairing_d=0.0
#pairing_s=0.025
#pairing_p=0.0
#
#  outfile=ldos_daghofer_p_${xi0}_${xi1}_${pairing_d}_${pairing_s}_${pairing_p}.jld
#  export OMP_NUM_THREADS=1
#  export MKL_NUM_THREADS=1
#
#  if [ ! -e "$outfile" ] ; then
#  julia ldos_daghofer.jl --nx 128 --ny 128 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile &
#  else
#      echo "file $outfile exists"
#  fi
#
#xi0=1.0
#xi1=1.0
#
#pairing_d=0.0
#pairing_s=0.025
#pairing_p=0.0
#
#  outfile=ldos_daghofer_p_${xi0}_${xi1}_${pairing_d}_${pairing_s}_${pairing_p}.jld
#  export OMP_NUM_THREADS=1
#  export MKL_NUM_THREADS=1
#
#  if [ ! -e "$outfile" ] ; then
#  julia ldos_daghofer.jl --nx 128 --ny 128 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile &
#  else
#      echo "file $outfile exists"
#  fi

xi0=0.1
xi1=0.1

pairing_d=0.05
pairing_s=0.025
for pairing_p in 0.0 1.0 2.0
do

    outfile=ldos_daghofer_v2_p_${xi0}_${xi1}_${pairing_d}_${pairing_s}_${pairing_p}.jld
    export OMP_NUM_THREADS=2
    export MKL_NUM_THREADS=2

    if [ ! -e "$outfile" ] ; then
        julia ldos_daghofer.jl --nx 128 --ny 128 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile &
    else
        echo "file $outfile exists"
    fi
done
