#!/bin/bash

RUN=1

NX=100
NY=100
PAIRING_D=0.1
PAIRING_S=0.1

COUNT=0

for SOC in 0.0 0.05 0.1 ; do
for XI0 in 0.001 1.0 2.0 5.0 10.0 ; do
for XI1 in 0.001 1.0 2.0 5.0 10.0 ; do
for PAIRING_P in 0.001 0.01 0.1 1.0 ; do

    COUNT=$((COUNT+1))
    OUTPUT_FILENAME="`printf 'Run%d-%04d.hdf5' ${RUN} ${COUNT}`"
    echo "${OUTPUT_FILENAME}"

    time julia ../../ldos_spinful_daghofer.jl \
        --nx ${NX} --ny ${NY} \
        --soc ${SOC} --xi0 ${XI0} --xi1 ${XI1} \
        --pairing_d ${PAIRING_D} \
        --pairing_s ${PAIRING_S} \
        --pairing_p ${PAIRING_P} \
        -o "${OUTPUT_FILENAME}"
done
done
done
done
