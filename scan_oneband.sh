#!/bin/bash

pairing_d=0.6
pairing_s=0.4
pairing_p="0.0 150.0"
xi0=8.0
xi1=8.0
outfile=testscan.jld

julia scan_oneband.jl --nx 128 --ny 128 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile
