#!/bin/bash

pairing_d=0.05
pairing_s=0.05
pairing_p=0.1
xi0=4.0
xi1=8.0
outfile=test.jld

julia run_daghofer.jl --nx 16 --ny 16 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile
