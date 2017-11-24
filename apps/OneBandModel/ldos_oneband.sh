#!/bin/bash

pairing_d=0.6
pairing_s=0.4
pairing_p=150.0
xi0=8.0
xi1=8.0
outfile=largeldos_nosoc.jld

julia ldos_oneband.jl --nx 256 --ny 256 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile
