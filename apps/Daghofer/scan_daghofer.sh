#!/bin/bash

pairing_d=0.05
pairing_s=0.025
pairing_p="0.0 0.1 1.0 2.0 4.0 8.0 16.0"
#xi0="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0"
#xi1="1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0"
xi0="4.0"
xi1="8.0"
outfile=/data/kmlee/Downloads/daghofer_scan_6.jld

julia scan_daghofer.jl --nx 64 --ny 64 --xi0 $xi0 --xi1 $xi1 --pairing_d $pairing_d --pairing_s $pairing_s --pairing_p $pairing_p --outfile $outfile
