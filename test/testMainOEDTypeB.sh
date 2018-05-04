#!/bin/bash
for ip in {0..3}
do
    julia -p 2 src/mainOEDtypeB.jl --matfile data/OED_rect.mat --exName /tmp/OED_B_rect --maxIter 2 --lambda 1e-2 --IPmode $ip --OEDtype 2 --nex 5
done
