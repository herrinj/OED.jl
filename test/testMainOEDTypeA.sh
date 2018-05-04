#!/bin/bash
for ip in {0..3}
do
    julia -p 2 src/mainOEDtypeA.jl --matfile data/OED_rect.mat --exName /tmp/OED_A_rect --maxIter 2 --alphaMin 1e0 --alphaMax 1e2 --nAlpha 2 --IPmode $ip --nex 2
done
