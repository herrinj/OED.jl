#!/bin/bash
for ip in {0..3}
do
    julia src/mainOEDtypeA.jl --matfile data/OED_pent.mat --exName results/OED_A_pent/OED_A_pent --maxIter 20 --alphaMin 1e-1 --alphaMax 1e2 --nAlpha 20 --IPmode $ip
done
