#!/bin/bash
for ip in {0..3}
do
    julia src/mainOEDtypeB.jl --matfile data/OED_rect.mat --exName results/OED_B_rect/OED_B_rect --maxIter 20 --lambda 1e-1 --IPmode $ip --OEDtype 2 --l 2
done
