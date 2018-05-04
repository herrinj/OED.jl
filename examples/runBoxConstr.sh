#!/bin/bash
for a in {1..100..10}
do
    julia mainOED.jl --matfile OED_A_rect.mat --exName OED_A_rect --maxIter 60 --alpha $a --IPmode 3
done
