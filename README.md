# Introduction
Repository for GenX code, inputs, and results for the paper "Robust capacity expansion planning in hydro-dominated power systems: A Nordic case study"

## Contents
- GenX - Modified version of GenX 0.3.3 with additional features for handling finite fuel supplies, model flows between zones at the same time resolution as other input values, constraints on the maximum instantaneous consumption of flexible loads and hydro reservoir level constraints.
- input - Three GenX inputs for the weather years 1991, 2002, and 2008.
- results - Results from running GenX on the three input files.


## Required software versions:
- Julia: 1.8.5
- Gurobi: 9.5.2

## Julia packages:
```julia
[6e4b80f9] BenchmarkTools v1.5.0
[336ed68f] CSV v0.10.14
[9961bab8] Cbc v1.2.0
[e2554f3b] Clp v1.1.0
[aaaa29a8] Clustering v0.14.4
[861a8166] Combinatorics v1.0.2
[a93c6f00] DataFrames v1.6.1
[864edb3b] DataStructures v0.18.20
[b4f34e82] Distances v0.10.11
[2e9cd046] Gurobi v1.3.0
[87dc4568] HiGHS v1.1.4
[4076af6c] JuMP v1.23.0
[b8f27783] MathOptInterface v1.31.1
[731186ca] RecursiveArrayTools v2.38.10
[82193955] SCIP v0.11.14
[2913bbd2] StatsBase v0.33.21
[ddb6d928] YAML v0.4.12
[ade2ca70] Dates
[37e2e46d] LinearAlgebra
[9a3f8284] Random
[10745b16] Statistics
```
