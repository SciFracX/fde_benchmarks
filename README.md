# Fractional Differential Equation Solver Benchmarks

This repository contains a cross-language benchmark suite for fractional ordinary differential equations (FODEs). It compares runtime and numerical error across Julia, MATLAB, and Python implementations on several representative problems:

- Single-term linear FODE
- Single-term nonlinear FODE
- Stiff FODE system
- Multi-term FODE

The benchmark outputs are stored as CSV files and then aggregated into publication-style performance profiles (time vs error) using plotting scripts.

## What Is Benchmarked

Each benchmark sweeps the same step-size family:

- dt = 2^(-3), 2^(-4), 2^(-5), 2^(-6), 2^(-7)

For every solver and dt:

- Runtime is measured (wall-clock time)
- Error is computed against an analytic or manufactured reference solution
- One CSV file is generated per solver

Most CSV files follow the two-column convention:

- time: execution time in seconds
- error: numerical error norm

## Repository Layout

Core benchmark scripts are under benchmarks/benchmark_fode:

- benchmark_scipt_Julia: Julia benchmark generators
- benchmark_scipt_MATLAB: MATLAB benchmark generators and solver functions
- benchmark_scipt_Python: Python benchmark generators (pycaputo)

Generated benchmark data is under:

- benchmarks/data

Cross-language plotting and comparison scripts are under:

- benchmarks/fode_wpd_*.jl

## Detailed File Guide

### 1) Julia benchmark generators

Location: benchmarks/benchmark_fode/benchmark_scipt_Julia

- xue_benchmarks.jl
	- Noncommensurate 3D nonlinear system benchmark
	- Exports: FdeSolver_PECE, FractionalDiffEq_PECE, FractionalDiffEq_PITrap, FractionalDiffEq_PIRect
- linear_singleterm_benchmarks.jl
	- Single-term linear model with Mittag-Leffler reference
	- Exports: Linear_Single_* CSV series
- nonlinear_singleterm_benchmarks.jl
	- Single-term nonlinear manufactured-solution benchmark
	- Exports: Single_* CSV series
- stiff_benchmarks.jl
	- Stiff 3D coupled system with manufactured forcing and analytic solution
	- Exports: Stiff_* CSV series
- multiterms_benchmarks.jl
	- Multi-term Caputo benchmark (MTPECE/MTPIEX/MTPITrap/MTPIRect)
	- Exports: FractionalDiffEq_MT* CSV series
- large_benchmarks.jl, fdesovler_on_stiff.jl
	- Additional or exploratory benchmark scripts
- Project.toml
	- Present as local environment marker (currently empty)

### 2) MATLAB benchmark generators and local solvers

Location: benchmarks/benchmark_fode/benchmark_scipt_MATLAB

Main benchmark scripts:

- xue_benchmarks.m
- linear_singleterm_benchmarks.m
- nonlinear_singleterm_benchmarks.m
- stiff_benchmarks.m
- multiterms_benchmarks.m

Local method implementations used by the scripts include:

- fde_pi1_ex.m, fde_pi12_pc.m, fde_pi1_im.m, fde_pi2_im.m
- mt_fde_pi1_ex.m, mt_fde_pi12_pc.m, mt_fde_pi1_im.m, mt_fde_pi2_im.m
- flmm2.m, ml.m
- MFC2demo3c directory and data directory for supporting routines/data

### 3) Python benchmark generators

Location: benchmarks/benchmark_fode/benchmark_scipt_Python

- benchmark.py
	- Noncommensurate 3D nonlinear benchmark using pycaputo PECE
	- Exports: PYCAPUTO_PECE.csv
- linear_singleterm_benchmark.py
	- Linear single-term benchmark (uses FractionalDiffEq.mittleff via python-julia bridge)
	- Exports: Linear_Singleterm_PYCAPUTO_PECE.csv
- nonlinear_singleterm_benchmark.py
	- Nonlinear single-term manufactured benchmark
	- Exports: Singleterm_PYCAPUTO_PECE.csv

## Plotting and Cross-Framework Comparison Scripts

Location: benchmarks

- fode_wpd_nonstiff.jl: non-stiff noncommensurate comparison
- fode_wpd_stiff.jl: stiff comparison
- fode_wpd_singleterm.jl: single-term comparison
- fode_wpd_linear_singleterm.jl: linear + nonlinear single-term comparison panels
- fode_wpd_mtnonstiff.jl: multi-term comparison
- fode_wpd_nonstiff_and_stiff.jl: combined non-stiff and stiff figure
- fdde_wpd_nonstiff.jl: additional non-stiff figure variant

These scripts read CSV files from benchmarks/data and produce SVG figures (for example: singleterm_fode_benchmarks.svg, stiff_fode_benchmarks.svg, fode_general_benchmarks.svg).

## How To Run

### Julia benchmarks

1. Start Julia in the repository root.
2. Install required packages if needed: FractionalDiffEq, FdeSolver, DataFrames, CSV, BenchmarkTools, SpecialFunctions.
3. Run a generator script, for example:
	 - include("benchmarks/benchmark_fode/benchmark_scipt_Julia/stiff_benchmarks.jl")

### MATLAB benchmarks

1. Add benchmarks/benchmark_fode/benchmark_scipt_MATLAB to the MATLAB path.
2. Run one script, for example:
	 - stiff_benchmarks
3. CSV files are written to benchmarks/data (or current folder for some scripts, depending on path usage).

### Python benchmarks

1. Use Python 3 with numpy, pandas, pycaputo, matplotlib installed.
2. For linear_singleterm_benchmark.py, ensure python-julia bridge is configured and Julia package FractionalDiffEq is available.
3. Run, for example:
	 - python benchmarks/benchmark_fode/benchmark_scipt_Python/nonlinear_singleterm_benchmark.py

## Data Naming Convention

CSV prefixes indicate source and scenario:

- Stiff_*: stiff problem benchmarks
- Single_* or Singleterm_*: nonlinear single-term benchmarks
- Linear_Single_*: linear single-term benchmarks
- FractionalDiffEq_MT* and MATLAB_MT*: multi-term benchmarks
- PYCAPUTO_*: Python pycaputo benchmarks
- MATLAB_*: MATLAB benchmarks

Suffixes indicate solver family (for example PECE, PIEX, PITrap, PIRect, BDF, NewtonGregory, Trapzoid).

## Reproducibility Notes

- Step-size sweep is consistent across languages in the benchmark scripts.
- Error definitions are based on script-level analytic references and may differ in norm details between implementations.
- Runtime is sensitive to machine load and language startup overhead.
- Some scripts use absolute output paths; update paths if the repository is moved.

## Suggested Workflow

1. Run benchmark generators in each language to refresh CSV data.
2. Run Julia plotting scripts in benchmarks to create comparative figures.
3. Use benchmarks/data as the canonical source for tables and plots in reports.