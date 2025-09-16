# Simulation of Population Dynamics with Memory-Dependent Resistance

This Julia script simulates a two-population model (sensitive and resistant) with a memory-dependent conversion rate and a metabolic cost for resistance. It uses the `FdeSolver` package to solve fractional differential equations and `Plots` with the GR backend for visualization.

## Overview
- **Purpose**: Models the dynamics of sensitive (N0) and resistant (N1) populations under periodic pulses, incorporating a memory-dependent conversion rate (`α(T)`) and a metabolic cost (`c`) that affects the resistant population's growth rate.
- **Key Features**:
  - Simulates population growth with logistic dynamics and periodic kill pulses.
  - Includes a memory effect where the conversion rate depends on the number of past pulses.
  - Accounts for a cumulative metabolic cost reducing the resistant population's growth rate.
- **Output**: Generates a plot (`Memory_rise_and_collapse.png`) showing the abundance of sensitive, resistant, and total populations over time, with shaded pulse intervals.

## Requirements
- Julia (tested with v1.9+)
- Packages: `FdeSolver`, `Plots`

## Usage
1. Install required packages: `using Pkg; Pkg.add(["FdeSolver", "Plots"])`
2. Run the script to solve the model and generate the plot.
3. Check the output file `Memory_rise_and_collapse.png` for results.

## Parameters
- `tSpan`: Simulation time range (0 to 130).
- `h`: Step size (0.01).
- `μ`: Fractional order parameters [1, 1].
- `X0`: Initial conditions [0.2, 0.0] for sensitive and resistant populations.
- `par`: Model parameters [r0, r1, k, pulse, period, c].