# Load Dependencies
using AlgebraOfGraphics, BenchmarkTools, CairoMakie, CSV, DataFrames, Polynomials, Revise

# Type and Function Definitions
includet("definitions.jl")

# Read Data
op_file = "OperatingStress.csv"
res_file = "ResidualStress.csv"
(s_op, σ_op, s_op_units, σ_op_units) = read_fracture_stress(op_file)
(s_res, σ_res, s_res_units, σ_res_units) = read_fracture_stress(res_file)

# Compute Polynomial Fits Individually
# (This method still works okay on 1.10.0-beta1.)
op_results = polyfit(s_op, σ_op)
res_results = polyfit(s_res, σ_res)

# Compute Polynomails Fits Jointly
# (This method is unusably slow on 1.10.0-beta1 and August 8th nightly.)
op_results, res_results = polyfit([s_op, s_res], [σ_op, σ_res])

# Plot Results
op_plot = plotpolyfit(op_results; x_units=s_op_units, y_units=σ_op_units)
res_plot = plotpolyfit(res_results; x_units=s_res_units, y_units=σ_res_units)
