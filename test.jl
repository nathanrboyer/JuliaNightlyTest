# Load Dependencies
using Revise                          # for `includet` function
using CSV, DataFrames, Polynomials    # for reading and fitting functions
using AlgebraOfGraphics, CairoMakie   # for plotting funtions

# Type and Function Definitions
includet("definitions.jl")

# Read Data
(x1, y1) = read_data("Dataset1.csv")
(x2, y2) = read_data("Dataset2.csv")

# Compute Polynomial Fits Individually
#   This method still works okay on 1.10.0-beta1.
#   Function call takes roughly 5 minutes.
@time result1 = polyfit(x1, y1)
@time result2 = polyfit(x2, y2)

# Compute Polynomails Fits Jointly
#   This method is unusably slow on 1.10.0-beta1 and August 8th nightly.
#   Function call takes roughly 5 minutes on 1.9.2 and more than 1 day on nightly.
@time result1, result2 = polyfit([x1, x2], [y1, y2])

# Plot Results
plot1 = plotpolyfit(result1; xlabel="x", ylabel="y")
plot2 = plotpolyfit(result2; xlabel="x", ylabel="y")
