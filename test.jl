# Load Dependencies
using Revise                          # for `includet` function
using CSV, DataFrames, Polynomials    # for reading and fitting functions
using AlgebraOfGraphics, CairoMakie   # for plotting funtions

# Type and Function Definitions
includet("definitions.jl")

# Read Data
(x1, y1) = read_data("Dataset1.csv")
(x2, y2) = read_data("Dataset2.csv")

# Trim Data for Speed
#   Can add or remove points as necessary for computational load.
#   `polyfit`` will error and tell you if using too few points.
x1 = x1[1:2:101]
x2 = x2[1:2:101]
y1 = y1[1:2:101]
y2 = y2[1:2:101]

# Compute Polynomial Fits Individually
#   This method still works okay on 1.10.0-beta1.
#   Function takes roughly 5 minutes with full data set.
@time result1 = polyfit(x1, y1)
@time result2 = polyfit(x2, y2)

# Compute Polynomails Fits Jointly
#   This method is unusably slow on 1.10.0-beta1 and August 8th nightly.
#   Function takes roughly 5 minutes on 1.9.2 with full data set.
#   Function takes more than 1 day on nightly with full data set.
@time result1, result2 = polyfit([x1, x2], [y1, y2])

# Plot Results
plot1 = plotpolyfit(result1; xlabel="x", ylabel="y")
plot2 = plotpolyfit(result2; xlabel="x", ylabel="y")
