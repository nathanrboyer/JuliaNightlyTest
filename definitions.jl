# ---------------------------------------------------------------------------- #
#                                    Types                                     #
# ---------------------------------------------------------------------------- #

"""
    PolyFit(x, y, knots, segmentlengths, polys, RMSE)

Type to hold all data returned from the `polyfit` function.

# Fields
- `x::Vector{Float64}}`: Independent variable
- `y::Vector{Float64}`: Dependent variable
- `knots::Vector{Int}`: Indices in input data to split data for fitting
- `segmentlengths::Vector{Float64}`: Length of each segment computed from x
- `polys::Vector{Polynomial{Float64, :x}}`: Polynomials fit to the data in each segment
- `RMSE::Float64`: Root mean square error for the entire series of polynomial fits
"""
struct PolyFit #TODO: Replace with T<:Number and test with Unitful.jl
    x::Vector{Float64}
    y::Vector{Float64}
    knots::Vector{Int64}
    segmentlengths::Vector{Float64}
    polys::Vector{Polynomial{Float64, :x}}
    RMSE::Float64
end
Base.isless(x1::PolyFit, x2::PolyFit) = Base.isless(x1.RMSE, x2.RMSE) # Define how to sort objects of type PolyFit
Base.show(io::IO, ::MIME"text/plain", x::PolyFit) = print(io, "PolyFit Result Fields: $(join(fieldnames(typeof(x)),", "))")

"""
    FitBias(type=:none, factor=1)

Type to hold inputs for the `compute_fitweights` function.

Defaults `FitBias()` is no bias.
See `compute_fitweights` function help for more information.

# Fields
- `type::Symbol=:none`:
- `factor::Float64=1`:
"""
struct FitBias
    type::Symbol
    factor::Float64
end
FitBias() = FitBias(:none, 1)
Base.show(io::IO, ::MIME"text/plain", x::FitBias) = print(io, "Fit Bias Fields: $(join(fieldnames(typeof(x)),", "))")

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    extractunits(s::String; keep_parentheses::Bool = false)::String

Isolate and return the portion of the input `String` inside paratheses,
usually the units.
"""
function extractunits(s::AbstractString; keep_parentheses::Bool=false)::String
    units = match(r"\(.*\)", s).match
    if !(keep_parentheses)
        units = strip(units, ['(',')'])
    end
    return units
end

"""
    read_fracture_stress(file_location) -> (s, σ, s_units, σ_units)

Reads a .csv file containing normal stress along a path exported from ANSYS.

# Arguments
- `file_location::String`: location of the stress file to be read
- `s::Vector{Float64}`: location data along path from inside to outside of vessel
- `σ::Vector{Float64}`: normal stress data at `s` points along path
- `s_units::String`: units of the `s` column in ANSYS data export
- `σ_units::String`: units of the `σ` column in ANSYS data export
"""
function read_fracture_stress(file_location)
    df = DataFrame(CSV.File(file_location))
    s_column_name = only(names(df, r"S "))
    σ_column_name = only(names(df, r"Normal Stress "))
    s = df[:, s_column_name]
    σ = df[:, σ_column_name]
    s_units = extractunits(s_column_name)
    σ_units = extractunits(σ_column_name)
    return (s, σ, s_units, σ_units)
end

"""
    polyfit(x, y, npolys=4, degree=3; fitselection=1) -> polyfitresult
    polyfit(x, y, knots, degree=3) -> polyfitresult
    polyfit(x, y, knotsets, degree=3; fitselection=1) -> polyfitresult
    polyfit(xsets, ysets, npolys=4, degree=3; fitselection=1, sortby=0) -> polyfitresults
    polyfit(xsets, ysets, knots, degree=3; fitselection=1, sortby=0) -> polyfitresults
    polyfit(xsets, ysets, knotsets, degree=3; fitselection=1, sortby=0) -> polyfitresults

Finds the best polynomial fits to the input data.

# Arguments
- `x::Vector{Float64}`: independent variable of data to be fit, usually position s
- `y::Vector{Float64}`: dependent variable of data to be fit, usually stress σ
- `npolys::Int=4`: number of polynomials to use to fit data. will automatically test all possible `knots` using the `allknots` function
- `degree::Int=3`: degree of the polynomials to use for fitting. the default degree=3 means x^3 will be the highest power used.
- `fitbias::FitBias=FitBias(:none, 1)`: allows point weighting to be applied when computing polynomial fit. See `FitBias` and `compute_fitweights` for more info.
- `fitselection::Int=1`: option to choose knots other than the best fit. 1 means best fit, 2 means second best fit, etc. ranked by lowest RMSE
- `knots::Vector{Int}`: the x-indices where the polynomials start and end. must include the data endpoints e.g. `[1, 5, 11, length(x)]`
- `knotsets::Vector{Vector{Int}}`: vector of `knots` to be tested for best fit
- `xsets::Vector{Vector{Float64}}`: vector of independent variable data, usually `[s, s]`
- `ysets::Vector{Vector{Float64}}`: vector of dependent variable data, usually `[σ_operating, σ_residual]`
- `sortby::Int`: which data in `ysets` to use to compute RMSE for determining "best fit"
    - `sortby=0`: measures the sum of operating and residual RMSE, treating both datasets equally
    - `sortby=1`: best fit is measured against the first `y` data in `ysets`
    - `sortby=2`: best fit is measured against the second `y` data in `ysets`
- `polyfitresult::PolyFit`: polynomial fitting data is collected into a struct PolyFit
- `polyfitresults::Vector{PolyFit}`: a vector of PolyFit data, returned when input data is in sets. `knots` and `segmentlengths` will be the same across the sets.

"""
function polyfit(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, knots::AbstractVector{Int}, degree::Int=3; fitbias::FitBias=FitBias())

    # Check inputs
    if length(x) != length(y)
        throw(ArgumentError("x and y must be the same length"))
    end
    if minimum(diff(knots[:])) < degree
        throw(ArgumentError("knots are spaced too closely for the given degree. reduce degree or increase knot spacing. knots must be sorted from smallest to largest."))
    end
    if knots[1] != 1
        throw(ArgumentError("first knot location must always be 1."))
    end
    if knots[end] != length(x)
        throw(ArgumentError("last knot location must always be the number of elements in x and y."))
    end

    # Fit polynomials
    nknots=length(knots) # Number of knots
    npolys=nknots-1      # Number of polynomials
    n_points_to_fit=length(x)+npolys-1  # Total number of points for all polynomials including overlap
    RSS=0.0              # Residual sum of squares initialization
    polys = Vector{Polynomial{Float64, :x}}(undef, npolys)
    for n=1:npolys
        x_segment = x[knots[n]:knots[n+1]]
        y_segment = y[knots[n]:knots[n+1]]
        fitweights = compute_fitweights(length(x_segment), fitbias)
        polys[n]=fit(x_segment, y_segment, degree; weights=fitweights)
        RSS+=sum((polys[n].(x_segment) .- y_segment) .^ 2)
    end
    RMSE=sqrt(RSS/n_points_to_fit) # Root mean squared error

    # Determine length of each polynomial segment
    segmentlengths = diff(x[knots])

    # Combine data into a struct
    polyfitresult=PolyFit(x, y, knots, segmentlengths, polys, RMSE)

    return polyfitresult
end
function polyfit(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, knotsets::AbstractVector{<:AbstractVector{Int}}, degree::Int=3; fitbias::FitBias=FitBias(), fitselection::Int=1)
    polyfitresults=Array{PolyFit}(undef, size(knotsets))
    for (n, knots) in pairs(knotsets)
        polyfitresults[n] = polyfit(x, y, knots, degree; fitbias)
    end
    polyfitresult=partialsort!(polyfitresults, fitselection) #Sorts by polyfitresult.RMSE
    return polyfitresult
end
function polyfit(x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, npolys::Int=4, degree::Int=3; fitbias::FitBias=FitBias(), fitselection::Int=1)

    #Validate inputs
    if npolys<1
        throw(ArgumentError("npolys must be a positive integer"))
    end
    if degree<0
        throw(ArgumentError("degree must be a positive integer or zero"))
    end

    # Function Inputs
    data_length=length(x)        # Total number of data points
    min_segment_length=degree+1  # Minimum number of points in any segment
    nknots=npolys+1              # Number of knots

    # Compute all possible knot combinations
    knotsets=allknots(data_length, min_segment_length, nknots)

    # Fit polynomials
    polyfitresult = polyfit(x, y, knotsets, degree; fitbias, fitselection)

    return polyfitresult
end
function polyfit(xsets::AbstractArray{<:AbstractVector{<:Number}}, ysets::AbstractArray{<:AbstractVector{<:Number}}, knotsets::AbstractArray{<:AbstractVector{Int}}, degree::Int=3; fitbias::FitBias=FitBias(), fitselection::Int=1, sortby::Int=0)

    #Check inputs
    if (length(xsets) != length(ysets)) throw(ArgumentError("number of datasets in x and y must match")) end
    if (any(length.(xsets) .!= length(first(xsets))) || any(length.(ysets) .!= length(first(ysets))) || length(first(xsets)) != length(first(ysets))) throw(ArgumentError("datasets must be the same length")) end
    if (any(sum.([abs.(x .- first(xsets)) for x in xsets]) .> 1e-8)) throw(ArgumentError("x must be the same for all datasets")) end

    # Fit polynomials
    polyfitresultsets=[Array{PolyFit}(undef, size(xsets)) for n in knotsets]
    for (n, knots) in pairs(knotsets)
        for m in eachindex(xsets)
            polyfitresultsets[n][m] = polyfit(xsets[m], ysets[m], knots, degree; fitbias)
        end
    end

    # Sort polynomials fits
    if sortby in eachindex(xsets)
        permvector = partialsortperm(getfield.(getindex.(polyfitresultsets, sortby), :RMSE), fitselection) # Sorts by selected dataset's polyfitresultsets.RMSE
    else
        RSSsum = sum.([getfield.(polyfitresultsets[n], :RMSE) .^ 2 .* length.(getfield.(polyfitresultsets[n], :x)) for n in 1:length(knotsets)]) # Sorts by sum of datasets polyfitresultsets.RMSE (RMSE cannot be added directly)
        permvector = partialsortperm(RSSsum, fitselection)
    end
    polyfitresults = polyfitresultsets[permvector] # Sort polyfitdataset according to permutation vector and select fit (defaults to best fit)

    return polyfitresults # Array{PolyFit} the same size as xsets
end
polyfit(xsets::AbstractArray{<:AbstractVector{<:Number}}, ysets::AbstractArray{<:AbstractVector{<:Number}}, knots::AbstractVector{Int}, degree::Int=3; fitbias::FitBias=FitBias(), fitselection::Int=1, sortby::Int=0) = polyfit(xsets, ysets, [knots], degree; fitbias, fitselection, sortby)
function polyfit(xsets::AbstractArray{<:AbstractVector{<:Number}}, ysets::AbstractArray{<:AbstractVector{<:Number}}, npolys::Int=4, degree::Int=3; fitbias::FitBias=FitBias(), fitselection::Int=1, sortby::Int=0)

    #Validate inputs
    if npolys<1
        throw(ArgumentError("npolys must be a positive integer"))
    end
    if degree<0
        throw(ArgumentError("degree must be a positive integer or zero"))
    end

    # Function Inputs
    data_length=length(first(xsets)) # Total number of data points
    min_segment_length=degree+1      # Minimum number of points in any segment
    nknots=npolys+1                  # Number of knots

    # Compute all possible knot combinations
    knotsets=allknots(data_length, min_segment_length, nknots)

    # Fit polynomials
    polyfitresults = polyfit(xsets, ysets, knotsets, degree; fitbias, fitselection, sortby)

    return polyfitresults # Array{PolyFit} the same size as xsets
end
export polyfit

"""
    compute_fitweights(num_points, fitbias) -> fitweights
    compute_fitweights(num_points, fitbiastype, fitbiasfactor) -> fitweights

Computes the weights to feed into the `fit` function from the Polynomials.jl package.

# Arguments
- `num_points::Int`: number of points in the data to be fit
- `fitbias::FitBias`: struct containing the fit bias type and factor
- `fitbiastype::Symbol=:none`: type of bias to use when computing weights (Note: points away from the biased area may be far from the fit curve)
    - `:none`: equal weighting to all points
    - `:left`: first points in the dataset are weighted for closer fit
    - `:right`: last points in the dataset are weighted for closer fit
    - `:ends`: points near both ends of the dataset are weighted for closer fit
    - `:center`: points near the center of the dataset are weighted for closer fit
- `fitbiasfactor::Int=1`: exponentiation factor to apply to the linear weights (causes more severe bias)
- `fitweights::Vector{Int}`: output vector of weights with length `num_points` (to pass to Polynomials.jl `fit` function)
"""
function compute_fitweights(num_points::Int, fitbiastype::Symbol, fitbiasfactor::Real)
    N = BigInt(num_points) # BigInt conversion needed to prevent integer overflow upon exponentiation.
    type = fitbiastype
    factor = fitbiasfactor
    if type === :none
        fitweights = ones(N)
    elseif type === :left
        fitweights = reverse(1:N) .^ factor
    elseif type === :right
        fitweights = (1:N) .^ factor
    elseif type === :ends
        fitweights = vcat(reverse(1:2:N), 1:2:N) .^ factor
    elseif type === :center
        fitweights = vcat(1:2:N, reverse(1:2:N)) .^ factor
    else
        throw(ArgumentError("The specified `fitbiastype` has not been defined. Choose :none, :left, :right, :ends, or :center."))
    end
    return fitweights
end
compute_fitweights(num_points::Int, fitbias::FitBias) = compute_fitweights(num_points, fitbias.type, fitbias.factor)
export compute_fitweights

"""
    allknots(data_length::Int, min_segment_length::Int, nknots::Int) -> knotsets::Vector{Vector{Int}}

Computes all possible knot locations on a dataset.

# Arguments
- `data_length`: total number of points in the dataset
- `min_segment_length`: minimum segment length between knots
- `nknots`: total number of knots
- `knotsets`: vector collection of all possible knots
"""
function allknots(data_length::Int, min_segment_length::Int, nknots::Int)

    # Initialize first knot combination
    knots = [ones(Int, nknots-1)..., data_length]
    for k = 2:nknots-1
        knots[k] = knots[k-1] + min_segment_length - 1
    end
    knotsets=[copy(knots)]

    # Ensure first knot combination is possible
    if (nknots < 2) || (knots[nknots] - knots[nknots-1] < min_segment_length)
        throw(ArgumentError("no possible knot sets with the given parameters. add more datapoints or remove some knots"))
    end

    # Determine all remaining knot combinations. k tracks the current knot inside knots. i tracks the value of the current knot somewhere in 1:data_length.
    function incrementknots!(knotsets, k=2)
        imax = data_length - (min_segment_length - 1) * (nknots - k) # Maximum i location in dataset that the current knot can take
        while knots[k] < imax                                        # While current knot not yet at maximum location
            if k < nknots-1                                          # If current knot is not the last mobile knot
                incrementknots!(knotsets, k+1)                       # Change current knot to next knot
            end
            knots[k] += 1                                            # Increment current knot
            for kk = k+1 : nknots-1                                  # For all mobile knots downstream of current knot
                knots[kk] = knots[kk-1] + min_segment_length - 1     # Reset downstream knots to minimum location relative to current knot
            end
            push!(knotsets, copy(knots))                             # Add current knot to knot set
        end
        return knotsets
    end

    incrementknots!(knotsets)
    return knotsets
end
export allknots

"""
    plotpolyfit(polydata::PolyFit)

Plots `PolyFit` data and returns the plot as an `AlgebraOfGraphics` `Layers` object.
"""
function plotpolyfit(polydata::PolyFit; xlabel="r", ylabel="σ", x_units=nothing, y_units=nothing)
    # Extract Variables
    x = polydata.x
    y = polydata.y
    knots = polydata.knots
    polys = polydata.polys

    # Plot Labels
    full_x_label = x_units === nothing ? xlabel : xlabel * " (" * x_units * ")"
    full_y_label = y_units === nothing ? ylabel : ylabel * " (" * y_units * ")"
    polytitle = "Polynomial Fits"
    polyentry = "Region"
    datatitle = "FEA Data"
    dataentry = "Stress Values"

    # Create DataFrame
    df = DataFrame(; x, y)
    temp_polyfit = eltype(y)[]
    temp_polynumber = String[]
    for (i, p) in pairs(polys)
        indexrange = knots[i]:knots[i+1]
        append!(temp_polyfit, p.(df.x[indexrange]))
        append!(temp_polynumber, fill(polyentry * " $i", indexrange))
    end
    for i in 2:length(knots)-1
        index = knots[i] + i - 2
        insert!(df, index, df[index,:])
    end
    df.y_label = fill(dataentry, nrow(df))
    df.y_polyfit = temp_polyfit
    df.y_polyfit_label = temp_polynumber

    # Plot Data
    plot = data(df) * (mapping(:x => full_x_label, :y_polyfit => full_y_label, color = :y_polyfit_label => polytitle) * visual(Lines, linewidth=15)
                        + mapping(:x => full_x_label, :y => full_y_label, label = :y_label => datatitle) * visual(Scatter))
    draw(plot) |> display
    return plot
end
