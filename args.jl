
const defFloat = Float64
const defInt = Int64
const Δt = defFloat( 1e-1 )

using Plots
using DelimitedFiles
using Measures
using LaTeXStrings

using StatsBase
using LinearAlgebra

# Predetermined file folders.
# For saving figures.
figurefolder = "~/Documents/papers/ant-spatial/figures/"

# Define a custom color palette and figure defaults.
colorlist = [
    :indianred,       # CD5C5C
    :cornflowerblue,  # 6699FF
    :darkorchid,      # 9932CC
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
    :indianred,       # CD5C5C
    :cornflowerblue,  # 6699FF
    :darkorchid,      # 9932CC
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
]
linestylelist = [:solid, :solid, :solid, :solid, :solid, :dash, :dash, :dash, :dash, :dash]
default(
    guidefont=font(10, "Computer Modern"),  # axis labels
    tickfont=font(10, "Computer Modern"),   # tick labels
    legendfont=font(10, "Computer Modern"), # legend
    titlefont=font(10, "Computer Modern"),  # title
    palette=colorlist,
    left_margin=5pt, right_margin=5pt, bottom_margin=5pt, top_margin=5pt,
)

function saveplot(plt, filename::String; background=:transparent)
    plot!( plt; background=background, legend_background_color=:white )
    savefig( plt,  filename )
    plot!( plt; background=:white )
    return plt
end

function readexpdata(foldername::String, filelist::Vector{String}; w::defInt=2
    )::Vector{Vector{defFloat}}
    adata = [readdlm( foldername*file, ',', defFloat )[1,:] for file ∈ filelist]
    return [sma( length( alist ), w, alist ) for alist ∈ adata]
end

function savemisshapendata(filename::String, xdata)
    open(filename, "w") do io
        for sublist in xdata
            println(io, join(string.(sublist), ","))
        end
    end
    return filename
end

function readmisshapendata(dtype::DataType, filename::String)
    data = Vector{Vector{dtype}}()
    open( filename, "r" ) do io
        for line in eachline(io)
            values = split(line, ",")
            parsed = map(x -> parse(dtype, x), values)
            push!(data, parsed)
        end
    end
    return data
end

function sma(T::Int, w::Int, alist::Vector{Float64})
    pad = div(w, 2)
    return [mean( alist[max( 1, i - pad ):min( T, i + pad )] ) for i ∈ 1:T]
end

function fitpower(x::Vector{defFloat}, y::Vector{defFloat})::Vector{defFloat}
    X = hcat( ones( length( x ) ), log10.( x ) )
    return X \ log10.( y )
end

function fitconfint(x::Vector{Float64}, y::Vector{Float64}, b0::Float64, φ::Float64; α=0.05)
    X = hcat( ones( length( x ) ), log10.( x ) )
    Y = log10.( y )
    β = [b0, φ]

    # Residuals
    r = Y - X*β
    n, p = size( X )
    σ̂2 = sum( r.^2 ) / (n - p)

    # Covariance matrix of β
    covβ = σ̂2 * inv( transpose( X ) * X )
    se = sqrt.( diag( covβ ) )

    # t-value for confidence interval
    tval = quantile( TDist(n - p), 1 - α/2 )

    # Confidence intervals
    ci = hcat( β .- tval .* se, β .+ tval .* se )
    return ci
end
