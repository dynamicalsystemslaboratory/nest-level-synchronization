
include( "vars.jl" )

using StaticArrays

# Pair-wise distance function.
function pairwisedist(L::defFloat, xlist::Vector{defFloat}, ylist::Vector{defFloat})::Vector{defFloat}
    # Initialize distance matrix variables.
    N = length( xlist )
    D = Vector{defFloat}( undef, round( defInt, N*(N-1)/2 ) )

    # Iterate through, calculating distance.
    for i ∈ 1:N-1
        for j ∈ i+1:N
            # Horizontal and vertical difference (with periodic boundary).
            dx = xlist[i] - xlist[j];  dx -= round( dx/L ) * L
            dy = ylist[i] - ylist[j];  dy -= round( dy/L ) * L

            # Distance.
            k = round( defInt, (i - 1)*(2N - i)/2 ) + (j - i)
            D[k] = √(dx*dx + dy*dy)
        end
    end

    # Return distance matrix.
    return D
end

# Outline the movement duration.
function createwedge(r::defFloat, α::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θ::defFloat;
    N::defInt=5)::Vector{SVector{2,defFloat}}
    # Line 1: From start to "left" boundary.
    ℓ1 = SVector{2,defFloat}.( [[xlist[1], ylist[1]],
        [xlist[1] + r*cos(θ - α), ylist[1] + r*sin(θ - α)],
        [xlist[end] + r*cos(θ - α), ylist[end] + r*sin(θ - α)]] )

    # Line 2: Field of view arc into N points.
    δθlist = range( θ - α*(1 - 1/N), θ + α*(1 - 1/N), N-2 )
    ℓ2 = [[xlist[end] + r*cos( δθ ), ylist[end] + r*sin( δθ )] for δθ in δθlist]

    # Line 3: From "right" boundary to start.
    ℓ3 = [[xlist[end] + r*cos( θ + α ),ylist[end] + r*sin( θ + α )],
        [xlist[1] + r*cos( θ + α ), ylist[1] + r*sin( θ + α )],
        [xlist[1], ylist[1]]]

    # Return closed polygon object.
    return [ℓ1..., ℓ2..., ℓ3...]
end

# Compute area of a simple polygon using the Shoelace formula.
function shoelace(plist::Vector{SVector{2,Float64}})::defFloat
    N = length( plist )

    s = defFloat( 0.0 )
    for i in 1:N
        j = i == N ? 1 : i + 1
        s += plist[i][1]*plist[j][2] - plist[j][1]*plist[i][2]
    end

    return abs( s )/2
end

# Compute the unique area coverage.
function visioncoverage(r::defFloat, α::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θlist::Vector{defFloat};
    verbose::Bool=false)
    # Extract points in time when the orientation was updated.
    tlist = findall( [(abs.( diff( θlist ) ) .> 1e-6)..., true] )

    # Collect field of vision polygons from state.
    pdata = [createwedge( r, α, xlist[t1:t2], ylist[t1:t2], θlist[t1+1] )
        for (t1, t2) ∈ zip( tlist[1:end-1], tlist[2:end] )]

    # Compute cumulative area over the activity stint.
    A = sum( shoelace.( pdata ) )  # TODO: Not accounting for holes/intersections.
    return verbose ? (A, pdata) : A
end

# Compute area per duration of movement in a given time-series.
function durationarea(zlist::Vector{State})::Vector{defFloat}
    # Unpack state variable to use with geometry functions.
    xlist = [z.x[1] for z ∈ zlist]
    ylist = [z.y[1] for z ∈ zlist]
    θlist = [z.θ[1] for z ∈ zlist]
    φlist = [z.φ[1] for z ∈ zlist]
    ψlist = [z.ψ[1] for z ∈ zlist]

    # Extract movement durations as individual lists.
    tlist = findall( ψlist .== 1 )
    t̂list = findall( φlist .== 1 )
    wdatalist = Vector{Any}()
    for (t1, t2) ∈ zip( [1; tlist[1:end-1]], tlist[1:end] )
        any( [t1 ≤ t̂ ≤ t2 for t̂ ∈ t̂list] ) || push!(
            wdatalist, (xlist[t1:t2-2] .- L/2, ylist[t1:t2-2] .- L/2, θlist[t1:t2-2]) )
    end

    # Compute area coverage numerically.
    return [visioncoverage( r, α, wdata... ) for wdata ∈ wdatalist]
end

# Critical spontaneous deactivation rate.
function criticalρ(params::Nondim, N::defInt, μ::defFloat)::defFloat
    num = 2*μ*params.β*sin( params.α )
    den = 1/(N*params.η) + params.τ
    return √(num/den)
end

# Critical spontaneous deactivation rate.
function criticalμ(params::Nondim, N::defInt)::defFloat
    num = params.ρ^2*(1/(N*params.η) + params.τ)
    den = 2*params.β*sin( params.α )
    return num/den
end

println( "Loaded geom.jl class file." )
