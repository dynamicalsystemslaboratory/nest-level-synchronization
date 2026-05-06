
include( "vars.jl" )

using StaticArrays

# Compute the phase coherence.
function computephase(Δφ::Vector{defInt})::Vector{defFloat}
    # Extract the transition points.
    Nt = length( Δφ ) + 1
    ι = [1;findall( Δφ .== 1 );Nt]

    # Iterate through each index, setting phase angle.
    θ = vcat( [2*π*(0.0:(1.0/((ι2-1) - ι1)):1.0) for (ι1, ι2) ∈ zip( ι[1:end-1], ι[2:end] )]... )

    # Return the phase angle list.
    @assert length( θ ) == Nt-1 "Phase angle list is not the correct length: $(length( θ ))."
    return θ
end

function euler(θ::defFloat)::Complex{defFloat}
    return cos( θ ) + im*sin( θ )
end

function order(θ::Vector{defFloat})
    # Return the coherence.
    return abs( sum( euler.( θ ) )/length( θ ) )
end

# Pair-wise distance function.
function pairwisedist(L::defFloat, xlist::Vector{defFloat}, ylist::Vector{defFloat}; n::defInt=-1)::Vector{defFloat}
    # Initialize distance matrix variables.
    N = length( xlist )

    # Consider subset of full distance set.
    N̂ = n < 0 ? N : n
    ιlist = sort( sample( 1:N, N̂; replace=false ) )
    D = Vector{defFloat}( undef, round( defInt, N̂*(N̂-1)/2 ) )

    # Iterate through, calculating distance.
    ι = 1
    for (k, i) ∈ enumerate( ιlist )
        for j ∈ ιlist[k+1:end]
            # Horizontal and vertical difference (with periodic boundary).
            dx = xlist[i] - xlist[j];  dx -= round( dx/L ) * L
            dy = ylist[i] - ylist[j];  dy -= round( dy/L ) * L

            # Distance.
            D[ι] = √(dx*dx + dy*dy);  ι += 1
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

# Compute estimate of the stationary area.
function stationaryarea(params::Params)::defFloat
    # Proportion of static area unique.
    # TODO: Scale by the probability that the agent leaves its current interaction radius.
    ξ0 = (params.r/params.s)^(1/params.ϕ)
    g = isnothing( params.Ξ ) ? 1.0 : (cdf( params.Ξ, Inf ) - cdf( params.Ξ, ξ0 ))
    # g = isnothing( λ ) ? 1.0 : (exp( -λ*α ) - exp( -λ*π ))/(1 - exp( -λ*π ))

    # Offset term.
    A0 = params.α*params.r^2/(params.ρ*params.ξ)
    return g*A0
end

# Compute estimate of the dynamic area.
function dynamicarea(params::Params)::defFloat
    # Proportion of dynamic area unique.
    f = isnothing( params.λ ) ? 1.0 : 1 - (exp( -params.λ*π )*(exp( params.λ*α ) - 1)/(1 - exp( -params.λ*π )))

    # Scaling term.
    A = 2*params.r*params.s*(params.ξ^params.ϕ)/params.ρ*(params.α < π ? sin( params.α ) : 1.0)
    return f*A
end

# Analytically estimated area coverage.
function areaestimate(params::Params)::defFloat
    # Return expected total area coverage.
    return dynamicarea( params ) + stationaryarea( params )
end

# Critical transition in area coverage.
function areacriticalξ(params::Params)::defFloat
    return (params.α*params.r/(2*params.s*params.ϕ*sin( params.α )))^(1/(params.ϕ + 1))
end

function criticalξ(params::Params, μ::defFloat)::defFloat
    # if params.ρ/params.β
    return (params.ρ/(2*params.β*params.γ*params.r*μ*params.s*sin( params.α )))^(1/params.ϕ)
end

# Activity duty cycle.
function dutycycle(params::Params, N::defInt)::defFloat
    return 1/(1 + params.ρ*(1/(N*params.η) + params.τ))
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
