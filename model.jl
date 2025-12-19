
include( "args.jl" )
include( "vars.jl" )
include( "network.jl" )

# --------------------- Spatial Model Functions --------------------- #
function periodicturn(params::Nondim, N::defInt, θlist::Vector{defFloat})::Vector{defFloat}
    return θlist .+ rand( params.Ω, N )
end

function periodicstep(N::defInt, L::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θlist::Vector{defFloat},
    vlist::Vector{defFloat}; δt::defFloat=Δt
    )::Tuple{Vector{defFloat},Vector{defFloat}}

    # Compute new positions.
    xnew = xlist .+ δt.*vlist.*cos.( θlist )
    ynew = ylist .+ δt.*vlist.*sin.( θlist )

    # Return wrapped positions.
    return mod.( xnew, L ), mod.( ynew, L )
end

function drawduration(params::Nondim, N::defInt)::Vector{defFloat}
    return rand( params.Ξ, N )
end

function computespeed(params::Nondim, ψlist::Vector{defFloat})::Vector{defFloat}
    return (params.s).*ψlist.^(params.ϕ)
end

# --------------------- Activity Model Functions -------------------- #
function nextstate(φ::defInt)::defInt  # Transition order: [0,1,2] = [active,inactive,refractory].
    return mod( φ + 2, 3 )
end

function extractgain(φ::defInt, τ::defInt, params::Nondim;
    δt::defFloat=Δt)::Vector{defFloat}
    if φ == 0
        # If active, return γ and ρ scaled by the time-step.
        return ([params.γ, params.ρ].*δt)
    elseif φ == 1
        # If inactive, return β and η scaled by the time-step.
        return ([params.β, params.η].*δt)
    elseif φ == 2
        # If refractory, check if exceed delay condition.
        return [(τ == round( defInt, params.τ/δt ))*1.0]
    end
    @assert false "Current activity state is φ = $(φ) but state should be constrained φ ∈ {0, 1, 2}."
end

function extractcount(a::Vector{defInt}, φ::defInt, φlist::Vector{defInt})::Vector{defInt}
    φ == 2 && (return [1]::Vector{defInt})  # If 2 (refractory) return 1.
    # Otherwise, count number of active connections and return.
    return [sum( φlist[a] .== 0 ), 1]
end

function prob(μlist::Vector{defInt}, ωlist::Vector{defFloat})::Vector{defFloat}
    return ((1 .- ωlist).^μlist)
end

function switch(a::Vector{defInt}, φ::defInt, φlist::Vector{defInt},
    τ::defInt, params::Nondim; δt::defFloat=Δt)::Bool
    μlist = extractcount( a, φ, φlist )
    ωlist = extractgain( φ, τ, params; δt=δt )
    return (rand() < (1 - prod( prob( μlist, ωlist ) )))
end

# ------------------ Agent-based Model Simulation ------------------- #
function step!(N::defInt, L::defFloat, params::Nondim, z::State, ẑ::State;
    δt::defFloat=Δt, distance=radialdistance,
    A::Union{Nothing,SparseMatrixCSC{defInt,defInt}}=nothing)
    # Compute adjacency matrix unless given.
    A = (A === nothing) ? proximity( N, L, params.r, params.α, z.x, z.y, z.θ; distance=distance ) : A

    # Transition agents between active, inactive, and refractive states.
    @inbounds @simd for i ∈ 1:N
        check = switch( A[i,:].nzind, z.φ[i], z.φ, z.τ[i], params; δt=δt )
        ẑ.φ[i] = check ? nextstate( z.φ[i] ) : z.φ[i]
        ẑ.τ[i] = check ? 0 : z.τ[i] + 1
    end

    # Extract indeces of active agents.
    ia = (ẑ.φ .== 0)
    ib = (ẑ.φ .== 1) .| (ẑ.φ .== 2)
    @assert sum( ia ) + sum( ib ) == N "Not all agents are being tracked."

    # Update movement durations appropriately.
    ẑ.ψ[ia] = z.ψ[ia] .- 1
    ẑ.ψ[ib] .= 0.0;  ẑ.v[ib] .= 0.0

    # If inactive/refractory, remain fixed in place.
    ẑ.x[ib], ẑ.y[ib], ẑ.θ[ib] = z.x[ib], z.y[ib], z.θ[ib]

    # From the active population, isolate those that have non-zero movement duration.
    im = ia .& (z.ψ .> 0)
    iz = ia .& (z.ψ .≤ 0)
    @assert sum( im ) + sum( iz ) == sum( ia ) "Not all active agents are being tracked."

    # If active and in a duration, use previous origientation and speed.
    ẑ.θ[im] = z.θ[im];  ẑ.v[im] = z.v[im]

    # If active but not in a duration, update orientation.
    ẑ.θ[iz] = periodicturn( params, sum( iz ), z.θ[iz] )

    # If active but not in a duration, draw new duration of movement.
    ψ = drawduration( params, sum( iz ) )
    @assert sum( ψ .< 0 ) == 0 "At least one duration length is negative."
    ẑ.ψ[iz] = round.( defInt, ψ./δt );  ẑ.v[iz] = computespeed( params, ψ )

    # If active and in a duration, move using velocity.
    ẑ.x[ia], ẑ.y[ia] = periodicstep( sum( ia ), L, z.x[ia], z.y[ia], ẑ.θ[ia], ẑ.v[ia]; δt=δt )

    # Check that positions are bounded.
    inbounds = all( 0 .< ẑ.x .< L ) && all( 0 .< ẑ.y .< L )
    if !inbounds
        println( t, " ", ia, " ", ib, " " )
        println( z.θ, ẑ.θ )
        println( z.x, ẑ.x )
        println( z.y, ẑ.y )
    end
    @assert inbounds "Agents have escaped the boundary."

    return A
end

function simulate(N::defInt, L::defFloat, α::defFloat, params::Nondim, T::Real, z0::State;
    δt::defFloat=Δt, δt̂::Union{Nothing,defFloat}=nothing, distance=radialdistance)
    # Temporal variable(s).
    Nt = round( defInt, T/δt + 1 )
    Nt̂ = round( defInt, 1/(isnothing( δt̂ ) ? 1 : δt̂) )
    tlist = Set( 1:Nt̂:Nt )

    # Initialize save state variables.
    zlist = Vector{State}( undef, length( tlist ) + 1 )
    zlist[1] = copystate( z0 )

    # Initialize leap-frog variables.
    z = copystate( z0 );  ẑ = copystate( z0 )

    # Run simulation using step function.
    t̂ = 2
    for t ∈ 1:Nt
        # Iterate state variable.
        step!( N, L, params, z, ẑ; δt=δt, distance=distance )

        # Save data at increments of δt̂.
        (t ∈ tlist) && (zlist[t̂] = copystate( ẑ );  t̂ += 1)

        # Swap contents.
        tmp = z;  z = ẑ;  ẑ = tmp
    end

    return defFloat.( sort( collect( tlist ) ) ), zlist
end

# Helper for saving data.
function statecomposition(N::defInt, z::State)::Vector{defFloat}
    return [sum( z.φ .== φ )/N for φ ∈ 0:2]
end

# Helper for computing who transitioned in post.
function transitioned(z::State, ẑ::State; φ::defInt=0)
    return (z.φ .!= φ) .& (ẑ.φ .== φ)
end


println( "Loaded model.jl class file." )
