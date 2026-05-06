
using Base.Threads
println( "Number of threads: ", nthreads() )

include( "args.jl" )
include( "vars.jl" )
include( "network.jl" )

# --------------------- Spatial Model Functions --------------------- #
function periodicturn(params::Nondim, N::defInt, θ::Union{defFloat,Vector{defFloat}}
    )::Union{defFloat,Vector{defFloat}}
    return N == 1 ? θ + rand( params.Ω ) : θ .+ rand( params.Ω, N )
end

function periodicstep(N::defInt, L::defFloat,
    x::Union{defFloat,Vector{defFloat}},
    y::Union{defFloat,Vector{defFloat}},
    θ::Union{defFloat,Vector{defFloat}},
    v::Union{defFloat,Vector{defFloat}};
    δt::defFloat=Δt
    )::Tuple{Union{defFloat,Vector{defFloat}},Union{defFloat,Vector{defFloat}}}

    # Compute new positions.
    x̂ = x .+ δt.*v.*cos.( θ )
    ŷ = y .+ δt.*v.*sin.( θ )

    # Return wrapped positions.
    return (mod.( x̂, L ), mod.( ŷ, L ))
end

function drawduration(N::defInt, params::Nondim)::Union{defFloat,Vector{defFloat}}
    return N == 1 ? rand( params.Ξ ) : rand( params.Ξ, N )
end

function agnosticspeed(a::Vector{defInt},
    ψ::Union{defFloat,Vector{defFloat}},
    params::Nondim)::Union{defFloat,Vector{defFloat}}
    return (params.s).*ψ.^(params.ϕ)
end

# ------------------- Dynamics Function Structure ------------------- #
mutable struct Dynamics
    speed::Function
    step::Function
    distance::Function
end

function Dynamics(; speed::Function=agnosticspeed,
    step::Function=periodicstep, distance::Function=radialdistance)::Dynamics
    return Dynamics( speed, step, distance )
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
    φ == 2 && (return Vector{defInt}( [1] ))  # If 2 (refractory) return 1.
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
function torus(L::defFloat, z::State)::Bool
    # Check that all agents are inside the torus.
    return all( 0 .< z.x .< L ) && all( 0 .< z.y .< L )
end

function inactive!(i::defInt, z::State, ẑ::State)::Bool
    # Check if the state is inactive.
    !((ẑ.φ[i] == 1) || (ẑ.φ[i] == 2)) && (return false)

    # Zero velocity and bout duration.
    ẑ.ψ[i] = 0.0;  ẑ.v[i] = 0.0

    # Immobile position.
    ẑ.x[i], ẑ.y[i], ẑ.θ[i] = z.x[i], z.y[i], z.θ[i]

    # Return update boolean.
    return true
end

function active!(i::defInt, L::defFloat, params::Nondim,
    a::Vector{defInt}, z::State, ẑ::State;
    δt::defFloat=Δt, dvar::Dynamics=Dynamics())
    # Check if active.
    !(ẑ.φ[i] == 0) && (return false)

    # Decrement bout duration by 1.
    ẑ.ψ[i] = z.ψ[i] - 1

    # Check if a new bout duration is needed.
    if 0 < ẑ.ψ[i]
        # Constant heading and velocity.
        ẑ.θ[i] = z.θ[i];  ẑ.v[i] = z.v[i]
    else
        # Draw new heading.
        ẑ.θ[i] = periodicturn( params, 1, z.θ[i] )

        # Draw new bout duration.
        ψ = drawduration( 1, params )
        ẑ.ψ[i] = round( defInt, ψ./δt )

        # Compute new speed.
        ẑ.v[i] = dvar.speed( a, ψ, params )
    end

    # Update agent position.
    ẑ.x[i], ẑ.y[i] = dvar.step( 1, L, z.x[i], z.y[i], ẑ.θ[i], ẑ.v[i]; δt=δt )

    # Return update boolean.
    return true
end

function safety(N::defInt, L::defFloat, z::State; bound::Function=torus)::Bool
    # Check that all agents have an activity state.
    ia = z.φ .== 0
    ib = (z.φ .== 1) .| (z.φ .== 2)
    @assert sum( ia ) + sum( ib ) == N "Not all agents are being tracked."

    # Check that active agents have a movement duration.
    @assert sum( ia .& (z.ψ .> 0) ) + sum( ia .& (z.ψ .≤ 0) ) == sum( ia ) "Not all active agents are being tracked."

    # Check legitimacy of movement duration.
    @assert sum( z.ψ[ia] .< 0 ) == 0 "At least one active duration length is negative."

    # Check boundary constraints.
    inbounds = bound( L, z )
    if !inbounds
        println( ia, " ", ib, " " )
        println( z.θ )
        println( z.x )
        println( z.y )
    end
    @assert inbounds "Agents have escaped the boundary."
    return true
end

function step!(N::defInt, L::defFloat, params::Nondim, z::State, ẑ::State;
    δt::defFloat=Δt, A::Union{Nothing,SparseMatrixCSC{defInt,defInt}}=nothing,
    dvar::Dynamics=Dynamics(), bound::Function=torus, safe::Bool=false
    )::SparseMatrixCSC{defInt,defInt}
    # Compute adjacency matrix unless given.
    A = isnothing( A ) ? proximity( N, L, params.r, params.α, z.x, z.y, z.θ; distance=dvar.distance ) : A

    # Transition agents between active, inactive, and refractive states.
    @threads for i ∈ 1:N
        @inbounds begin

        # Transition between activity states.
        transition = switch( A[i,:].nzind, z.φ[i], z.φ, z.τ[i], params; δt=δt )
        ẑ.φ[i] = transition ? nextstate( z.φ[i] ) : z.φ[i]
        ẑ.τ[i] = transition ? 0 : z.τ[i] + 1

        # Update the movement states.
        check = inactive!( i, z, ẑ ) || active!( i, L, params, A[i,:].nzind, z, ẑ; δt=δt, dvar=dvar )

        # If neither active states reached, trigger assertion.
        @assert check "Agent is neither active nor inactive (φ = $(ẑ.φ[i]))."

        end
    end

    # Check that model is well-behaved.
    safe && safety( N, L, ẑ; bound=bound )
    return A
end

function simulate(N::defInt, L::defFloat, α::defFloat, params::Nondim, T::Real, z0::State;
    δt::defFloat=Δt, δt̂::Union{Nothing,defFloat}=nothing, distance::Function=radialdistance)
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
        step!( N, L, params, z, ẑ; δt=δt, dvar=Dynamics(; distance=radialdistance ) )

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
