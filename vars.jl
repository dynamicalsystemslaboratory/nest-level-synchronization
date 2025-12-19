
include( "args.jl" )

using Distributions
using JSON3, DelimitedFiles

# Global default interaction radius.
const r = defFloat( 1.0 )

# ----------------- Initialization Helpers ----------------- #
function computeμ(μ::defFloat, σ::defFloat)
    return log( μ ) - σ^2/2
end

function initlognorm(; m::Float64=1/10, c::Float64=1/2)
    σ = √log( 1 + c^2 )
    μ  = computeμ( m, σ )
    return LogNormal( μ, σ )
end

function createdistributions(; ξ::defFloat=1/10, c::defFloat=1/2,
    ω::defFloat=0.0, λ::defFloat=1/2)::Tuple{Distribution, Distribution}
    # Event duration.
    Ξ = initlognorm(; m=ξ, c=c )

    # Change in angle.
    Ω = truncated( Laplace( ω, λ ), -π, π )

    return Ξ, Ω
end

# ----------------- Parameter Structure -------------------- #
mutable struct Params
    # Activity transition variables.
    β::defFloat
    γ::defFloat
    ρ::defFloat
    η::defFloat
    τ::defFloat

    # Movement distribution variables.
    ξ::defFloat
    c::defFloat
    ω::defFloat
    λ::defFloat
    ϕ::defFloat
    s::defFloat

    # Intereaction variables.
    r::defFloat
    α::defFloat

    # Distributions for noise parameters.
    Ξ::Distribution
    Ω::Distribution
end

Params(; β::defFloat=1.0, γ::defFloat=1/10,
    ρ::defFloat=0.0, η::defFloat=0.0, τ::defFloat=75.0,
    ξ::defFloat=0.1, c::defFloat=0.5, ω::defFloat=0.0, λ::defFloat=0.5,
    ϕ::defFloat=0.25, s::defFloat=10.0,
    r::defFloat=r, α::defFloat=π/4) = Params( β, γ, ρ, η, τ, ξ, c, ω, λ, ϕ, s, r, α,
        createdistributions(; ξ=ξ, c=c, ω=ω, λ=λ )... )

function adjparams(params::Params; k::defFloat=Inf,
    β::defFloat=Inf, γ::defFloat=Inf,
    ρ::defFloat=Inf, η::defFloat=Inf, τ::defFloat=Inf,
    ξ::defFloat=Inf, c::defFloat=Inf, ω::defFloat=Inf, λ::defFloat=Inf,
    ϕ::defFloat=Inf, s::defFloat=Inf,
    r::defFloat=Inf, α::defFloat=Inf)
    return Params(;
        β=isfinite( β ) ? β : params.β,
        γ=isfinite( γ ) ? γ : params.γ,
        ρ=isfinite( ρ ) ? ρ : params.ρ,
        η=isfinite( η ) ? η : params.η,
        τ=isfinite( τ ) ? τ : params.τ,
        ξ=isfinite( ξ ) ? ξ : params.ξ,
        c=isfinite( c ) ? c : params.c,
        ω=isfinite( ω ) ? ω : params.ω,
        λ=isfinite( λ ) ? λ : params.λ,
        ϕ=isfinite( ϕ ) ? ϕ : params.ϕ,
        s=isfinite( s ) ? s : params.s,
        r=isfinite( r ) ? r : params.r,
        α=isfinite( α ) ? α : params.α,
    )
end

function saveparams(filename::String, params::Params)
    tmp = Dict(
        :ρ => params.ρ, :η => params.η,
        :β => params.β, :γ => params.γ, :τ => params.τ,
        :ξ => params.ξ, :c => params.c, :ω => params.ω,
        :λ => params.λ, :ϕ => params.ϕ, :s => params.s,
        :r => params.r, :α => params.α
    )
    open( filename, "w" ) do io
        JSON3.write( io, tmp )
    end
end

function loadparams(filename::String)::Params
    tmp = JSON3.read(open(filename, "r"), Dict{String,defFloat})

    return Params(;
        β=tmp["β"], γ=tmp["γ"], ρ=tmp["ρ"], η=tmp["η"], τ=tmp["τ"],
        ξ=tmp["ξ"], c=tmp["c"], ω=tmp["ω"], λ=tmp["λ"], ϕ=tmp["ϕ"], s=tmp["s"],
        r=tmp["r"], α=tmp["α"]
    )
end

# ------------ Non-dimensional Scale Parameters ------------ #
struct Scale
    # Characteristic length and time scales.
    L::defFloat
    T::defFloat
end

function savescale(filename::String, scale::Scale)
    tmp = Dict( :L => scale.L, :T => scale.T )
    open( filename, "w" ) do io
        JSON3.write( io, tmp )
    end
end

function loadscale(filename::String)::Scale
    tmp = JSON3.read(open(filename, "r"), Dict{String,defFloat})
    return Scale( tmp["L"], tmp["T"] )
end

# ------------ Non-dimensional Model Parameters ------------ #
struct Nondim
    # Activity transition variables.
    β::defFloat
    γ::defFloat
    ρ::defFloat
    η::defFloat
    τ::defFloat

    # Movement distribution variables.
    ξ::defFloat
    c::defFloat
    ω::defFloat
    λ::defFloat
    ϕ::defFloat
    s::defFloat

    # Intereaction variables.
    r::defFloat
    α::defFloat

    # Distributions for noise parameters.
    Ξ::Distribution
    Ω::Distribution
end

function Nondim(params::Params; scale::Scale=Scale( 1.0, 0.1 ))::Nondim
    return Nondim(
        params.β*scale.T, params.γ*scale.T, params.ρ*scale.T,
        params.η*scale.T, params.τ/scale.T, params.ξ/scale.T,
        params.c, params.ω, params.λ, params.ϕ,
        params.s*scale.T^(params.ϕ + 1)/scale.L,
        params.r/scale.L, params.α,
        createdistributions(; ξ=params.ξ/scale.T, c=params.c, ω=params.ω, λ=params.λ )... )
end

# ---------------- State Variable Structure ---------------- #
mutable struct State
    x::Vector{defFloat}
    y::Vector{defFloat}
    θ::Vector{defFloat}
    φ::Vector{defInt}
    τ::Vector{defInt}
    ψ::Vector{defInt}
    v::Vector{defFloat}
end

function State(N::defInt, L::defFloat; A::defInt=1, B::defInt=N-A, iτ::defInt=-1)::State
    random = A < 0
    A = random ? rand( defInt.( 1:N ) ) : A
    B = random ? (A < N ? rand( defInt.( 1:N-A ) ) : 0) : B
    R = N - A - B
    return State(
        L*rand( defFloat, N ),
        L*rand( defFloat, N ),
        2*π*rand( defFloat, N ) .- π,
        [zeros( defInt, A ); ones( defInt, B ); fill( defInt( 2 ), R )],
        iτ < 0 ? zeros( defInt, N ) : rand( defInt.( 1:iτ ), N ),
        zeros( defInt, N ),
        zeros( defFloat, N ) )
end

function State(state::Matrix{Float64})
    return State(
        state[1,:],
        state[2,:],
        state[3,:],
        round.( defInt, state[4,:] ),
        round.( defInt, state[5,:] ),
        round.( defInt, state[6,:] ),
        state[7,:] )
end

function copystate(z::State)::State
    return State(
        copy( z.x ),
        copy( z.y ),
        copy( z.θ ),
        copy( z.φ ),
        copy( z.τ ),
        copy( z.ψ ),
        copy( z.v ) )
end

function State(N::defInt, L::defFloat,
    x::defFloat=L/2, y::defFloat=L/2, θ::defFloat=π/2, φ::defInt=0)::State
    # Create initial conditions.
    return State(
        [x; L*rand( defFloat, N-1 )],
        [y; L*rand( defFloat, N-1 )],
        [θ; 2*π*rand( defFloat, N-1 ) .- π],
        [φ; ones( defInt, N-1 )],
        zeros( defInt, N ),
        zeros( defInt, N ),
        zeros( defFloat, N ) )
end

# Save given state.
function savestate(filename::String, z::State)
    # Unpack state variable.
    rawstate = [z.x, z.y, z.θ, z.φ, z.τ, z.ψ, z.v]
    return writedlm( filename, rawstate )
end

# Load state from file.
function loadstate(filename::String)::State
    # Unpack state variable.
    rawstate = readdlm( filename )
    return State( rawstate )
end

# Initial state function - capabale of importing from existing state file.
function initialstate(N::defInt, L::defFloat; A::defInt=1, file::Union{Nothing,String}=nothing)
    if isnothing(file)
        return State( N, L; A=A )
    else
        return loadstate( file )
    end
end
