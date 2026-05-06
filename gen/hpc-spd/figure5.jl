
include( "../args.jl" )
include( "../model.jl" )
include( "../geom.jl" )
include( "../recur.jl" )

include( "../speed-list-heat.jl" )

using ProgressMeter

base = "../data/results/"

println( "Initializing simulation parameters: $(ARGS)." )


# Number of replicates.
δt = 1e-2

# Social parameters.
ρ = 0.01;  β = 0.1

# System parameters.
μ = 1.0
N = parse( defInt, ARGS[1] )
ℓ = √(N/μ)

# Activity transition variables.
η = 0.0002;  γ = 0.01;  τ = 750.0

# Movement variables.
ξ = 1.0;  λ = 0.50;  ϕ = 0.25

# Sensing area parameter.
α = (1/4)π

# Get speed value.
s = slist[parse( defInt, ARGS[2] )]

# Generate parameter variable.
nondim = Nondim(; ρ=ρ, η=η, β=β, γ=γ, τ=τ, ξ=ξ, λ=λ, ϕ=ϕ, s=s, α=α )
println( "nondim. parameters: ", nondim )

# Data folder name.
folder = findfolder( N, μ, nondim; base=base )

# Run simulation under each environment parameter.
T = 100_000
Tload = 5000  # If applicable.

# Compute simulation time-step.
Nt = round( defInt, T/δt );  tlist = 1:Nt
nt = round( defInt, 10/(2*δt) );  tsave = Set( 1:nt:Nt )

# Frequency of adjacency calculation.
δt̂ = round( defInt, 0.01/δt );

# Initialize progress meter.
P = Progress( Nt, "Running simulation." )

println( "Beginning simulation." )

# Initialize list to save agent states.
φdata = Matrix{defInt}( undef, length( tsave ) + 1, N )

# Initialize agent states.
z = initialstate( N, ℓ; A=1, file=folder*"steps/state_T-$(Tload)_m-1.txt" )
ẑ = copystate( z )

# Initialize adjacency and saved state.
A = proximity( N, ℓ, nondim.r, nondim.α, z.x, z.y, z.θ )

# Run simulation.
t̂ = 2
for t ∈ tlist
    # Update the adjacency matrix.
    (t % δt̂) == 0 && (global A = proximity( N, ℓ, nondim.r, nondim.α, z.x, z.y, z.θ ))

    # Step simulation.
    step!( N, ℓ, nondim, z, ẑ; A=A, δt=δt )

    # Save state if in appropriate subset.
    t ∈ tsave && (φdata[t̂,:] = ẑ.φ; global t̂ += 1)

    # Swap contents.
    tmp = z;  global z = ẑ;  global ẑ = tmp

    # Update progress meter.
    next!( P )
end

println( "Saving data." )

# Save data.
writedlm( folder*"state-composition_T-$(T)_M-1.txt", φdata )
