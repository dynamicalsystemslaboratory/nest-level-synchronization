
include( "../args.jl" )
include( "../model.jl" )
include( "../geom.jl" )
include( "../recur.jl" )

include( "rho-list.jl" )
include( "beta-list.jl" )

using ProgressMeter

base = "../data/results/"

println( "Initializing simulation parameters: $(ARGS)." )

# Number of replicates.
δt = 1e-2

# Social parameters.
ρ = ρlist[parse( defInt, ARGS[1] )]

# System parameters.
N = 100

# Activity transition variables.
η = 0.0002;  γ = 0.01;  τ = 750.0

# Movement variables.
ξ = 1.0;  λ = 0.50;  ϕ = 0.25
s = 10.0

# Sensing area parameter.
α = (1/4)π

# Get speed value.
μ = 1.0
ℓ = √(N/μ)
println( (ρ, Nβ) )

# Generate parameter variable.
nondimlist = [Nondim(; ρ=ρ, η=η, β=β, γ=γ, τ=τ, ξ=ξ, λ=λ, ϕ=ϕ, s=s, α=α )
    for β ∈ βlist]

# Data folder name.
folderlist = [findfolder( N, μ, nondim; base=base, warn=false )
    for nondim ∈ nondimlist]

# Run simulation under each environment parameter.
T = 100_000

# Compute simulation time-step.
Nt = round( defInt, T/δt );  tlist = 1:Nt
nt = round( defInt, 10/(2*δt) );  tsave = Set( 1:nt:Nt )

# Frequency of adjacency calculation.
δt̂ = round( defInt, 0.1/δt );

# Initialize progress meter.
P = Progress( Nβ*Nt, "Running simulation." )

println( "Beginning simulation." )

# Initialize list to save agent states.
zlist = Vector{State}( undef, Nβ )
φdatalist = Vector{Matrix{defInt}}( undef, Nβ )
@threads for i ∈ 1:Nβ
    # Initialize agent states.
    nondim = nondimlist[i]
    z = initialstate( N, ℓ; A=1 )
    ẑ = copystate( z )

    # Initialize state data.
    φdatalist[i] = Matrix{defInt}( undef, length( tsave ) + 1, N )
    φdatalist[i][1,:] = ẑ.φ

    # Initialize adjacency and saved state.
    A = proximity( N, ℓ, nondim.r, nondim.α, z.x, z.y, z.θ )

    # Run simulation.
    t̂ = 2
    for t ∈ tlist
        # Update the adjacency matrix.
        (t % δt̂) == 0 && (A = proximity( N, ℓ, nondim.r, nondim.α, z.x, z.y, z.θ ))

        # Step simulation.
        step!( N, ℓ, nondim, z, ẑ; A=A, δt=δt )

        # Save state if in appropriate subset.
        t ∈ tsave && (φdatalist[i][t̂,:] = ẑ.φ; t̂ += 1)

        # Swap contents.
        tmp = z;  z = ẑ;  ẑ = tmp

        # Update progress meter.
        next!( P )
    end

    # Save last simulation state.
    zlist[i] = z
end

println( "Saving data." )

# Save data.
for i ∈ 1:Nβ
    # Save final steps.
    if !isdir( folderlist[i] )
        mkpath( folderlist[i] )
    end
    # Save the composite state variable.
    writedlm( folderlist[i]*"state-composition_T-$(T)_M-1.txt", φdatalist[i] )
end
