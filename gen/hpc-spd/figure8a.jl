
include( "../args.jl" )
include( "../model.jl" )
include( "../geom.jl" )
include( "../recur.jl" )

include( "../speed-list-heat.jl" )

using ProgressMeter

loadsteps = length( ARGS ) > 2 ? parse( Bool, ARGS[3] ) : true
savesteps = !loadsteps
loadsteps, savesteps

base = "../data/results/"

println( "Initializing simulation parameters: $(ARGS)." )


# Number of replicates.
M = 50
δt = 1e-3

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
folder = findfolder( N, μ, nondim; base=base, warn=loadsteps )

# Run simulation under each environment parameter.
T = 5000
Tload = 5000  # If applicable.

# Compute simulation time-step.
Nt = round( defInt, T/δt );  tlist = 1:Nt
nt = round( defInt, 10/(2*δt) );  tsave = Set( 1:nt:Nt )

# Frequency of adjacency calculation.
δt̂ = round( defInt, 0.01/δt );

# Initialize progress meter.
P = Progress( M*Nt, "Running simulation."; dt=1.0 )

println( "Beginning simulation." )

# Initialize list and run optimization.
xdata = [Matrix{defFloat}( undef, length( tsave ) + 1, 3 ) for _ ∈ 1:M]
zlist = Vector{State}( undef, M )
@threads for m ∈ 1:M
    # If steps are already saved, use as initial state.
    file = loadsteps ? folder*"steps/state_T-$(Tload)_m-$(m).txt" : nothing

    # Initialize agent states.
    z = initialstate( N, ℓ; A=1, file=file )
    ẑ = copystate( z )

    # Initialize adjacency and saved state.
    A = proximity( N, ℓ, nondim.r, nondim.α, z.x, z.y, z.θ )
    xdata[m][1,:] = statecomposition( N, ẑ )

    # Run simulation.
    t̂ = 2
    for t ∈ tlist
        # Update the adjacency matrix.
        (t % δt̂) == 0 && (A = proximity( N, ℓ, nondim.r, nondim.α, z.x, z.y, z.θ ))

        # Step simulation.
        step!( N, ℓ, nondim, z, ẑ; A=A, δt=δt )

        # Save state if in appropriate subset.
        if t ∈ tsave
            xdata[m][t̂,:] = statecomposition( N, ẑ )
            t̂ += 1
        end

        # Swap contents.
        tmp = z;  z = ẑ;  ẑ = tmp

        # Update progress meter.
        next!( P )
    end

    # Save last simulation state.
    zlist[m] = z
end

println( "Simulation finished, computing determinism." )

# Compute determinism metric and related statistics.
ε0 = 10.0.^(round( log10( 1/N ) ):0.05:-1)
εlist = Vector{defFloat}( undef, M )
Rlist = Vector{RecurrenceMap}( undef, M )
ςlist = Vector{defFloat}( undef, M )
@threads for m ∈ 1:M
    εlist[m], Rlist[m] = calibraterecur( xdata[m]; δxlist=ε0, RR=0.025 )
    ςlist[m] = determinism( Rlist[m]; ℓ0=10 )
end

println( "Determinism computed, saving data." )

# Create folder and save data.
if !isdir( folder*"/steps/" )
    mkpath( folder*"/steps/" )
end

if loadsteps
    writedlm( folder*"activity_T-$(round( defInt, T ))_M-$(M).txt", [xlist[:,1] for xlist ∈ xdata] )
    writedlm( folder*"inactivity_T-$(round( defInt, T ))_M-$(M).txt", [xlist[:,2] for xlist ∈ xdata] )
    writedlm( folder*"refractory_T-$(round( defInt, T ))_M-$(M).txt", [xlist[:,3] for xlist ∈ xdata] )
    writedlm( folder*"determinism_T-$(round( defInt, T ))_M-$(M).txt", ςlist )
end

for m ∈ 1:M
    savestate( folder*"steps/state_T-$(T)_m-$(m).txt", zlist[m] )
end
