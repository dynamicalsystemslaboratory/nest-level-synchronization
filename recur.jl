
include( "args.jl" )

using SparseArrays
using LinearAlgebra

mutable struct RecurrenceMap
    t::Vector{defInt}
    R::SparseMatrixCSC{Bool,defInt}
end

function Base.length(R::RecurrenceMap)
    return length( R.t )
end

function Base.sum(R::RecurrenceMap)
    return sum( R.R )
end

function proximity(x::Vector{defFloat}, x̂::Vector{defFloat})::defFloat
    return norm( x - x̂ )
end

function cosine(x::Vector{defFloat}, x̂::Vector{defFloat}; ε::defFloat=1e-12)::defFloat
    nx = norm( x )
    nx̂ = norm( x̂ )

    if nx < ε || nx̂ < ε
        return Inf
    end

    return 1/2*(1 - dot( x, x̂ )/(norm( x )*norm( x̂ )))
end

function createbuckets(xlist::Matrix{defFloat}, δx::defFloat)
    bdict = Dict{defInt, Vector{defInt}}()
    for t ∈ 1:size( xlist )[1]
        b = floor( defInt, norm( xlist[t,:] )/δx )
        push!( get!( bdict, b, Vector{defInt}() ), t )
    end
    return bdict
end

function recurrence(xlist::AbstractArray{defFloat};
    δx::defFloat=1e-2, step::defInt=1, distance::Function=proximity)::RecurrenceMap
    T = size( xlist )[1]
    tlist = 1:step:T;  Nt = length( tlist )

    # Make sparse index list.
    rlist, clist = Vector{defInt}(), Vector{defInt}()

    # Initialize buckets for efficiency.
    xsub = xlist[tlist,:]
    buckets = createbuckets( xsub, δx )

    # Compare only within same or neighboring buckets
    for (b, indices) in buckets
        nlist = vcat(
            get( buckets, b-1, Vector{defInt}() ),
            get( buckets, b, Vector{defInt}() ),
            get( buckets, b+1, Vector{defInt}() ) )

        for i in indices, j in nlist
            i < j && distance( xsub[i,:], xsub[j,:] )^2 < δx^2 || continue

            push!( rlist, i );  push!( clist, j )
            push!( rlist, j );  push!( clist, i )
        end
    end

    return RecurrenceMap( tlist, sparse( rlist, clist, ones( Bool, length( rlist ) ), Nt, Nt ) )
end

function calibraterecur(xlist::AbstractArray{defFloat};
    RR::defFloat=0.025, δxlist::AbstractArray{defFloat}=1.0./(10:10:100), args...
    )::Tuple{defFloat, RecurrenceMap}
    # For each threshold compile the recurrence map.
    Rlist = [recurrence( xlist; δx=δx, args... ) for δx ∈ δxlist]

    # Compute the recurrence rate of each map.
    RRlist = recurrate.( Rlist )

    # Return the closest threshold value and the corresponding map.
    i = argmin( abs.( RR .- RRlist ) )
    return δxlist[i], Rlist[i]
end

function recurrate(R::RecurrenceMap)::defFloat
    return sum( R.R )/length( R.t )^2
end

function diagk(R::RecurrenceMap, k::Int)
    N = size( R.R,1 )
    return k < 0 ? [R.R[i-k,i] for i in 1:(N+k)] : [R.R[i,i+k] for i in 1:(N-k)]
end

function countsegments(rlist::Vector{Bool})
    ℓlist = Vector{defInt}()

    ℓ = 0
    for r ∈ rlist
        if r
            ℓ += 1
        else
            ℓ == 0 || push!( ℓlist, ℓ )
            ℓ = 0  # Reset current length.
        end
    end

    # Catch the last segment.
    ℓ > 0 && push!( ℓlist, ℓ )

    return ℓlist
end

function extractlinelengths(R::RecurrenceMap, ℓ0::defInt)::Dict{defInt,defInt}
    # Initialize dictionary.
    N = length( R )
    ℓd = Dict{defInt,defInt}()

    # Count line segments.
    for k in -(N - ℓ0):(N - ℓ0)
        # Extract the k-th diagonal.
        rlist = diagk( R, k )
        sum( rlist ) > 0 || continue

        # Count the number of contiguous segments.
        ℓlist = countsegments( rlist )
        length( ℓlist ) > 0 || continue

        # Add the segment lengths to the dictionary.
        for ℓ ∈ ℓlist
            ℓd[ℓ] = get( ℓd, ℓ, 0 ) + 1
        end
    end

    # Return the compiled dictionary.
    return ℓd
end

function determinism(R::RecurrenceMap; ℓd::Union{Nothing,Dict{defInt,defInt}}=nothing,
    ℓ0::defInt=10, verbose::Bool=false)
    # Count line lengths if necessary.
    ℓd = isnothing( ℓd ) ? extractlinelengths( R, ℓ0 ) : ℓd

    # Compute the determinism and return: DET = ∑ℓ*N(ℓ)/∑R.
    ς = sum( [ℓ*n for (ℓ, n) in ℓd if ℓ0 ≤ ℓ] )/sum( R )
    return verbose ? (ς, ℓd) : ς
end

function calibratedeter(R::RecurrenceMap; ℓ0::defInt=2,
    ℓd::Union{Nothing,Dict{defInt,defInt}}=nothing, verbose::Bool=false)
    # Unpack the default determinism and line lengths.
    ℓd = isnothing( ℓd ) ? determinism( R; ℓ0=2, verbose=true )[2] : ℓd

    # Extract the number of minima to check.
    ℓmax = maximum( keys( ℓd ) )

    # Create line length distribution.
    p = Vector{defFloat}( [haskey( ℓd, ℓ ) ? ℓd[ℓ] : 0 for ℓ ∈ 1:ℓmax] )
    p ./= sum( p )

    # Compute Otsu weights.
    ℓlist = ℓ0:ℓmax-1
    wblist = [sum( p[1:ℓmin-1] ) for ℓmin ∈ ℓlist]
    wflist = [sum( p[ℓmin:ℓmax] ) for ℓmin ∈ ℓlist]

    # Compute expected value of partitions.
    μblist = [sum( (1:ℓmin-1).*p[1:ℓmin-1] ) for ℓmin ∈ ℓlist]./wblist
    μflist = [sum( (ℓmin:ℓmax).*p[ℓmin:ℓmax] ) for ℓmin ∈ ℓlist]./wflist

    # Compute the between-class variance.
    σlist = wblist.*wflist.*(μblist .- μflist).^2

    # Return the optimal line length and corresponding determinism.
    i = argmax( σlist )
    ℓ = ℓlist[i]
    ς = determinism( R; ℓd=ℓd, ℓ0=ℓ )
    return verbose ? (ℓ, ς, ℓd) : (ℓ, ς)
end

function doubledeter(R::RecurrenceMap; verbose::Bool=false)
    # Extract the threshold using first Otsu.
    (ℓ1, ς1, ℓd1) = calibratedeter( R; verbose=true )

    # Unpack the line lengths that are above the threshold.
    ℓmax = maximum( keys( ℓd1 ) )
    ℓlist = ℓ1:ℓmax
    ℓd2 = Dict{defInt,defInt}( [ℓ => (haskey( ℓd1, ℓ ) ? ℓd1[ℓ] : 0) for ℓ ∈ ℓlist] )

    # Recompute the threshold.
    ℓ2, ς2 = calibratedeter( R; ℓ0=ℓ1, ℓd=ℓd2 )
    return verbose ? ((ℓ2, ς2), (ℓ1, ς1), ℓd1) : (ℓ2, ς2)
end

function correlate(R::RecurrenceMap, τ::defInt)
    N = length( R )
    return sum( [R.R[i,i+τ] for i ∈ 1:N-τ] )/(N - τ)
end

function correlate(R::RecurrenceMap; τlist::Vector{defInt}=R.t)
    return [correlate( R, τ ) for τ ∈ τlist]
end

println( "Loaded recur.jl class file." )
