
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

function recurrence(xlist::Matrix{defFloat}; δx::defFloat=1e-3, step::defInt=1, distance=proximity)::RecurrenceMap
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
    ℓ == 0 || push!(ℓlist, ℓ)

    return ℓlist
end

function determinism(R::RecurrenceMap; ℓ0::defInt=10)
    N = length( R )
    ℓd = Dict{defInt,defInt}()

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

    # Compute the determinism and return: DET = ∑ℓ*N(ℓ)/∑R.
    return sum( [ℓ*n for (ℓ, n) in ℓd if ℓ0 ≤ ℓ] )/sum( R )
end

function correlate(R::RecurrenceMap, τ::defInt)
    N = length( R )
    return sum( [R.R[i,i+τ] for i ∈ 1:N-τ] )/(N - τ)
end

function correlate(R::RecurrenceMap; τlist::Vector{defInt}=R.t)
    return [correlate( R, τ ) for τ ∈ τlist]
end

println( "Loaded recur.jl class file." )
