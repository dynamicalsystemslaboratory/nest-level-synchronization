
include( "args.jl" )

using SparseArrays

function anglediff(θ1::defFloat, θ2::defFloat)
    return atan( sin( θ1 - θ2 ), cos( θ1 - θ2 ) )
end

function rawdistance(L::defFloat, r::defFloat, α::defFloat,
    x::Vector{defFloat}, y::Vector{defFloat}, θ::Vector{defFloat})::Bool
    # Compute raw Euclidian distance.
    dx = diff( x )[1]
    dy = diff( y )[1]

    # Return Euclidian distance check.
    return dx^2 + dy^2 < r^2
end

function periodicdistance(L::defFloat, r::defFloat, α::defFloat,
    x::Vector{defFloat}, y::Vector{defFloat}, θ::Vector{defFloat})::Bool
    # Compute raw Euclidian distance.
    dx = abs( diff( x )[1] )
    dy = abs( diff( y )[1] )

    # Adjust for periodic boundary.
    dx = dx > L/2 ? L - dx : dx
    dy = dy > L/2 ? L - dy : dy

    # Return distance values.
    return dx^2 + dy^2 < r^2
end

function radialdistance(L::defFloat, r::defFloat, α::defFloat,
    x::Vector{defFloat}, y::Vector{defFloat}, θ::Vector{defFloat})::Bool
    # Compute raw Euclidian distance.
    dx = abs( diff( x )[1] )
    dy = abs( diff( y )[1] )

    # Adjust for periodic boundary.
    dx = dx > L/2 ? L - dx : dx
    dy = dy > L/2 ? L - dy : dy

    # First check if agents are within periodic interaction radius.
    dx^2 + dy^2 < r^2 || return false

    # Compute the angle between the agents.
    δα = atan( dy, dx )

    # Compute angle difference between each agent's heading and δα.
    dθ1 = abs( anglediff( θ[1], δα ) )
    dθ2 = abs( anglediff( θ[2], δα + π ) )

    # Return true if either are facing each other within threshold α.
    return (dθ1 < α) || (dθ2 < α)
end

function createbuckets(N::defInt, r::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat})::Dict{Tuple{defInt,defInt}, Vector{defInt}}
    # Grid index for each agent.
    ilist = floor.( defInt, xlist./r .+ 1)
    jlist = floor.( defInt, ylist./r .+ 1)

    # Create dictionary of bucket-index combinations.
    bdict = Dict{Tuple{defInt,defInt}, Vector{defInt}}()
    for i in 1:N
        # Pack the agent index into a single variable.
        key = (ilist[i], jlist[i])

        # For each agent, save its index to the appropriate grid bucket.
        push!( get!( bdict, key, Vector{defInt}() ), i )
    end

    return bdict
end

function proximity(N::defInt, L::defFloat, r::defFloat, α::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θlist::Vector{defFloat};
    distance=radialdistance)::SparseMatrixCSC{defInt,defInt}
    # Partition agent list into grid-based buckets.
    bdict = createbuckets( N, r, xlist, ylist )

    # Initialize list of column/row indeces for connection tracking.
    rlist, clist = Vector{defInt}(), Vector{defInt}()

    # For each unique grid ID contained in the buckets dictionary.
    for (key, ilist) in bdict
        # Unpack current grid bucket.
        i, j = key

        # For each grid bucket adjacent to the current bucket.
        for di in -1:1, dj in -1:1
            # Get current grid bucket to check.
            nkey = (i + di, j + dj)

            # If grid box is empty, skip loop iteration.
            haskey( bdict, nkey ) || continue

            # Unpack neighbors in adjacent grid bucket.
            jlist = bdict[nkey]
            for i in ilist, j in jlist
                # Avoid double-computations.
                i < j || continue

                # Skip pairs that are not within distance.
                distance( L, r, α, xlist[[i,j]], ylist[[i,j]], θlist[[i,j]] ) || continue

                # Update index of connections list.
                push!( rlist, i );  push!( clist, j )
                push!( rlist, j );  push!( clist, i )
            end
        end
    end

    # Return the adjacency graph as a space matrix variables.
    return sparse( rlist, clist, ones( defInt, length( rlist ) ), N, N )
end

function erdosrenyiadjacency(N::defInt, L::defFloat, r::defFloat,
    xlist::Vector{defFloat}, ylist::Vector{defFloat}, θlist::Vector{defFloat},
    k::defFloat=1.0)
    # Initialize list of column/row indeces for connection tracking.
    rlist, clist = defInt[], defInt[]

    # For each agent pair, form connects with probability k/N.
    for i ∈ 1:N-1
        for j ∈ i+1:N
            if rand() < k/N
                # Update index of connections list.
                push!( rlist, i );  push!( clist, j )
                push!( rlist, j );  push!( clist, i )
            end
        end
    end

    # Use simple probability to compute adjacency.
    return sparse( rlist, clist, ones( defInt, length( rlist ) ), N, N )
end

function intadjacency(N::defInt, L::defFloat, r::defFloat, α::defFloat,
    xdata::Matrix{defFloat}, ydata::Matrix{defFloat}, θdata::Matrix{defFloat};
    distance=radialdistance)::SparseMatrixCSC{defInt,defInt}
    # Length of time to consider.
    tlist = 1:size( xdata, 1 )

    # Compute the integrated network in the time frame of the first peak.
    Â = spzeros( defInt, N, N )
    for t ∈ tlist
        # Compute the instantaneous adjacency.
        A = proximity( N, L, r, α, xdata[t,:], ydata[t,:], θdata[t,:]; distance=distance )

        # Update integrated adjacency matrix.
        Â = max.( Â, A )
    end

    # Return integrated networks.
    return Â
end

function conditionaladj(N::defInt, A::SparseMatrixCSC{defInt,defInt}, z::State)::SparseMatrixCSC{defInt,defInt}
    # Extract list of active agents.
    ia = findall( z.φ .== 0 )
    ib = findall( z.φ .== 1 )

    # Return the connections which are also active.
    Â = spzeros( defInt, N, N )
    Â[ia,ib] = A[ia,ib];  Â[ib,ia] = A[ib,ia]
    return Â
end

function conditionaladj(N::defInt, A::SparseMatrixCSC{defInt,defInt}, z::State; φ::defInt=0
    )::SparseMatrixCSC{defInt,defInt}
    # Extract list of active agents.
    i = findall( z.φ .== φ )

    # Return the connections by active agents.
    Â = spzeros( defInt, N, N )
    Â[i,:] .= A[i,:]
    return Â
end

function dirconditionaladj(N::defInt, A::SparseMatrixCSC{defInt,defInt}, z::State; φ1::defInt=0, φ2::defInt=1
    )::SparseMatrixCSC{defInt,defInt}
    # Extract list of active and inactive agents.
    i1 = findall( z.φ .== φ1 )
    i2 = findall( z.φ .== φ2 )

    # Return only the connections between active and inactive agents.
    Â = spzeros( defInt, N, N )
    Â[i1,i2] = A[i1,i2]
    return Â
end

function conditionaladj_3(N::defInt, A::SparseMatrixCSC{defInt,defInt}, z::State; φ1::defInt=0, φ2::defInt=1
    )::SparseMatrixCSC{defInt,defInt}
    # Extract list of active and inactive agents.
    i1 = findall( z.φ .== φ1 )
    i2 = findall( z.φ .== φ2 )

    # Return only the connections between active and inactive agents.
    Â = spzeros( defInt, N, N )
    Â[i1,i2] = A[i1,i2];  Â[i2,i1] = A[i2,i1]
    return Â
end

# Compute the depth of a node in a tree (assumes tree structure).
function depth(i::defInt, A::SparseMatrixCSC{defInt,defInt};
    memo::Dict{Int,Int}=Dict{Int,Int}(), f=minimum)::defInt
    # Return cached value if already computed
    haskey( memo, i ) && return memo[i]

    # Extract parents and initialize parent depth variable.
    jlist = findall( j -> j == 1, A[:,i] )

    # If empty, return 1, else, return the maximum memo depth. Update memo variable.
    d = length( jlist ) == 0 ? 1 : f( depth( j, A; memo=memo ) for j ∈ jlist ) + 1
    memo[i] = d

    # Return depth.
    return d
end

function depthN(N::defInt, A::SparseMatrixCSC{defInt,defInt}; f=minimum)::Vector{defInt}
    # Initialize memo variable.
    memo = Dict{defInt,defInt}()
    return [depth( i, A; memo=memo, f=f ) for i ∈ 1:N]
end

println( "Loaded network.jl class file." )
