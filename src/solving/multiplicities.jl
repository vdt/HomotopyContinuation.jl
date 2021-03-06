## The structure for sorting the solutions is a k search tree.
## First, we cluster the points by ordering them w.r.t. the the absolute value of their entries and w.r.t. the arguments of the entries.
## Then, we compare the points pairwise in each cluster.
mutable struct SearchTree
    node::Vector{Int}
    left::Union{Nothing, SearchTree}
    right::Union{Nothing, SearchTree}
end

SearchTree(i::Int) = SearchTree([i], nothing, nothing)


function cluster_by_norm!(current::SearchTree, i, vectors, nvars, τ)
    added = false
    j = 1
    while !added && j ≤ nvars
        Δ = abs(vectors[i][j]) - abs(vectors[current.node[1]][j])
        if Δ < (- τ)
            added = true
            if current.left === nothing
                current.left = SearchTree(i)
            else
                cluster_by_norm!(current.left, i, vectors, nvars, τ)
            end
        elseif Δ > τ
            added = true
            if current.right === nothing
                current.right = SearchTree(i)
            else
                cluster_by_norm!(current.right, i, vectors, nvars, τ)
            end
        end
        j += 1
    end

    if !added
        push!(current.node, i)
    end
end

function cluster_by_angle!(current::SearchTree, i, vectors, nvars, τ)
    added = false
    j = 1
    α = β = 0.0
    while j ≤ nvars
        if abs(vectors[i][j]) > 0.1 && abs(vectors[current.node[1]][j]) > 0.1
            α = angle(vectors[i][j])
            β = angle(vectors[current.node[1]][j])
            j+=1
            break
        end
        j+=1
    end

    while !added && j ≤ nvars
        if abs(vectors[i][j]) > 0.1 && abs(vectors[current.node[1]][j]) > 0.1
            Δ = abs(angle(vectors[i][j]) - α) - abs(angle(vectors[current.node[1]][j]) - β)
            if Δ < (- τ)
                added = true
                if current.left === nothing
                    current.left = SearchTree(i)
                else
                    cluster_by_angle!(current.left, i, vectors, nvars, τ)
                end
            elseif Δ > τ
                added = true
                if current.right === nothing
                    current.right = SearchTree(i)
                else
                    cluster_by_angle!(current.right, i, vectors, nvars, τ)
                end
            end
        end
        j += 1
    end

    if !added
        push!(current.node, i)
    end
end

function push_for_identifying_multiplicities!(current::SearchTree, i, vectors, τ, distance::F) where {F<:Function}
    # This compares with the distance function
    if distance(vectors[i], vectors[current.node[1]]) < τ
        push!(current.node, i)
    else
        if current.left === nothing
            current.left = SearchTree(i)
        else
            push_for_identifying_multiplicities!(current.left, i, vectors, τ, distance)
        end
    end
end


function Base.foreach(f::F, t::SearchTree) where {F<:Function}
    f(t.node)
    if t.left !== nothing
        foreach(f, t.left)
    end
    if t.right !== nothing
        foreach(f, t.right)
    end

    nothing
end


"""
    multiplicities(vectors, tol, distance)

Returns an array of arrays of integers. Each vector v in vectors contains all indices i,j such that V[i] and V[j] have distance at most tol.
"""
function multiplicities(vectors::Vector{<:AbstractVector{T}}, tol, distance::F) where {T<:Number, F<:Function}
    nvars = length(vectors[1])
    # root0 is for sorting by norm
    n = length(vectors)
    root0 = SearchTree(1)
    for i in 2:n
        cluster_by_norm!(root0, i, vectors, nvars, tol)
    end

    mults::Vector{Vector{Int}} = Vector{Vector{Int}}()
    foreach(root0) do m0
        if length(m0) == 1
            return nothing
        end

        #root1 is for sorting by argument.
        root1 = SearchTree(m0[1])
        for i in Iterators.drop(m0, 1)
            cluster_by_angle!(root1, i, vectors, nvars, tol)
        end

        foreach(root1) do m1
            if length(m1) == 1
                return nothing
            end

            #root2 is for sorting by Fubini-Study distance
            root2 = SearchTree(m1[1])
            for i in Iterators.drop(m1, 1)
                push_for_identifying_multiplicities!(root2, i, vectors, tol, distance)
            end
            foreach(root2) do m2
                if length(m2) == 1
                    return nothing
                end

                push!(mults, m2)
                nothing
            end
            nothing
        end
    end
    mults::Vector{Vector{Int}}
end
