struct TropicalPolynomial{T}
    exponents::Matrix{Int32}
    weights::Vector{T}
end

struct TropicalPolynomialSystem{T}
    polys::Vector{TropicalPolynomial{T}}
end
Base.length(T::TropicalPolynomialSystem) = length(T.polys)

function tropicalize_straight_line(g::MPPoly, f::MPPoly, variables)
    m, n₁, n₂ = length(variables), MP.nterms(f), MP.nterms(g)
    E = zeros(Int32, m, n₁ + n₂)
    weights = zeros(Int32, n₁ + n₂)
    for (k, t) in enumerate(MP.terms(f)), i in 1:m
        E[i, k] = MP.degree(t, variables[i])
        weights[k] = Int32(0)
    end
    for (k, t) in enumerate(MP.terms(g)), i in 1:m
        E[i, k + n₁] = MP.degree(t, variables[i])
        weights[k + n₁] = Int32(1)
    end
    TropicalPolynomial(E, weights)
end
function tropicalize_straight_line(g, f)
    vars = variables(f)
    polys = map(expand(g), expand(f)) do gᵢ, fᵢ
        tropicalize_straight_line(gᵢ, fᵢ, vars)
    end
    TropicalPolynomialSystem(polys)
end

"""
    ApproximateResult

Containing the smallest two evaluated values and its indices.

## Fields
* `min₁::Float64`
* `min₂::Float64`
* `i₁::Int`
* `i₂::Int`
"""
struct TropicalApproximationResult
    min₁::Float64
    min₂::Float64
    i₁::Int
    i₂::Int
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::TropicalApproximationResult) = x
function Base.show(io::IO, x::TropicalApproximationResult)
    print(io, "(min₁: $(x.min₁), min₂: $(x.min₂), i₁: $(x.i₁), i₂: $(x.i₂))")
end

function approximate_evaluate(T::TropicalPolynomialSystem, w)
    out = Vector{TropicalApproximationResult}(undef, length(T))
    approximate_evaluate!(out, T, w)
end

function approximate_evaluate!(out, T::TropicalPolynomialSystem, w)
    for (i, Tᵢ) in enumerate(T.polys)
        out[i] = approximate_evaluate(Tᵢ, w)
    end
    out
end

function approximate_evaluate(T::TropicalPolynomial, w::AbstractVector{S}) where {S<:AbstractFloat}
    m, n = size(T.exponents)
    min₁ = min₂ = Inf
    i₁ = i₂ = 1
    for j in 1:n
        vⱼ = S(T.weights[j])
        for i in 1:m
            vⱼ += T.exponents[i,j] * w[i]
        end
        if vⱼ < min₁
            min₁, min₂ = vⱼ, min₁
            i₁, i₂ = j, i₁
        elseif vⱼ < min₂
            min₂ = vⱼ
            i₂ = j
        end
    end
    TropicalApproximationResult(min₁, min₂, i₁, i₂)
end

function initial_system(T::TropicalPolynomialSystem, approximation_result)
    A = zeros(Int, length(T.polys), size(T.polys[1].exponents, 1))
    b = zeros(Int, length(T.polys))
    initial_system!(A, b, T, approximation_result)
    A, b
end
function initial_system!(A, b, T::TropicalPolynomialSystem, approximation_result)
    j = 1
    for k = 1:length(T.polys)
        p = T.polys[k]
        c₁ = approximation_result[k].i₁
        c₂ = approximation_result[k].i₂
        for i in 1:size(A, 2)
            A[k, i] = p.exponents[i, c₁] - p.exponents[i, c₂]
        end
        b[k] = p.weights[c₂] - p.weights[c₁]
    end
    nothing
end
