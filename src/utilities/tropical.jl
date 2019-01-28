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


function approximate_evaluate(T::TropicalPolynomialSystem, w)
    out = Vector{NTuple{2, eltype(w)}}(undef, length(T))
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
    for j in 1:n
        vⱼ = S(T.weights[j])
        for i in 1:m
            vⱼ += T.exponents[i,j] * w[i]
        end
        if vⱼ < min₁
            min₁, min₂ = vⱼ, min₁
        elseif vⱼ < min₂
            min₂ = vⱼ
        end
    end
    min₁, min₂
end
