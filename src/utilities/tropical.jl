using LinearAlgebra

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
    TropicalEvaluationResult{T}

Containing the smallest two evaluated values and its indices.

## Fields
* `min₁::T`
* `min₂::T`
* `i₁::Int`
* `i₂::Int`
"""
struct TropicalEvaluationResult{T}
    min₁::T
    min₂::T
    i₁::Int
    i₂::Int
end
Base.show(io::IO, ::MIME"application/prs.juno.inline", x::TropicalEvaluationResult) = x
function Base.show(io::IO, x::TropicalEvaluationResult)
    print(io, "(min₁: $(x.min₁), min₂: $(x.min₂), i₁: $(x.i₁), i₂: $(x.i₂))")
end

function sortevaluate(T::TropicalPolynomialSystem, w, m=nothing)
    map(T.polys) do Tᵢ

        vals = map(1:size(Tᵢ.exponents, 2)) do i
            evaluate_term(Tᵢ, w, m, i)
        end
        p = sortperm(vals)
        vals[p], p
    end
end

function evaluate(T::TropicalPolynomialSystem, w, m=nothing)
    out = Vector{TropicalEvaluationResult}(undef, length(T))
    evaluate!(out, T, w, m)
end
function evaluate!(out, T::TropicalPolynomialSystem, w, m=nothing)
    for (i, Tᵢ) in enumerate(T.polys)
        out[i] = evaluate(Tᵢ, w, m)
    end
    out
end

function evaluate(T::TropicalPolynomial, w::AbstractVector, m=nothing)
    v₁ = evaluate_term(T, w, m, 1)
    v₂ = evaluate_term(T, w, m, 2)

    if v₁ < v₂
        min₁, min₂ = v₁, v₂
        i₁, i₂ = 1, 2
    else
        min₁, min₂ = v₂, v₁
        i₁, i₂ = 2, 1
    end
    nterms = size(T.exponents, 2)
    for j in 3:nterms
        vⱼ = evaluate_term(T, w, m, j)
        if vⱼ < min₁
            min₁, min₂ = vⱼ, min₁
            i₁, i₂ = j, i₁
        elseif vⱼ < min₂
            min₂ = vⱼ
            i₂ = j
        end
    end
    TropicalEvaluationResult(min₁, min₂, i₁, i₂)
end

"""
    evaluate_term(T::TropicalPolynomial, w, j)
Compute the j-th term in the minimum of `T(w)`.
"""
Base.@propagate_inbounds function evaluate_term(T::TropicalPolynomial, w, j)
    vⱼ = T.weights[j] + T.exponents[1, j] * w[1]
    for i in 2:size(T.exponents, 1)
        vⱼ += T.exponents[i,j] * w[i]
    end
    vⱼ
end

"""
    evaluate_term(T::TropicalPolynomial, w, m, j)

Compute the j-th term in the minimum of `T(w//m)`. The result is scaled by `m`.
"""
Base.@propagate_inbounds function evaluate_term(T::TropicalPolynomial, w, m::Integer, j)
    vⱼ = T.weights[j] * m + T.exponents[1, j] * w[1]
    for i in 2:size(T.exponents, 1)
        vⱼ += T.exponents[i,j] * w[i]
    end
    vⱼ
end
evaluate_term(T::TropicalPolynomial, w, ::Nothing, j) = evaluate_term(T, w, j)


function is_zero(T::TropicalPolynomialSystem, w::AbstractVector, m=nothing)
    for Tᵢ in T.polys
        is_zero(Tᵢ, w, m) || return false
    end
    true
end
function is_zero(T::TropicalPolynomial, w::AbstractVector, m=nothing)
    is_zero(evaluate(T, w, m))
end
is_zero(R::TropicalEvaluationResult) = R.min₁ == R.min₂

function is_approximate_zero(T::TropicalPolynomialSystem, w::AbstractVector; tol=1e-1)
    for Tᵢ in T.polys
        is_approximate_zero(T, w, tol) || return false
    end
    true
end
function is_approximate_zero(Rs::Vector{<:TropicalEvaluationResult}; tol=1e-1)
    for R in Rs
        is_approximate_zero(R, tol) || return false
    end
    true
end
function is_approximate_zero(T::TropicalPolynomial, w::AbstractVector, tol)
    is_approximate_zero(evaluate(T, w, m), tol)
end
function is_approximate_zero(R::TropicalEvaluationResult, tol)
    R.min₂ - R.min₁ ≤ tol
end

function has_totaldegree_term(T::TropicalPolynomial, R::TropicalEvaluationResult)
    if T.weights[R.i₁] == 1
        return true
    elseif T.weights[R.i₂] == 1
        return true
    end
    return false
end
function has_totaldegree_term(T::TropicalPolynomialSystem, R::Vector{<:TropicalEvaluationResult})
    for (Tᵢ, Rᵢ) in zip(T.polys, R)
        has_totaldegree_term(Tᵢ, Rᵢ) && return true
    end
    false
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

function affine_initial_system(T::TropicalPolynomialSystem, approximation_result)
    A = zeros(Int, length(T.polys), size(T.polys[1].exponents, 1) - 1)
    b = zeros(Int, length(T.polys))
    affine_initial_system!(A, b, T, approximation_result)
    A, b
end
function affine_initial_system!(A, b, T::TropicalPolynomialSystem, approximation_result)
    j = 1
    for k = 1:length(T.polys)
        p = T.polys[k]
        c₁ = approximation_result[k].i₁
        c₂ = approximation_result[k].i₂
        h = p.exponents[end, c₁] - p.exponents[end, c₂]
        for i in 1:size(A, 2)
            A[k, i] = p.exponents[i, c₁] - p.exponents[i, c₂] + h
        end
        b[k] = p.weights[c₂] - p.weights[c₁]
    end
    nothing
end

function is_zero_of_system(A, b, w, m)
    for i=1:size(A,1)
        rᵢ = - m * b[i]
        for j=1:size(A, 2)
            rᵢ += A[i,j] * w[j]
        end
        rᵢ == 0 || return false
    end
    true
end
function best_w_m(T, approx_results, v, m_max)
    A, b = initial_system(T, approx_results)
    best_w = nothing
    best_m = -1
    best_m_diff = Inf
    for m=1:m_max
        # @show m
        w = round.(Int, m .* v)
        # @show evaluate(T, w .// m)
        # Check 1 verify that w is a zero of A w = m b
        # Check 2: verity that it is also a zero of trop(H)
        system_zero = is_zero_of_system(A, b, w, m)
        if system_zero
            @show m
            best_m = max(best_m, 0)
        end
        w_m_is_zero = system_zero && is_zero(T, w, m)

        if w_m_is_zero
            m_diff = 0.0
            for i=1:length(v)
                m_diff += abs(w[i] / m - v[i])
            end
            if m_diff < best_m_diff
                best_w = w
                best_m_diff = m_diff
                best_m = m
            end
        end
    end
    best_w, best_m, best_m_diff
end


function correct_val(T, approx_results, v)
    A, b = initial_system(T, approx_results)
    v′ = orthogonal_projection(v, A, b)
    v′ .-= minimum(v′)
end

function orthogonal_projection(v, A, b, m)
    p = qr(A, Val(true)) \ (m .* b)
    p + orthogonal_projection(m .* v - p, A)
end

function orthogonal_projection(v, A, b)
    p = qr(A, Val(true)) \ b
    p + orthogonal_projection(v - p, A)
end

function orthogonal_projection(v, A)
    N = nullspace(A)
    x = (N' * N) \ (N' * v)
    N * x
end
