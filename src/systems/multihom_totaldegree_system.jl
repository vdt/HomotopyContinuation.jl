export MultiHomogenousTotalDegreeSystem

"""
    MultiHomogenousTotalDegreeSystem(polynomials, vars) <: AbstractSystem

Create a tottal degree system
"""
struct MultiHomTotalDegreeSystem <: AbstractSystem
    D::Matrix{Int}
    C::Matrix{Float64}
end

struct MultiHomTotalDegreeSystemCache{N,T} <: AbstractSystem
    B::Matrix{T}
    Ĝ::Matrix{T}
    R::Matrix{T}
    S::Matrix{T}
    ranges::NTuple{N, UnitRange{Int}}
    homvars::NTuple{N, Int}
end

Base.size(F::MultiHomTotalDegreeSystem{N}) where N = (size(F.D, 2), sum(size(F.D)))

function cache(::MultiHomTotalDegreeSystem, x::PVector{M}) where M
    B =
end

function evaluate!(u, F::MultiHomTotalDegreeSystem, z::ProjectiveVectors.PVector{M}, cache::MultiHomTotalDegreeSystemCache{M}) where M
    D, C = F.D, F.C
    B, ranges, homvars = cache.B, cache.ranges, cache.homvars

    n = size(B, 1)
    # Compute all bᵢⱼ and store in B
    # Since B is column major we store bᵢⱼ in B[j, i]
    for i=1:n
        for j=1:M
            if D[j, i] != 0
                bᵢⱼ = -zero(B[j, i])
                for k in ranges[i]
                    bᵢⱼ += C[k, i] * z[k]
                end
                B[j, i] = bᵢⱼ
            end
        end
    end

    for i=1:n
        gᵢ = one(eltype(u))
        for j=1:M
            dᵢⱼ, bᵢⱼ = D[j, i], B[j, i]
            if dᵢⱼ != 0
                gᵢ *= bᵢⱼ^dᵢⱼ - z[homvars[j]]^dᵢⱼ
            end
        end
        u[i] = gᵢ
    end

    u
end
function evaluate(F::TotalDegreeSystem, x, cache::SystemNullCache)
    u = similar(x, size(F, 1))
    evaluate!(u, F, x, cache)
    u
end

function jacobian!(U, F::MultiHomTotalDegreeSystem, z::ProjectiveVectors.PVector{M}, cache::MultiHomTotalDegreeSystemCache{M}) where M
    evaluate_and_jacobian!(nothing, U, F, z, cache)
    U
end
function jacobian(F::TotalDegreeSystem, x, cache::SystemNullCache)
    U = similar(x, size(F))
    jacobian!(U, F, x, cache)
    U
end

function evaluate_and_jacobian!(u, U, F::TotalDegreeSystem, x, ::SystemNullCache)
    D, C = F.D, F.C
    B, Ĝ, R, S, ranges, homvars = cache.B, cache.Ĝ, cache.R, cache.S, cache.ranges, cache.S

    n = size(B, 1)

    for i=1:n
        # Compute all bᵢⱼ and store in B
        # Since B is column major we store bᵢⱼ in B[j, i]
        for j=1:M
            if D[j, i] != 0
                bᵢⱼ = -zero(B[j, i])
                for k in ranges[j]
                    bᵢⱼ += C[k, i] * z[k]
                end
                B[j, i] = bᵢⱼ
            end
        end

        # Compute all ĝᵢⱼ and store in Ĝ
        # Since Ĝ is column major we store ĝᵢⱼ in Ĝ[j, i]
        for j=1:M
            dᵢⱼ, bᵢⱼ = D[j, i], B[j, i]
            if dᵢⱼ == 0
                Ĝ[j, i] = one(eltype(Ĝ))
            else
                Ĝ[j, i] = bᵢⱼ^dᵢⱼ - z[homvars[j]]^dᵢⱼ
            end
        end

        # Accumulate subproducts forward
        R[1, i] = rᵢⱼ_prev = Ĝ[1, i]
        for j=2:M
            if D[j, i] != 0 # otherwise Ĝ[j, i] = 1
                r[j, i] = rᵢⱼ_prev = rᵢⱼ_prev * Ĝ[j, i]
            end
            if u !== nothing
                u[i] = rᵢⱼ_prev
            end
        end

        # Accumulate subproducts backward
        S[M, i] = sᵢⱼ_prev = Ĝ[M, i]
        for j=M-1:-1:1
            if D[j, i] != 0 # otherwise Ĝ[j, i] = 1
                S[j, i] = sᵢⱼ_prev = sᵢⱼ_prev * Ĝ[j, i]
            end
        end

        # Compute partial derivatives
        for j=1:M
            dᵢⱼ = D[j, i]
            for k in ranges[j]
                c = C[k, i]
                if iszero(c) || iszero(dᵢⱼ)
                    U[i, k] = zero(eltype(U))
                else
                    u_ik = dᵢⱼ * B[j,i]^(dᵢⱼ - 1)
                    if j > 1
                        u_ik *= R[j-1,i]
                    end
                    if j < M
                        u_ik *= S[j+1, i]
                    end
                    U[i, k] = u_ik
                end
            end
        end
    end
    nothing
end
