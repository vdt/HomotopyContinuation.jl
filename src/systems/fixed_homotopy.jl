export FixedHomotopy

"""
    FixedHomotopy(H, t) <: AbstractSystem

Fix a homotopy `H(x,t)` at `t`
"""
struct FixedHomotopy{Hom<:AbstractHomotopy, T} <: AbstractSystem
    H::Hom
    t::T
end

struct FixedHomotopyCache{HC<:AbstractHomotopyCache} <: AbstractSystemCache
    cache::HC
end

cache(FH::FixedHomotopy, x) = FixedHomotopyCache(cache(FH.H, x, FH.t))

Base.size(F::FixedHomotopy) = size(F.H)

function evaluate!(u, F::FixedHomotopy, x, c::FixedHomotopyCache)
    evaluate!(u, F.H, x, F.t, c.cache)
end
evaluate(F::FixedHomotopy, x, c::FixedHomotopyCache) = evaluate(F.H, x, F.t, c.cache)
function jacobian!(U, F::FixedHomotopy, x, c::FixedHomotopyCache)
    jacobian!(U, F.H, x, F.t, c.cache)
end
function jacobian(F::FixedHomotopy, x, c::FixedHomotopyCache)
    jacobian(F.H, x, F.t, c.cache)
end
function evaluate_and_jacobian!(u, U, F::FixedHomotopy, x, c::FixedHomotopyCache)
    evaluate_and_jacobian!(u, U, F.H, x, F.t, c.cache)
end
function evaluate_and_jacobian(F::FixedHomotopy, x, c::FixedHomotopyCache)
    evaluate_and_jacobian(F.H, x, F.t, c.cache)
end
