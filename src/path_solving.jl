module PathSolving

import MultivariatePolynomials
import ..Problems
import ..PathTracking

using ..Utilities

const MP = MultivariatePolynomials


struct TropicalSystem
    exponents::Vector{Matrix{Int16}}
end

function TropicalSystem(f::Vector{<:MP.AbstractPolynomialLike})
    vars = Utilities.variables(f)
    TropicalSystem(map(fᵢ -> exponent_matrix(fᵢ, vars), f))
end
function exponent_matrix(f::MP.AbstractPolynomialLike, variables)
    m, n = length(variables), MP.nterms(f)
    E = zeros(Int16, m, n)
    for (k, t) in enumerate(MP.terms(f))
        for i in 1:m
            E[i, k] = MP.degree(t, variables[i])
        end
    end
    E
end

function evaluate!(out, T::TropicalSystem, x)
    for (i, E) in enumerate(T.exponents)
        out[i] = tropical_evaluation(E, x)
    end
end
function tropical_evaluation(E::Matrix, x)
    m, n = size(E)
    max₁ = max₂ = -Inf
    min₁ = min₂ = Inf
    @inbounds for j in 1:n
        vⱼ = 0.0
        for i in 1:m
            vⱼ += E[i,j] * x[i]
        end
        if vⱼ > max₁
            max₁, max₂ = vⱼ, max₁
        elseif vⱼ > max₂
            max₂ = vⱼ
        end
        if vⱼ < min₁
            min₁, min₂ = vⱼ, min₁
        elseif vⱼ < min₂
            min₂ = vⱼ
        end
    end
    max₁, max₂, min₁, min₂
end


mutable struct State
    segment::ComplexSegment
    s::Float64
    val::Vector{Float64}
    tropical_evaluation::Vector{NTuple{4,Float64}}
end

function State(tracker, t₁::Number, t₀::Number)
    segment = ComplexSegment(t₁, t₀)
    s = 1.0
    val = zeros(length(tracker.state.x))
    tropical_evaluation = Vector{NTuple{4,Float64}}(undef, length(val) - 1)
    State(segment, s, val, tropical_evaluation)
end

struct PathSolver{Tracker<:PathTracking.PathTracker}
    tracker::Tracker
    state::State
    tropical_system::TropicalSystem
end

function PathSolver(prob::Problems.Projective, x; kwargs...)
    tracker = PathTracking.PathTracker(prob, x, 1.0, 0.0; kwargs...)
    state = State(tracker, 1.0, 0.0)
    tropical_system = TropicalSystem(prob.target_system)
    PathSolver(tracker, state, tropical_system)
end


function pathsolver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, Problems.supported_keywords)
    prob, startsolutions = Problems.problem_startsolutions(args...; supported...)
    solver = PathSolver(prob, Utilities.start_solution_sample(startsolutions); rest...)
    solver, collect(startsolutions)
end

function track(solver::PathSolver, x, t₁, t₀)
    tracker = solver.tracker

    PathTracking.setup!(tracker, x, t₁, t₀)

    while tracker.state.status == PathTracking.Status.tracking
        PathTracking.step!(tracker)
        PathTracking.check_terminated!(tracker)
        t = real(PathTracking.currt(solver.tracker))
        if !tracker.state.last_step_failed
            update_valuation!(solver)
            update_tropical_system!(solver)
            println(t)
            println(solver.state.val)
            display(solver.state.tropical_evaluation)
        end

        if tracker.state.cond > 1e12
            tracker.state.status = PathTracking.Status.terminated_singularity
        end
    end
    if tracker.state.status == PathTracking.Status.success
        PathTracking.refine!(tracker)
    end
    PathTracking.PathTrackerResult(tracker)
end

function valuation!(val, z, ż, t)
    for i in eachindex(ż)
        val[i] = t * re_dot(z[i], ż[i]) / abs2(z[i])
    end
    val
end

"Compute `real(z * conj(ż))`."
re_dot(z, ż) = begin x, y = reim(z); ẋ, ẏ = reim(ż); x * ẋ + y * ẏ end

function update_valuation!(solver::PathSolver)
    t = real(PathTracking.currt(solver.tracker))
    valuation!(solver.state.val, solver.tracker.state.x, solver.tracker.state.ẋ, t)
end

function update_tropical_system!(solver::PathSolver)
    evaluate!(solver.state.tropical_evaluation, solver.tropical_system, solver.state.val)
end



end
