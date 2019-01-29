export pathsolver_startsolutions

mutable struct PathSolverState
    segment::ComplexSegment
    s::Float64
    val::Vector{Float64}
    tropical_approximation_results::Vector{TropicalApproximationResult}
end

function PathSolverState(tracker, tropical_system, t₁::Number, t₀::Number)
    segment = ComplexSegment(t₁, t₀)
    s = 1.0
    val = zeros(length(tracker.state.x))
    tropical_approximation_results = approximate_evaluate(tropical_system, val)
    PathSolverState(segment, s, val, tropical_approximation_results)
end

struct PathSolver{Tracker<:PathTracker}
    tracker::Tracker
    state::PathSolverState
    tropical_system::Union{Nothing, TropicalPolynomialSystem{Int32}}
end

function PathSolver(prob::ProjectiveProblem, x; kwargs...)
    tracker = PathTracker(prob, x, 1.0, 0.0; kwargs...)
    state = PathSolverState(tracker, prob.tropical_system, 1.0, 0.0)
    PathSolver(tracker, state, prob.tropical_system)
end


function pathsolver_startsolutions(args...; kwargs...)
    supported, rest = splitkwargs(kwargs, problem_startsolutions_supported_keywords)
    prob, startsolutions = problem_startsolutions(args...; supported...)
    solver = PathSolver(prob, start_solution_sample(startsolutions); rest...)
    solver, collect(startsolutions)
end

function track(solver::PathSolver, x, t₁, t₀)
    tracker = solver.tracker

    setup!(tracker, x, t₁, t₀)

    while tracker.state.status == PathTrackerStatus.tracking
        step!(tracker)
        check_terminated!(tracker)
        t = real(currt(solver.tracker))
        if !tracker.state.last_step_failed
            update_valuation!(solver)
            update_tropical_system!(solver)
            println(t)
            println(solver.state.val)
            display(solver.state.tropical_approximation_results)
        end

        if tracker.state.cond > 1e12
            tracker.state.status = PathTrackerStatus.terminated_singularity
        end
    end
    if tracker.state.status == PathTrackerStatus.success
        refine!(tracker)
    end
    PathTrackerResult(tracker)
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
    t = real(currt(solver.tracker))
    valuation!(solver.state.val, solver.tracker.state.x, solver.tracker.state.ẋ, t)
end

function update_tropical_system!(solver::PathSolver)
    approximate_evaluate!(solver.state.tropical_approximation_results, solver.tropical_system, solver.state.val)
end
