using HomotopyContinuation, PolynomialTestSystems, ProjectiveVectors
using LinearAlgebra
import MultivariatePolynomials
const MP = MultivariatePolynomials

f = equations(cyclic(7))
# solve(f, seed=130793, threading=false)

P, starts = PathSolving.pathsolver_startsolutions(f, seed=130793)
PathSolving.track(P, starts[12], 1.0, 0.0001)
PathSolving.track(P, starts[111], 1.0, 1e-12)
PathSolving.track(P, starts[222], 1.0, 0)

PathSolving.update_tropical_system!(P)
P.state







P






















x = [-3.9756e-12, -3.77352e-12, 1.0, -6.6696e-12, -6.40806e-12, 1.0, 1.0, 0.42754]



P


P.tracker.state.cond


g = equations(bacillus_subtilis())


Q, Q_starts = PathSolving.pathsolver_startsolutions(g, seed=869148)
PathSolving.track(Q, Q_starts[1194], 1.0, 0.0)

solve(f, patch=AffinePatches.ScaledOrthogonalPatch())



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
function exponent_matrices(f; parameters=nothing)
    vars = Utilities.variables(f; parameters=parameters)
    # If there are parameters we can have duplicate exponent matrices
    map(fᵢ -> exponent_matrix(fᵢ, vars), f)
end


tropical_evaluation(EM::Vector{<:Matrix}, x) = tropical_evaluation.(EM, Ref(x))
function tropical_evaluation(E::Matrix, x)
    m, n = size(E)
    max₁ = max₂ = -Inf

    @inbounds for j in 1:n
        vⱼ = 0.0
        for i in 1:m
            vⱼ += E[i,j] * x[i]
        end
        if vⱼ > max₁
            max₂ = max₁
            max₁ = vⱼ
        elseif vⱼ > max₂
            max₂ = vⱼ
        end
    end
    max₁, max₂
end

EM = exponent_matrices(Utilities.homogenize(f))
x = [7.08972e-7, 0.333331, 0.333332, 0.333333, -7.09829e-7, 0.666665, 0.666664, 0.333191]
tropical_evaluation(EM, x)





tracker, starts = żtartsolutions(f, seed=869148, tol=1e-5,
    patch=AffinePatches.ScaledOrthogonalPatch(), predictor=Predictors.Heun())
S = collect(starts)
track(tracker, S[1194], 1.0, 0.0)


tracker.state






tracker.state

tracker.cache.Jac


R = PathTracking.track(tracker, S[23], 1.0, 0.95)


S = collect(starts)

for (k, s) in enumerate(S)
    R = PathTracking.track(tracker, s, 1.0, 0.0)
    if R.returncode == PathTracking.Status.success
        @show k
    end
end




Utilities.disable_rowscaling!(tracker.cache.Jac)
tracker.cache
4.925646704945219e-5 + -4.119487552933787e-5

1.8722276032945295e-5 + -1.8712276032895048e-5
1.8673383000988864e-6 + -1.8573383000486388e-6
0.00016023730934799918 + -6.455605280897281e-5
tracker.state.x.data + t * tracker.state.ẋ
y = ProjectiveVectors.affine_chart(tracker.state.x)

PathTracking.track(tracker, y, 0.001, 0.0)

ws = begin
    x = tracker.state.x
    ẋ = tracker.state.ẋ
    t * (real.(x .* conj.(ẋ)) ./ abs2.(x))
end

rationalize.(ws, tol=2//50)


y0001 = ProjectiveVectors.affine_chart(tracker.state.x)
y001 = ProjectiveVectors.affine_chart(tracker.state.x)
y01 = ProjectiveVectors.affine_chart(tracker.state.x)
tracker.state.ẋ


@polyvar x y z
f = (x-1*y)*(x-1.01*y)
g = x - y - z
tracker, starts = żtartsolutions([f, g], predictor=Predictors.Heun())
S = collect(starts)


r1 = PathTracking.track(tracker, S[1], 1.0, 0.0)
r2 = PathTracking.track(tracker, S[2], 1.0, 0.0)





@show x1 x2

solve([x - y])


@polyvar x y z

f = [x^2+y^2+z^2, 2x-y+4z]
Problems.problem_startsolutions(f, homvar=y)



tracker, starts = pathtracker_startsolutions(f, homvar=y)


PathTracking.track(tracker, first(starts), 1.0, 0.2)









ProjectiveVectors.embed!(tracker.state.x, first(starts))






typeof(tracker.state.x)

using HomotopyContinuation, PolynomialTestSystems#, HCTools
using LinearAlgebra, FixedPolynomials


A = rand(5, 3)
b = rand(5)
x = zeros(3)

ldiv!(qr(A), copy(b))


qr(A) \ b

tracker, starts = pathtracker_startsolutions([f; L₁], [f; L₂], [[1, 1, 1, 1]])

PathTracking.track(tracker, first(starts), 1.0, 0.0)


tracker.state



F = equations(cyclic(5))


solve(F, seed=405709, threading=false, corrector_maxiters=10, predictor=Predictors.Euler())
solver, starts = Solving.solver_startsolutions(F, seed=405709,
    predictor=Predictors.Euler(), maxiters=100_000, tol=5e-5)

S = collect(starts)

R = solve(solver, S, threading=false)
R[69]
solver.tracker.state.x.homvar
r1 = PathTracking.track(solver.tracker, S[23], 1.0, 0.5)


solve(F, seed=32421)


tracker, starts = pathtracker_startsolutions(F, seed=32421)
tracker

r1 = PathTracking.track(solver.tracker, S[29], 1.0, 0.1)
r2 = PathTracking.track(solver.tracker, S[47], 1.0, 0.1)

r2 = PathTracking.track(solver.tracker, S[17], 1.0, 0.1)
#

# julia> ProjectiveVectors.affine(r1.x)
# 5-element Array{Complex{Float64},1}:
#  -0.7624382933488744 - 0.9646173838511025im
#   0.8972052090854287 - 0.38549396294581384im
#  -0.8252912542049768 + 0.13519867071765057im
#   0.2124762303572757 - 1.2064444376492391im
#   0.3650447368234911 + 1.0264009950243023im
#
# julia> ProjectiveVectors.affine(r2.x)
# 5-element Array{Complex{Float64},1}:
#     0.05505395895765863 - 0.7638544607309543im
#      0.4646407326543956 - 0.26402447889706065im
#     -0.7263710532336596 - 0.944689888434731im
#      0.2903769794651848 + 0.08467072763300884im
#  -0.0007456614329823807 + 1.0458704946876713im





r2 = PathTracking.track(solver.tracker, S[8], 1.0, 0.4)

r1 = PathTracking.track(solver.tracker, S[71], 1.0, 0.1)
println("---")
r24 = PathTracking.track(solver.tracker, S[95], 1.0, 0.1)


ProjectiveVectors.affine(r1.x)
ProjectiveVectors.affine(r2.x)




map(r -> r.solution, atinfinity(R))



























































, BenchmarkTools
using Gadfly, LinearAlgebra

f = equations(rps10())
tracker, s = pathtracker_startsolutions(f, seed=5233,
        predictor=Predictors.Ralston(), corrector=Correctors.Newton(), maxiters=200)
starts = collect(s)


function plot_step_sizes(tracker, x, t₁::Real, t₀::Real; kwargs...)
    PathTracking.setup!(tracker, x, t₁, t₀; kwargs...)
    last_t_dt = nothing
    t_dt_values = Vector{NTuple{2, Float64}}()
    rejections = Vector{NTuple{2, Float64}}()
    t = real(PathTracking.currt(tracker))
    Δt = abs(PathTracking.currΔt(tracker))
    last_t_dt = (t, Δt)
    PathTracking.step!(tracker)
    for _ in tracker
        t = real(PathTracking.currt(tracker))
        Δt = abs(PathTracking.currΔt(tracker))
        if t == last_t_dt[1] # still same t -> last step was rejected
            push!(rejections, last_t_dt)
        else
            push!(t_dt_values, last_t_dt)
        end
        last_t_dt = (t, Δt)
    end

    l1 = layer(x=map(first, t_dt_values), y=map(last, t_dt_values), Stat.step(direction=:vh),
        Geom.line)
    l2 = layer(x=map(first, rejections), y=map(last, rejections), Geom.point, color=[colorant"red"], shape=[Shape.xcross])
    Gadfly.with_theme(:default) do
        plot(l2, l1, Guide.xlabel("t"), Guide.ylabel("|Δt|"))
    end
end

function trackall(tracker, starts)
    for (i, s) in enumerate(starts)
        # @show i
        retcode = PathTracking.track!(tracker, s, 1.0, 0.2)
        if retcode != PathTracking.done
            @show i, retcode
        end
    end
end

@code_typed trackall(tracker, starts)
@time trackall(tracker, starts[1:12])

plot_step_sizes(tracker, starts[123], 1.0, 0.2, initial_steplength=0.1)

PathTracking.track(tracker, starts[6], 1.0, 0.2, initial_steplength=0.10)
Homotopies.jacobian(tracker.cache.homotopy, tracker.state.x, 1.0)

x = starts[123]
@benchmark PathTracking.track!($tracker, $(x), $1.0, $0.1)


tracker.homotopy
using BenchmarkTools

v = rand(ComplexF64, 20)


function totuple(v)
    map(vᵢ -> round(real(vᵢ);  digits=5), v)
end

function round_array!(v)
    for vᵢ in v
         round(real(vᵢ)
    map(vᵢ -> round(real(vᵢ);  digits=5), v)
end


function generic_lufact!(A::StridedMatrix{T}, ::Val{Pivot} = Val(true);
                         check::Bool = true) where {T,Pivot}
    m, n = size(A)
    minmn = ifelse(m < n, m, n)
    info = 0
    ipiv = Vector{Int}(undef, minmn)
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if Pivot
                amax = zero(real(T))
                for i = k:m
                    absi = abs2(A[i,k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            end
            ipiv[k] = kp
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                end
                # Scale first column
                @fastmat Akkinv = inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            elseif info == 0
                info = k
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    check
    return LU{T,typeof(A)}(A, ipiv, convert(Int, info))
end

using BenchmarkTools

A = rand(ComplexF64, 12, 12)

@benchmark generic_lufact!(X) setup=(X=copy($A)) evals=1 samples=100_000







u, v = rand(ComplexF64, 2)
@benchmark (@fastmath $u * $v)

@btime totuple(v)

u = tuple(w...)

@btime hash($u)
@btime hash($w)
@btime hash($(float.(w)))
