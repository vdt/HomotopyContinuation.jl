using LinearAlgebra, HomotopyContinuation, PolynomialTestSystems
const HC = HomotopyContinuation
import MultivariatePolynomials
const MP = MultivariatePolynomials

f = equations(cyclic(7))


P, starts = pathsolver_startsolutions(f, seed=130793)

track(P, starts[2], 1.0, 1e-12)
track(P, starts[12], 1.0, 0.0001)

track(P, starts[111], 1.0, 1e-12)
track(P, starts[222], 1.0, 0)

A, b = HC.initial_system(P.tropical_system, P.state.tropical_approximation_results)

HC.best_w_m(P.tropical_system,
       P.state.tropical_approximation_results,
       P.state.val,
       15)

P.tropical_system




rank(A)

orthogonal_projection(P.state.val, A, b, 22)









P.state.val .* 8




@which nullspace(Rational.(A))
QR = qr(A)
QR.R
QR.Q


@which QR \ b

out = zeros(8)
ldiv!(out, QR, b)

A \ (7 .* b)
v = [5.29693e-5, 6.86643e-5, 0.99931, -0.000273215, -0.000240976, 0.999717, 0.999532, 0.428303]

nullspace(A)

HC.evaluate(P.tropical_system, P.state.val)

function orthogonal_projection(v, A, b, m)
    p = A \ (m .* b)
    p + orthogonal_projection(m .* v - p, A)
end

function orthogonal_projection(v, A)
    N = nullspace(A)
    x = (N' * N) \ (N' * v)
    N * x
end

orthogonal_projection(v, A, b, 7)

nullspace(A)


lu_A.L
lu_A.U
lu_A.P
b

PathSolving.update_tropical_system!(P)
P.state



g = equations(bacillus_subtilis())
Q, starts2 = pathsolver_startsolutions(g, seed=130793)
track(Q, starts2[2], 1.0, 3[]e-5)

A, b = HC.initial_system(Q.tropical_system, Q.state.tropical_approximation_results)

lu(Rational.(A), check=false)
S = nullspace(Rational.(A))
p = A \ b
v = [0.475455, 0.00165342, -0.00130019, 1.00497, 0.00375947, -0.00045538, 2.00403, 1.00625, 1.00116, 2.0027, 1.00114]

r = [S A'] \ (v - p)

p + S * r[1:2]





using RowEchelon


rref_with_pivots((A))
rank(A)

A * [1//2, 0, 0, 1, 0, 0, 2, 1, 1, 2, 1] - b
