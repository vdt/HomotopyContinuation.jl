using LinearAlgebra, HomotopyContinuation, PolynomialTestSystems
const HC = HomotopyContinuation
import MultivariatePolynomials
const MP = MultivariatePolynomials

f = equations(cyclic(7))


P, starts = pathsolver_startsolutions(f, seed=130793)
track(P, starts[12], 1.0, 0.0001)

track(P, starts[111], 1.0, 1e-12)
track(P, starts[222], 1.0, 0)

PathSolving.update_tropical_system!(P)
P.state



g = equations(bacillus_subtilis())
Q, starts2 = pathsolver_startsolutions(g, seed=130793)
track(Q, starts2[2], 1.0, 5e-6)

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
