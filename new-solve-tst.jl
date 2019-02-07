using LinearAlgebra, HomotopyContinuation, PolynomialTestSystems
const HC = HomotopyContinuation
import MultivariatePolynomials
const MP = MultivariatePolynomials

f = equations(cyclic(7))

P, starts = pathsolver_startsolutions(f, seed=130793)

R = map(starts) do s
    track(P, s, 1.0, 0.0, debug=false)
end

count(r -> r[2] == :failed, R)

track(P, starts[9], 1.0, 1e-3)
track(P, starts[57], 1.0, 0.0)
using AbstractAlgebra

C, d = HC.affine_initial_system(P.tropical_system, P.state.tropical_approximation_results)
A, b = HC.initial_system(P.tropical_system, P.state.tropical_approximation_results)
v = [0.0, 0.997296, 0.997062, 0.00356295, 0.0, 0.994626, 0.0, 0.138842]
v = [0.0359355, 0.0164393, 0.815228, -0.0687553, -0.0384354, 0.926298, 0.829523, 0.357542]
v2 = rationalize.(v, tol=1/1000)

orthogonal_projectio- 0.0192

aavec(v) = matrix(ZZ, reshape(v, length(v), 1))
A1 = matrix(ZZ, A[2:end,:])
b1 = matrix(ZZ, reshape(b[2:end], length(b)-1, 1))
S, U, V = snf_with_trafo(A1)
c = U * b1
S_inv, denom = inv(S[1:6, 1:6])

x1 = V * aavec([[div((S_inv * c)[i, 1], denom)  for i=1:size(c, 1)]; 0; 0])


H, U = hnf_cohen_with_trafo(A1)

bH
hnf_with_trafo(A1)




rank(A[2:end,:])

res = solve(f, seed=130793)

findall(map(r -> r.returncode == :at_infinity, res) .!= map(r -> r[2] == :at_infinity, R))

track(P, starts[24], 1.0, 1e-12)
track(P, starts[12], 1.0, 0.0001)

track(P, starts[111], 1.0, 1e-12)
track(P, starts[222], 1.0, 0)

A, b = HC.initial_system(P.tropical_system, P.state.tropical_approximation_results)

HC.best_w_m(P.tropical_system,
       P.state.tropical_approximation_results,
       P.state.val,
       15)

P.tropical_system

v = [0.48444, 0.00460564, 0.966807, -0.00667572, 0.731983, 0.116612, 0.369713, 0.381069]

v .* 64


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

v2 = [0.480437, 0.0147987, 0.840107, -0.0715342, 0.644538, 0.0472122, 0.212632, 0.310643]
v3 = [0.464675, 0.00947736, 0.926048, -0.0228519, 0.701272, 0.098716, 0.348311, 0.360828]
nullspace(A)

orthogonal_projection(v2, A, b)

64 .* orthogonal_projection(v3, A, b)

A


HC.evaluate(P.tropical_system, P.state.val)


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
