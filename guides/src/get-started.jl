# **This guide assumes that you have Julia 1.0 installed and running.**

# If this is not the case follow the instructions at [julialang.org](https://julialang.org/downloads/). <h3 class="section-head" id="h-installation"><a href="#h-installation">Installation</a></h3>
#
# HomotopyContinuation.jl is available through the Julia package manager by
# ```julia-repl
# pkg> add HomotopyContinuation
# ```
# you can enter the Julia package manager by pressing `]` in the REPL.
#
#
# <h3 class="section-head" id="h-basic-usage"><a href="#h-basic-usage">Basic usage</a></h3>
# HomotopyContinuation.jl aims at having easy-to-understand top-level commands.
# For instance, suppose we want to solve the following polynomial system
# ```math
# f=\begin{bmatrix}x^2+2y \\\\ y^2-2 \end{bmatrix}
# ```
#
#
# This can be accomplished as follows
#
using HomotopyContinuation
@polyvar x y; # declare the variables x and y
#md solve([x^2+2y, y^2-2], report_progress=false) #hide
result = solve([x^2+2y, y^2-2])

#
# Let us see what is the information that we get. Four paths were attempted to be solved, four of which were completed successfully. Since we tried to solve an affine system, the algorithm checks whether there are solutions at infinity: in this case there are none. With *solutions at infinity* we mean solutions of the [homogenization](https://en.wikipedia.org/wiki/Homogeneous_polynomial#Homogenization) of the system which are no solutions of the affine system. None of the solutions are singular and two of them are real. To access the first solution in the array we write
#

## Assume we are only interested in the *real* solutions. We can obtain these by
real(result)

# where we can see that there are 2 real solutions, ``(2^{\frac34},-\sqrt{2})`` and ``(-2^{\frac34}, -\sqrt{2})``. We also can look into more detail into the first result by
#

real(result)[1]

# The returncode tells us that the pathtracking was successfull.
# What do the other entries of that table tell us? Let us consider the most relevant:
#
#   * `solution`: the zero that is computed (here it is ``(2^{\frac34},-\sqrt{2})``).
#   * `singular`: boolean that shows whether the zero is singular.
#   * `residual`: the computed value of ``|f([1.68179+1.27738e-17im, -1.41421-1.18454e-17im])|``.
#   * `windingnumber`: This is a lower bound on the multiplicity of the solution. A windingnumber greater than 1 indicates that the solution is singular.
#
# To extract the solutions you can do
#

solutions(result)

#
# or if you only want real solutions
#

realsolutions(result)
