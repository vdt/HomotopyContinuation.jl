module HomotopyContinuation

    import DynamicPolynomials
    import FixedPolynomials
    import LinearAlgebra
    import MultivariatePolynomials
    import Printf
    import ProgressMeter
    import ProjectiveVectors
    import Random
    import StaticArrays
    import StaticPolynomials
    import TreeViews

    import DoubleFloats: Double64
    import DynamicPolynomials: @polyvar
    import ProjectiveVectors: PVector
    import StaticArrays: SVector, @SVector
    import Test: @test

    const FP = FixedPolynomials
    const MP = MultivariatePolynomials
    const SP = StaticPolynomials

    export @polyvar

    include("utilities.jl")
    include("parallel.jl")
    include("affine_patches.jl")

    include("systems_and_homotopies.jl")
    include("input.jl")
    include("problems.jl")
    include("totaldegree.jl")
    include("predictors.jl")
    include("correctors.jl")

    include("path_tracking.jl")
    include("endgaming.jl")

    include("solving.jl")
    include("solve.jl")
    include("monodromy.jl")

    include("path_solving.jl")

    import LinearAlgebra: issuccess
    export issuccess
end
