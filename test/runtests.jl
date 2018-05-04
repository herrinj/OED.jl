println("=== Tests for OED.jl  ===")
include("testOEDTypeA.jl")
include("testUnconstrainedOED.jl")
include("testEqConstrOED.jl")
include("testIneqConstrOED.jl")

include("testSolveKKT.jl")
include("testIPQP.jl")

include("testOEDTypeB.jl")
include("testRotation2D.jl")
include("testLinearInter.jl")
include("testGetLinearInterMatrixTranspose.jl")

include("testRunOEDTypeA.jl")
include("testRunOEDTypeB.jl")
println("=== [OED.jl] : all tests passed ===")
