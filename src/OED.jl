"""
OED.jl - Optimal Experimental Design for Inverse Problems with State Constraints

Codes used to generate the results presented in

@article{ruthotto2017optimal,
  title={Optimal Experimental Design for Inverse Problems with State Constraints},
  author={Ruthotto, Lars and Chung, Julianne and Chung, Matthias},
  journal={SIAM Journal on Scientific Computing},
  year={2018}
}

Highlights:
- We solve the OED problem involving constraints on the states
                         (bound and linear equality constraints in inversion)
- This code parallelizes ove the training data using jInv's parallelization.
"""
module OED

	using jInv.ForwardShare
	using jInv.Utils
	using jInv.Mesh
	using jInv.InverseSolve
	using LinearOperators
	using MAT

	# descriptions of OED problems
	include("AbstractOEDParam.jl")
	include("oedTypeA.jl")
	include("oedTypeB.jl")
	include("unconstrainedOED.jl")
	include("eqconstrainedOED.jl")
	include("ineqconstrainedOED.jl")
	include("sensMat.jl")

	# optimization routines
	include("ipqp.jl")
	include("regFun.jl")
	include("models.jl")

	# codes to set up OED problem B
	include("rotation2D.jl")
	include("linearInter.jl")
	include("getLinearInterMatrix.jl")
	include("getLinearInterMatrixTranspose.jl")

	# main drivers
	include("runOEDTypeA.jl")
	include("runOEDTypeB.jl")
end
