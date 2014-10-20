##########################################################################
# IARegularizedWMMSEEval
#
# Evaluation environment for the IARegularizedWMMSE project
# https://gitr.sys.kth.se/rabr5411/IARegularizedWMMSE.jl
##########################################################################

module IARegularizedWMMSEEval

using CoordinatedPrecoding, Base.LinAlg.BLAS
import Gurobi

export
    # algorithms.jl
    IARegularizedWMMSE

include("algorithms/IARegularizedWMMSE.jl")

end
