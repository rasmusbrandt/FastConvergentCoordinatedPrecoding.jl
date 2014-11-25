##########################################################################
# IntRankRegularizedWSR
#
# Evaluation environment for the IntRankRegularizedWSR project
# https://gitr.sys.kth.se/rabr5411/IntRankRegularizedWSR.jl
##########################################################################

module IntRankRegularizedWSR

using CoordinatedPrecoding
import Lumberjack

export
    # precoding
    LogDetHeuristic, NuclearNormHeuristic

include("precoding/LogDetHeuristic.jl")
include("precoding/NuclearNormHeuristic.jl")

include("logging.jl")
include("utils.jl")

end
