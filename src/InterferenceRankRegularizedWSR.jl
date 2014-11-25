##########################################################################
# InterferenceRankRegularizedWSR
#
# Evaluation environment for the InterferenceRankRegularizedWSR project
# https://gitr.sys.kth.se/rabr5411/InterferenceRankRegularizedWSR.jl
##########################################################################

module InterferenceRankRegularizedWSR

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
