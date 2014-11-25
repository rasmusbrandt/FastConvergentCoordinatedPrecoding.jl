##########################################################################
# InterferenceRankRegularizedWSR
#
# Evaluation environment for the InterferenceRankRegularizedWSR project
# https://gitr.sys.kth.se/rabr5411/InterferenceRankRegularizedWSR.jl
##########################################################################

module InterferenceRankRegularizedWSR

using CoordinatedPrecoding
import Lumberjack, Convex, ECOS, SCS

export
    # precoding
    LogDetHeuristic, NuclearNormHeuristic

include("precoding/LogDetHeuristic.jl")
include("precoding/NuclearNormHeuristic.jl")
include("precoding/precoding.jl")

include("logging.jl")

end
