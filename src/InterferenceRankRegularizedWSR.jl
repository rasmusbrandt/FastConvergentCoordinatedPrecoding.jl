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
    LogDetHeuristic, NuclearNormHeuristicLinearized, NuclearNormHeuristicMosek

include("precoding/LogDetHeuristic.jl")
include("precoding/NuclearNormHeuristicLinearized.jl")
include("precoding/NuclearNormHeuristicMosek.jl")

##########################################################################
# Logging defaults
let
    console = Lumberjack._lumber_mill.timber_trucks["console"]
    Lumberjack.configure(console; mode = "warn")
    file = Lumberjack.add_truck(Lumberjack.LumberjackTruck("default.log", "info"), "default")
end

end
