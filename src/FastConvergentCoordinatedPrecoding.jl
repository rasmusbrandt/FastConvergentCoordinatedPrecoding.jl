##########################################################################
# FastConvergentCoordinatedPrecoding.jl
#
# Evaluation environment for
# R. Brandt, M. Bengtsson, "Fast-Convergent Distributed Coordinated Precoding
# for TDD Multicell MIMO Systems", IEEE Int. Workshop Computational Advances
# in Multi-Sensor Adaptive Process. (CAMSAP'15), 2015. To appear.
##########################################################################

module FastConvergentCoordinatedPrecoding

using CoordinatedPrecoding
import Lumberjack, Convex, Mosek
using Compat

export
    # proposed algorithms
    LogDetHeuristic, NuclearNormHeuristic,

    # baselines
    Ghauch2015_Turbo,
    Papailiopoulos2011_RCRM,
    Du2013_ReweightedRCRM,
    Du2013_ReweightedRCRMl2Reg

include("precoding/LogDetHeuristic.jl")
include("precoding/NuclearNormHeuristic.jl")

include("precoding/Ghauch2015_Turbo.jl")
include("precoding/GeneralizedRCRM.jl")

##########################################################################
# Logging defaults
let
    console = Lumberjack._lumber_mill.timber_trucks["console"]
    Lumberjack.configure(console; mode = "warn")
    file = Lumberjack.add_truck(Lumberjack.LumberjackTruck("default.log", "info"), "default")
end

end
