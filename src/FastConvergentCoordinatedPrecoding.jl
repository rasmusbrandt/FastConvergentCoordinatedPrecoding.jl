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

##########################################################################
# Hermitian methods, previously from CoordinatedPrecoding.jl
import Base: +, -, .*, logdet, diag
+(A::Hermitian{Complex128}, B::Hermitian{Complex128}) = Hermitian(A.S + B.S)
+(B::Matrix{Float64}, A::Hermitian{Complex128}) = +(A, B)
+(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S + B
+(A::Hermitian{Complex128}, B::Matrix{Complex128}) = A.S + B

-(A::Hermitian{Complex128}, B::Matrix{Complex128}) = A.S - B
-(B::Array{Complex128, 2}, A::Hermitian{Complex128}) = -(A, B)
-(A::Hermitian{Complex128}, B::Matrix{Float64}) = A.S - B

.*(a::Float64, B::Hermitian{Complex128}) = Hermitian(a.*B.S)

logdet(A::Hermitian{Complex128}) = logdet(A.S)
diag(A::Hermitian{Complex128}) = diag(A.S)

end
