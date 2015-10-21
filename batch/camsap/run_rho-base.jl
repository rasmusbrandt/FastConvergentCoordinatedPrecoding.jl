#!/usr/bin/env julia

##########################################################################
# run_rho.jl
#
# Performance vs. rho
##########################################################################

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using JLD, Compat

##########################################################################
# Interference channel
simulation_params = @compat Dict(
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 3,
    "P_dBm" => 30.,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,

        Shi2011_WMMSE,
        Du2013_ReweightedRCRM,
        Eigenprecoding
    ],
    "precoding_methods_nosweep" => [
        Shi2011_WMMSE,
        Du2013_ReweightedRCRM,
        Eigenprecoding
    ],
    "aux_precoding_params" => Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 4,

        "delta" => 1.,
    ),
    "independent_variable" => ((n, v) -> set_aux_precoding_param!(n, v, "rho"), logspace(-1, 2, 25)),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 2, 4]),
    ]
)
