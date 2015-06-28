#!/usr/bin/env julia

##########################################################################
# run_rho.jl
#
# Performance vs. rho
##########################################################################

require("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Interference channel
simulation_params = [
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 3,
    "P_dBm" => 30.,
    "d" => 1,
    "Ndrops" => 100, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,

        Shi2011_WMMSE,
        Du2013_ReweightedRCRM,
        Eigenprecoding
    ],
    "dont_sweep_precoding_methods" => [
        Shi2011_WMMSE,
        Du2013_ReweightedRCRM,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 4,

        "delta" => 1.,
    ],
    "independent_variable" => ((n, v) -> set_aux_precoding_param!(n, v, "rho"), logspace(-1, 2, 20)),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 2, 4]),
    ]
]
