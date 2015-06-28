#!/usr/bin/env julia

##########################################################################
# run_SNR-max_iters.jl
#
# Performance vs. SNR and max_iters
##########################################################################

require("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# Interference channel
simulation_params = [
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 3,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,

        Shi2011_WMMSE,
        Du2013_ReweightedRCRM,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "turbo_iters" => 4,

        "rho" => 10.,
        "delta" => 1.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -10:4:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
    ]
]
