#!/usr/bin/env julia

##########################################################################
# run_SNR-turbo_iters.jl
#
# Performance vs. SNR and turbo_iters
##########################################################################

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using JLD, Compat

##########################################################################
# General settings
srand(8123123)

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "SNR-turbo_iters",
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
        "max_iters" => 4,

        "rho" => 10.,
        "delta" => 1.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, -10:5:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 2, 4]),
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"])
raw_results, _ = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
