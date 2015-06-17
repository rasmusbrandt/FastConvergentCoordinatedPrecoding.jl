#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance vs. SNR
##########################################################################

include("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(867123)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "SNR_$(start_time)",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 3,

        "rho" => 3e-2,
        "delta" => 1.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, 0:3:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [2, 3, 4]),
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
