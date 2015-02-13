#!/usr/bin/env julia

##########################################################################
# run_turbo_iters.jl
#
# Performance vs. turbo_iters
##########################################################################

include("../../src/InterferenceRankRegularizedWSR.jl")
using InterferenceRankRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(867123)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "turbo_iters_$(start_time)",
    "I" => 2, "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => -9.8,
    "d" => 2,
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

        "rho" => 3e-2,
        "delta" => 1.,
    ],
    "independent_variable" => ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), 1:5),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
    ]
]
network =
    setup_itu_r_inh_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
