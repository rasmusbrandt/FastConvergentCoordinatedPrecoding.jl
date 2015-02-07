#!/usr/bin/env julia

##########################################################################
# run_delta.jl
#
# Performance vs. delta
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
    "name" => "delta_$(start_time)",
    "I" => 2, "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => -9.8,
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
        "initial_precoders" => "dft",
        "stop_crit" => 0.,
        "max_iters" => 3,
        "turbo_iters" => 10,

        "rho" => 3e-2,
    ],
    "independent_variable" => ((n, v) -> set_aux_precoding_param!(n, v, "delta"), logspace(-1, 3, 100)),
]
network =
    setup_itu_r_inh_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["name"]) results")
save("$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)