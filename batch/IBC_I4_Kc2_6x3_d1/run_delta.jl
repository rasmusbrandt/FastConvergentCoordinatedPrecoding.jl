#!/usr/bin/env julia

##########################################################################
# run_delta.jl
#
# Performance vs. delta
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
    "simulation_name" => "delta_$(start_time)",
    "I" => 4, "Kc" => 2, "N" => 3, "M" => 6,
    "P_dBm" => 30.,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Du2013_ReweightedRCRM,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 3,
        "turbo_iters" => 5,

        "rho" => 10.,
    ],
    "independent_variable" => ((n, v) -> set_aux_precoding_param!(n, v, "delta"), logspace(-1, 3, 10)),
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results, _ = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)