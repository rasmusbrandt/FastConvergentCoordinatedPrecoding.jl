#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance as a function of transmit power.
##########################################################################

include("src/IntRankRegularizedWMMSE.jl")
using IntRankRegularizedWMMSE, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

precoding_settings = [
    "stop_crit" => 0,
    "max_iters" => 200,
    "a" => "b",
]

##########################################################################
# Interference channel
simulation_params = [
    "name" => "$(start_time)-ic",
    "I" => 3, "Kc" => 1, "N" => 2, "M" => 2,
    "Ps_dBm" => 0:3:40,
    "d" => 1,
    "Ndrops" => 50, "Nsim" => 2,
    "precoding_methods" => [
        LogDetHeuristic,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ]
]
precoding_settings["user_priorities"] = ones(simulation_params["I"]*simulation_params["Kc"])
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

results = simulate_SNR(network, simulation_params, precoding_settings)

println("-- Saving $(simulation_params["name"]) results")
save("SNR_$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "precoding_settings", clean_precoding_settings_for_jld(precoding_settings),
     "results", results)
quit()
##########################################################################
# Interfering broadcast channel channel
simulation_params = [
    "name" => "$(start_time)-ibc",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "Ps_dBm" => 0:3:40,
    "d" => 1,
    "Ndrops" => 5, "Nsim" => 2,
    "precoding_methods" => [
        LogDetHeuristic,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ]
]
precoding_settings["user_priorities"] = ones(simulation_params["I"]*simulation_params["Kc"])
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])

results = simulate_SNR(network, simulation_params, precoding_settings)

println("-- Saving $(simulation_params["name"]) results")
save("SNR_$(simulation_params["name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "precoding_settings", clean_precoding_settings_for_jld(precoding_settings),
     "results", results)
