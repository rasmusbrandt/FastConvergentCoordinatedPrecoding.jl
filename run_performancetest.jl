#!/usr/bin/env julia

##########################################################################
# run_performancetest.jl
#
# Compares runtimes of different methods.
##########################################################################

include("src/InterferenceRankRegularizedWSR.jl")
using InterferenceRankRegularizedWSR, CoordinatedPrecoding

##########################################################################
# General settings
srand(8071232234)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Performance test
simulation_params = {
    "I" => 3, "Kc" => 1, "N" => 2, "M" => 2,
    "P_dBm" => 20.,
    "d" => 1,
    "Ntest" => 100,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristicLinearized,
        # NuclearNormHeuristicMosek,

        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "turbo_iters" => 10,

        "rho" => 1e-2,
        "delta" => 1.,
    ],
}
precoding_settings = {
    "stop_crit" => 20,

    "rho" => 1/30,
    "delta" => 1.,
}
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])

simulate_performance(network, simulation_params)
