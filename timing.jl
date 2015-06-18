#!/usr/bin/env julia

##########################################################################
# timing.jl
#
# Compares runtimes of different methods.
##########################################################################

include("src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding

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
        NuclearNormHeuristicLinearized,
        NuclearNormHeuristicMosek,

        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_receivers" => "eigendirection",
        "initial_precoders" => "dft",
        "stop_crit" => 0.,
        "max_iters" => 10,
        "turbo_iters" => 10,

        "rho" => 10.,
        "delta" => 1.,
    ],
}
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])

timing(network, simulation_params)
