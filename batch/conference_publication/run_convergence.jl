#!/usr/bin/env julia

##########################################################################
# run_convergence.jl
#
# Convergence as a function of number of outer and turbo iterations.
##########################################################################

require("../../src/DoFRegularizedWSR.jl")
using DoFRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(873232123)

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "convergence",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => 20.,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,

        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
    ],
    "aux_precoding_params" => [
        "initial_receivers" => "eigendirection", "initial_precoders" => "dft",
        "stop_crit" => 0.,
        "max_iters" => 9,
        "turbo_iters" => 2,

        "rho" => 1e-1,
        "delta" => 1.,
    ],
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
