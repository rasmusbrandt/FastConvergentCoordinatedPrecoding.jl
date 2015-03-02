#!/usr/bin/env julia

##########################################################################
# run_convergence.jl
#
# Convergence as a function of number of outer and turbo iterations.
##########################################################################

include("../../src/DoFRegularizedWSR.jl")
using DoFRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(867123)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "convergence_$(start_time)",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => 18.2,
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
        "max_iters" => 20,

        "rho" => 1.,
        "delta" => 1.,
    ],
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 10]),
    ]
]
network =
    setup_triangular3site_network(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
