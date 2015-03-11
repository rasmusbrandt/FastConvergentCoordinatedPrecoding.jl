#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance as a function of transmit power.
##########################################################################

include("src/DoFRegularizedWSR.jl")
using DoFRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(973472333)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "SNR_$(start_time)-ic",
    "I" => 3, "Kc" => 1, "N" => 3, "M" => 3,
    "d" => 2,
    "Ndrops" => 10, "Nsim" => 1,
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
    "independent_variable" => (set_transmit_powers_dBm!, 0:3:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
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

##########################################################################
# Interfering broadcast channel channel
simulation_params = [
    "simulation_name" => "SNR_$(start_time)-ibc",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
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

        "rho" => 1/30,
        "delta" => 1.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, 0:3:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
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

##########################################################################
# Triangular3Site network
simulation_params = [
    "simulation_name" => "SNR_$(start_time)-triangular3site",
    "Kc" => 2, "N" => 2, "M" => 4,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
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

        "rho" => 1/30,
        "delta" => 1.,
    ],
    "independent_variable" => (set_transmit_powers_dBm!, 0:3:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
    ]
]
network =
    setup_triangular3site_network(
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        no_streams=simulation_params["d"])
raw_results = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
