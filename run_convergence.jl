#!/usr/bin/env julia

##########################################################################
# run_convergence.jl
#
# Convergence as a function of number of iterations.
##########################################################################

include("src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using HDF5, JLD

##########################################################################
# General settings
srand(867123)
start_time = strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "convergence_$(start_time)-ic",
    "I" => 3, "Kc" => 1, "N" => 3, "M" => 3,
    "P_dBm" => 30.,
    "d" => 2,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Papailiopoulos2011_RCRM,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 100,

        "rho" => 10.,
        "delta" => 1.,
    ],
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 10]),
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate_precoding_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
##########################################################################
# Interfering broadcast channel
simulation_params = [
    "simulation_name" => "convergence_$(start_time)-ibc",
    "I" => 3, "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => 30.,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Papailiopoulos2011_RCRM,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 20,

        "rho" => 10.,
        "delta" => 1.,
    ],
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 10]),
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate_precoding_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)

##########################################################################
# Triangular3Site network
simulation_params = [
    "simulation_name" => "convergence_$(start_time)-triangular3site",
    "Kc" => 2, "N" => 2, "M" => 4,
    "P_dBm" => 30.,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Papailiopoulos2011_RCRM,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => [
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 20,

        "rho" => 10.,
        "delta" => 1.,
    ],
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 10]),
    ]
]
network =
    setup_triangular3site_network(
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        no_streams=simulation_params["d"])
raw_results = simulate_precoding_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
