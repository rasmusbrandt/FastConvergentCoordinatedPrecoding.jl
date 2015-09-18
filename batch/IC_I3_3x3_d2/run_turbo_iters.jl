#!/usr/bin/env julia

##########################################################################
# run_turbo_iters.jl
#
# Performance vs. turbo_iters
##########################################################################

include("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using JLD

##########################################################################
# General settings
srand(867123)
start_time = Libc.strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Interference channel
simulation_params = [
    "simulation_name" => "turbo_iters_$(start_time)",
    "I" => 3, "Kc" => 1, "N" => 3, "M" => 3,
    "P_dBm" => 30.,
    "d" => 2,
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

        "rho" => 10.,
        "delta" => 1.,
    ],
    "independent_variable" => ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), 1:5),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
    ]
]
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        num_streams=simulation_params["d"])
raw_results, _ = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
