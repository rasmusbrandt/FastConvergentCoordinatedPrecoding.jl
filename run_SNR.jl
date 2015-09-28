#!/usr/bin/env julia

##########################################################################
# run_SNR.jl
#
# Performance as a function of transmit power.
##########################################################################

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using JLD, Compat

##########################################################################
# General settings
srand(973472333)
start_time = Libc.strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Simulation
simulation_params = @compat Dict(
    "simulation_name" => "SNR_$(start_time)",
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 3,
    "d" => 1,
    "Ndrops" => 10, "Nsim" => 1,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Ghauch2015_Turbo,
        # Papailiopoulos2011_RCRM,
        # Du2013_ReweightedRCRM,
        # Du2013_ReweightedRCRMl2Reg,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "turbo_iters" => 4,

        "rho" => 10.,
        "delta" => 1.,
    ),
    "independent_variable" => (set_transmit_powers_dBm!, 0:10:30),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "max_iters"), [2, 3, 4]),
    ]
)
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        num_streams=simulation_params["d"])
raw_results, _ = simulate(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
