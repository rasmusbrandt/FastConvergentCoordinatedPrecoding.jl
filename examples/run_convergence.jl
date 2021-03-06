#!/usr/bin/env julia

##########################################################################
# run_convergence.jl
#
# Convergence as a function of number of iterations.
##########################################################################

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using JLD, Compat

##########################################################################
# General settings
srand(867123)
start_time = Libc.strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Simulation
simulation_params = @compat Dict(
    "simulation_name" => "convergence_$(start_time)",
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 3,
    "P_dBm" => 30.,
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
        "max_iters" => 10,

        "rho" => 10.,
        "delta" => 1.,
    ),
    "aux_independent_variables" => [
        ((n, v) -> set_aux_precoding_param!(n, v, "turbo_iters"), [1, 10]),
    ]
)
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        num_streams=simulation_params["d"])
raw_results = simulate_precoding_convergence(network, simulation_params)

println("-- Saving $(simulation_params["simulation_name"]) results")
save("$(simulation_params["simulation_name"]).jld",
     "simulation_params", clean_simulation_params_for_jld(simulation_params),
     "raw_results", raw_results)
