#!/usr/bin/env julia

##########################################################################
# timing.jl
#
# Compares runtimes of different methods.
##########################################################################

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using Compat

##########################################################################
# General settings
srand(8071232234)
start_time = Libc.strftime("%Y%m%dT%H%M%S", time())

##########################################################################
# Performance test
simulation_params = @compat Dict(
    "I" => 6, "Kc" => 1, "N" => 2, "M" => 3,
    "P_dBm" => 30.,
    "d" => 1,
    "Ntest" => 100,
    "precoding_methods" => [
        LogDetHeuristic,
        # NuclearNormHeuristic,

        Papailiopoulos2011_RCRM,
        Du2013_ReweightedRCRM,
        Du2013_ReweightedRCRMl2Reg,
        Shi2011_WMMSE,
        Gomadam2008_MaxSINR,
        Eigenprecoding
    ],
    "aux_precoding_params" => Dict(
        "initial_precoders" => "eigendirection",
        "stop_crit" => 0.,
        "max_iters" => 4,
        "turbo_iters" => 4,

        "rho" => 100.,
        "delta" => 1.,
    ),
)
network =
    setup_interfering_broadcast_channel(simulation_params["I"],
        simulation_params["Kc"], simulation_params["N"], simulation_params["M"],
        transmit_power=10^(simulation_params["P_dBm"]/10),
        num_streams=simulation_params["d"])

timing(network, simulation_params)
