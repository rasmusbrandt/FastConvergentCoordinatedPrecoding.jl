#!/usr/bin/env julia

##########################################################################
# plot_convergence.jl
#
# Plots convergence curves.
##########################################################################

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using HDF5, JLD, ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "file_name"
        help = "file name with results"
        required = true
end
parsed_args = parse_args(s)
data = load(parsed_args["file_name"])

##########################################################################
# Plot parameters
include("src/IntRankRegularizedWSR.jl")
using IntRankRegularizedWSR, CoordinatedPrecoding

if haskey(data["simulation_params"], "I")
  I = data["simulation_params"]["I"]
else
  error("simulation_params does not contain I!")
end
if haskey(data["simulation_params"], "Kc")
  Kc = data["simulation_params"]["Kc"]
else
  error("simulation_params does not contain Kc!")
end
K = I*Kc

plot_params = [
    "name_suffix" => "",

    "figsize" => (8,4),

    "objectives" => [
        "sumrate" => (r -> sum(r, 3:4), [ "xlabel" => "Iterations", "ylabel" => "Sum rate [bits/s/Hz]" ]),
    ],

    "precoding_methods" => {
        "Shi2011_WMMSE" => [
            ("logdet_rates", [ "key" => "b-", "legend" => "WMMSE" ]),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", [ "key" => "r-", "legend" => "MaxSINR" ]),
        ],

        "LogDetHeuristic" => [
            ("logdet_rates", [ "key" => "g-", "legend" => "IARegularized-WMMSE" ]),
        ],

        # "Eigenprecoding" => {
        #     ("intercell_tdma_rates", [ "key" => "c-", "legend" => "TDMA" ]),
        #     ("intracell_tdma_rates", [ "key" => "c-.", "legend" => "Intracell TDMA" ]),
        #     ("uncoord_rates", [ "key" => "k-", "legend" => "Uncoordinated" ]),
        # },
    },
]

##########################################################################
# Plot it
plot_convergence(
    data["results"],
    data["simulation_params"],
    plot_params)
