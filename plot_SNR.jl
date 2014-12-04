#!/usr/bin/env julia

##########################################################################
# plot_SNR.jl
#
# Plots SNR curves.
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
include("src/InterferenceRankRegularizedWSR.jl")
using InterferenceRankRegularizedWSR, CoordinatedPrecoding

plot_params = [
    "name_suffix" => "",

    "figsize" => (8,5),

    "objectives" => [
        "sumrate" => (r -> sum(r, 4:5), [ "xlabel" => "SNR [dB]", "ylabel" => "Sum rate [bits/s/Hz]" ]),
        "minrate" => (r -> minimum(sum(r, 5), 4), [ "xlabel" => "SNR [dB]", "ylabel" => "Min rate [bits/s/Hz]", ]),
    ],

    "precoding_methods" => {
        "LogDetHeuristic" => [
            ("logdet_rates", [ "key" => "g-", "legend" => "LogDetHeuristic" ]),
            ("utilities", [ "key" => "g:", "legend" => "LogDetHeuristic" ]),
        ],

        "NuclearNormHeuristic" => [
            ("logdet_rates", [ "key" => "y-", "legend" => "NuclearNormHeuristic" ]),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates", [ "key" => "b-", "legend" => "WMMSE" ]),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", [ "key" => "r-", "legend" => "MaxSINR" ]),
        ],

        "Eigenprecoding" => {
            ("intercell_tdma_logdet_rates", [ "key" => "c-", "legend" => "TDMA" ]),
            ("intracell_tdma_logdet_rates", [ "key" => "c-.", "legend" => "Intracell TDMA" ]),
            ("uncoord_logdet_rates", [ "key" => "k-", "legend" => "Uncoordinated" ]),
        },
    },
]

##########################################################################
# Plot it
plot_SNR(
    data["results"],
    data["simulation_params"],
    plot_params)
