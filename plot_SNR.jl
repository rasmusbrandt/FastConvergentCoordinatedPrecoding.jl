#!/usr/bin/env julia

##########################################################################
# plot_SNR.jl
#
# Plots SNR curves.
##########################################################################

include("src/InterferenceRankRegularizedWSR.jl")
using InterferenceRankRegularizedWSR, CoordinatedPrecoding

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using HDF5, JLD, ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "file_names"
        help = "file names with results"
        required = true
        nargs = '+'
end
parsed_args = parse_args(s)

##########################################################################
# Plot parameters
plot_params = [
    "name_suffix" => "",

    "figsize" => (8,5),

    "objectives" => [
        "sumrate" => (r -> sum(r, 5:6), [ "xlabel" => "SNR [dB]", "ylabel" => "Sum rate [bits/s/Hz]" ]),
        "minrate" => (r -> minimum(sum(r, 6), 5), [ "xlabel" => "SNR [dB]", "ylabel" => "Min rate [bits/s/Hz]", ]),
    ],

    "precoding_methods" => {
        "LogDetHeuristic" => [
            ("logdet_rates", [ "key" => "g-", "legend" => "LogDetHeuristic" ]),
            ("utilities", [ "key" => "g:", "legend" => "LogDetHeuristic" ]),
        ],

        "NuclearNormHeuristicLinearized" => [
            ("logdet_rates", [ "key" => "y-", "legend" => "NuclearNormHeuristicLinearized" ]),
        ],

        "NuclearNormHeuristicMosek" => [
            ("logdet_rates", [ "key" => "y:", "legend" => "NuclearNormHeuristicMosek" ]),
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
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = process(data["raw_results"], data["simulation_params"], plot_params)
    plot(processed_results, data["simulation_params"], plot_params)
end
