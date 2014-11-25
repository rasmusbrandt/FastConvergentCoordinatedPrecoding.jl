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
include("src/IntRankRegularizedWSR.jl")
using IntRankRegularizedWSR, CoordinatedPrecoding

plot_params = [
    "name_suffix" => "",

    "figsize" => (8,5),

    "objectives" => [
        "sumrate" => (r -> sum(r, 4:5), [ "xlabel" => "SNR [dB]", "ylabel" => "Sum rate [bits/s/Hz]" ]),
        "minrate" => (r -> minimum(sum(r, 5), 4), [ "xlabel" => "SNR [dB]", "ylabel" => "Min rate [bits/s/Hz]", ]),
    ],

    "precoding_methods" => {
        "LogDetHeuristic" => {
            ("MMSE_rates", [ "key" => "g-", "legend" => "LogDetHeuristic bound" ]),
        },
        
        "Shi2011_WMMSE" => {
            ("MMSE_rates", [ "key" => "b-", "legend" => "Shi2011_WMMSE bound" ]),
        },

        "Gomadam2008_MaxSINR" => {
            ("MMSE_rates", [ "key" => "r-", "legend" => "Gomadam2008_MaxSINR bound" ]),
        },

        "Eigenprecoding" => {
            ("intercell_tdma_MMSE_rates", [ "key" => "c-", "legend" => "TDMA bound" ]),
            ("intracell_tdma_MMSE_rates", [ "key" => "c-.", "legend" => "Intracell TDMA bound" ]),
            ("uncoord_MMSE_rates", [ "key" => "k-", "legend" => "Uncoordinated bound" ]),
        },
    },
]

##########################################################################
# Plot it
plot_SNR(
    data["results"],
    data["simulation_params"],
    plot_params)
