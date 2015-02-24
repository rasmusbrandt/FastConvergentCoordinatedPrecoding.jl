#!/usr/bin/env julia

##########################################################################
# plot_convergence.jl
#
# Plots convergence curves.
##########################################################################

include("../../src/InterferenceRankRegularizedWSR.jl")
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
    "plot_name" => "",

    "objective" => :sumrate,

    "figure" => [
        :figsize => (8,5),
        :dpi => 125,
    ],

    "axes" => [
        :xlabel => "Iterations",
        :ylabel => "Sum rate [bits/s/Hz]",
    ],

    "legend" => [
        :loc => "best",
        :fontsize => 8,
    ],

    "methods" => [
        "LogDetHeuristic" => [
            ("logdet_rates", [ :color => "g", :linestyle => "-", :label => "LogDetHeuristic" ]),
            ("utilities", [ :color => "g", :linestyle => "--",  :label => "LogDetHeuristic (utilities)" ]),
        ],

        "NuclearNormHeuristicLinearized" => [
            ("logdet_rates", [ :color => "y", :linestyle => "-", :label => "NuclearNormHeuristicLinearized" ]),
        ],

        "NuclearNormHeuristicMosek" => [
            ("logdet_rates", [ :color => "y", :linestyle => ":", :label => "NuclearNormHeuristicMosek" ]),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates", [ :color => "b", :linestyle => "-", :label => "WMMSE" ]),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", [ :color => "r", :linestyle => "-", :label => "MaxSINR" ]),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates", [ :color => "c", :linestyle => "-", :label => "TDMA" ]),
            ("intracell_tdma_logdet_rates", [ :color => "c", :linestyle => "-.",  :label => "Intracell TDMA" ]),
            ("uncoord_logdet_rates", [ :color => "k", :linestyle => "-", :label => "Uncoord. transm." ]),
        ],
    ]
]

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess_convergence(data["raw_results"], data["simulation_params"], plot_params)
    plot_convergence(processed_results, data["simulation_params"], plot_params)
end
