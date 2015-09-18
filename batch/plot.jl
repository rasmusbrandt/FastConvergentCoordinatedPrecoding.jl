#!/usr/bin/env julia

##########################################################################
# plot.jl
#
# Plots curves.
##########################################################################

include("../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding

##########################################################################
# Load data
#
# Do this before loading other code, otherwise the JLD module might crash!
using JLD, ArgParse
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
plot_params = @Compat.Dict(
    "plot_name" => "",

    "objective" => :sum,

    "figure" => @Compat.Dict(
        :figsize => (8,5),
        :dpi => 125,
    ),

    "axes" => @Compat.Dict(
        :ylabel => "Sum rate [bits/s/Hz]",
    ),

    "legend" => @Compat.Dict(
        :loc => "best",
        :fontsize => 6,
    ),

    "methods" => [
        "LogDetHeuristic" => [
            ("logdet_rates", @Compat.Dict(:color => "g", :linestyle => "-", :label => "LogDetHeuristic")),
        ],

        "NuclearNormHeuristic" => [
            ("logdet_rates", @Compat.Dict(:color => "y", :linestyle => ":", :label => "NuclearNormHeuristic")),
        ],

        "Papailiopoulos2011_RCRM" => [
            ("logdet_rates", @Compat.Dict(:color => "m", :linestyle => "-", :label => "Papailiopoulos2011_RCRM")),
        ],

        "Du2013_ReweightedRCRM" => [
            ("logdet_rates", @Compat.Dict(:color => "m", :linestyle => "--", :label => "Du2013_ReweightedRCRM")),
        ],

        "Du2013_ReweightedRCRMl2Reg" => [
            ("logdet_rates", @Compat.Dict(:color => "m", :linestyle => "-.", :label => "Du2013_ReweightedRCRMl2Reg")),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates", @Compat.Dict(:color => "b", :linestyle => "-", :label => "WMMSE")),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates", @Compat.Dict(:color => "r", :linestyle => "-", :label => "MaxSINR")),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates", @Compat.Dict(:color => "c", :linestyle => "-", :label => "TDMA")),
            ("intracell_tdma_logdet_rates", @Compat.Dict(:color => "c", :linestyle => "-.",  :label => "Intracell TDMA")),
            ("uncoord_logdet_rates", @Compat.Dict(:color => "k", :linestyle => "-", :label => "Uncoord. transm.")),
        ],
    ]
)

##########################################################################
# Plot it
for file_name in parsed_args["file_names"]
    data = load(file_name)
    processed_results = postprocess(data["raw_results"], data["simulation_params"], plot_params)
    plot(processed_results, data["simulation_params"], plot_params)
end
