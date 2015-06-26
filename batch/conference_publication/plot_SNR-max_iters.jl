#!/usr/bin/env julia

##########################################################################
# plot_SNR-max_iters.jl
#
# Plots SNR curves with varying number of max_iters.
##########################################################################

require("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using LaTeXStrings

##########################################################################
# Load data
using HDF5, JLD
data = load("SNR-max_iters.jld")

##########################################################################
# Perform post processing
postprocess_params = [
    "objective" => :sumrate,
    "methods" => [
        "LogDetHeuristic" => [
            ("logdet_rates",),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates",),
        ],

        "Du2013_ReweightedRCRM" => [
            ("logdet_rates",),
        ],
    ],
]
results, results_mean, results_var = postprocess(data["raw_results"], data["simulation_params"], postprocess_params)

##########################################################################
# Build figure
PyPlot.rc("lines", linewidth=1., markersize=6)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=8)
PyPlot.rc("xtick", labelsize=8)
PyPlot.rc("ytick", labelsize=8)
PyPlot.rc("legend", fancybox=true, fontsize=6)
PyPlot.rc("figure", figsize=(3.50,2.16), dpi=125)

fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

xvals = data["simulation_params"]["independent_variable"][2]

ax[:plot](xvals, results_mean["LogDetHeuristic"]["logdet_rates"][:,3], color="g", linestyle="-", label=L"TurboCP ($L_\text{OTA} = 3$)")
ax[:plot](xvals, results_mean["LogDetHeuristic"]["logdet_rates"][:,2], color="g", linestyle="--", label=L"TurboCP ($L_\text{OTA} = 2$)")
ax[:plot](xvals, results_mean["LogDetHeuristic"]["logdet_rates"][:,1], color="g", linestyle=":", label=L"TurboCP ($L_\text{OTA} = 1$)")
ax[:plot](xvals, results_mean["Shi2011_WMMSE"]["logdet_rates"][:,3], color="b", linestyle="-", label=L"WMMSE ($L_\text{OTA} = 3$)")
ax[:plot](xvals, results_mean["Shi2011_WMMSE"]["logdet_rates"][:,2], color="b", linestyle="--", label=L"WMMSE ($L_\text{OTA} = 2$)")
ax[:plot](xvals, results_mean["Shi2011_WMMSE"]["logdet_rates"][:,1], color="b", linestyle=":", label=L"WMMSE ($L_\text{OTA} = 1$)")
ax[:plot](xvals, results_mean["Du2013_ReweightedRCRM"]["logdet_rates"][:,3], color="r", linestyle="-", label=L"Reweighted RCRM ($L_\text{OTA} = 3$)")
ax[:plot](xvals, results_mean["Du2013_ReweightedRCRM"]["logdet_rates"][:,2], color="r", linestyle="--", label=L"Reweighted RCRM ($L_\text{OTA} = 2$)")
ax[:plot](xvals, results_mean["Du2013_ReweightedRCRM"]["logdet_rates"][:,1], color="r", linestyle=":", label=L"Reweighted RCRM ($L_\text{OTA} = 1$)")

ax[:set_ylim](0, 35)

ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Sum rate [bits/s/Hz]")

legend = ax[:legend](loc="upper left")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

##########################################################################
# Write file
fig[:savefig]("SNR-max_iters.eps")
