#!/usr/bin/env julia

##########################################################################
# plot_convergence.jl
#
# Plots convergence curves.
##########################################################################

require("../../src/InterferenceRankRegularizedWSR.jl")
using InterferenceRankRegularizedWSR, CoordinatedPrecoding
using LaTeXStrings

##########################################################################
# Load data
using HDF5, JLD
data = load("convergence.jld")

##########################################################################
# Perform post processing
postprocess_params = [
    "objective" => :sumrate,
    "methods" => [
        "LogDetHeuristic" => [
            ("logdet_rates",),
            ("utilities",),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates",),
        ],

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates",),
        ],
    ],
]
results, results_mean, results_var = postprocess_convergence(data["raw_results"], data["simulation_params"], postprocess_params)

##########################################################################
# Build figure (see http://matplotlib.org/users/customizing.html for options)
PyPlot.rc("lines", linewidth=1., markersize=6)
PyPlot.rc("font", size=8, family="serif", serif="Computer Modern Sans Serif")
PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble="\\usepackage{amsmath}")
PyPlot.rc("axes", linewidth=0.5, labelsize=8)
PyPlot.rc("xtick", labelsize=8)
PyPlot.rc("ytick", labelsize=8)
PyPlot.rc("legend", fancybox=true, fontsize=8)
PyPlot.rc("figure", figsize=(3.50,2.16), dpi=125)

fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

xvals = (1:data["simulation_params"]["aux_precoding_params"]["max_iters"]) - 1 # minus one due to the definition of 'over-the-air iterations'

ax[:plot](xvals, results_mean["LogDetHeuristic"]["logdet_rates"][:,1], color="g", linestyle="-", label="TurboCP (rate)")
ax[:plot](xvals, results_mean["LogDetHeuristic"]["utilities"][:,1], color="g", linestyle="", marker=".", label="TurboCP (objective)")
ax[:plot](xvals, results_mean["Shi2011_WMMSE"]["logdet_rates"][:,1], color="b", linestyle="-", label="WMMSE")
ax[:plot](xvals, results_mean["Gomadam2008_MaxSINR"]["logdet_rates"][:,1], color="r", linestyle="-", label="MaxSINR")

ax[:set_ylim](0, 30)

ax[:set_xlabel](L"Number of over-the-air iterations $L_\text{OTA}$")
ax[:set_ylabel]("Sum rate [bits/s/Hz]")

legend = ax[:legend](loc="lower right")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

##########################################################################
# Write file
fig[:savefig]("convergence.pdf")
