#!/usr/bin/env julia

##########################################################################
# plot_rho.jl
#
# Plots rho curves.
##########################################################################

require("../../src/InterferenceRankRegularizedWSR.jl")
using InterferenceRankRegularizedWSR, CoordinatedPrecoding
using LaTeXStrings

##########################################################################
# Load data
using HDF5, JLD
data = load("rho.jld")

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

        "Gomadam2008_MaxSINR" => [
            ("logdet_rates",),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates",),
            ("intracell_tdma_logdet_rates",),
            ("uncoord_logdet_rates",),
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
PyPlot.rc("legend", fancybox=true, fontsize=8)
PyPlot.rc("figure", figsize=(3.50,2.16), dpi=125)

fig = PyPlot.figure()
ax = fig[:add_axes]((0.11,0.15,0.95-0.11,0.95-0.15))

xvals = data["simulation_params"]["independent_variable"][2]

ax[:plot](xvals, results_mean["LogDetHeuristic"]["logdet_rates"][:,1], color="g", linestyle="-", label="TurboCP")
ax[:plot](xvals, results_mean["Eigenprecoding"]["intercell_tdma_logdet_rates"][:,1], color="c", linestyle="-", label="TDMA")
ax[:plot](xvals, results_mean["Shi2011_WMMSE"]["logdet_rates"][:,1], color="b", linestyle="-", label="WMMSE")
ax[:plot](xvals, results_mean["Gomadam2008_MaxSINR"]["logdet_rates"][:,1], color="r", linestyle="-", label="MaxSINR")
ax[:plot](xvals, results_mean["Eigenprecoding"]["uncoord_logdet_rates"][:,1], color="k", linestyle="-", label="Uncoordinated transmission")

ax[:set_xscale]("log")
ax[:set_ylim](0, 25)

ax[:set_xlabel](L"Regularization parameter $\rho$")
ax[:set_ylabel]("Sum rate [bits/s/Hz]")

legend = ax[:legend](loc="lower center")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

##########################################################################
# Write file
fig[:savefig]("rho.pdf")
