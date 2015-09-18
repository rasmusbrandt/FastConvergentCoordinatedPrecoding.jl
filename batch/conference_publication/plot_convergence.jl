#!/usr/bin/env julia

##########################################################################
# plot_convergence.jl
#
# Plots convergence curves.
##########################################################################

require("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using LaTeXStrings

##########################################################################
# Load data
using JLD, Compat
data = load("convergence.jld")

##########################################################################
# Perform post processing
postprocess_params = @Compat.Dict(
    "objective" => :sum,
    "methods" => @Compat.Dict(
        "LogDetHeuristic" => [
            ("logdet_rates",),
            ("utilities",),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates",),
        ],

        "Du2013_ReweightedRCRM" => [
            ("logdet_rates",),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates",),
            ("uncoord_logdet_rates",),
        ],
    ),
)
results, results_mean, results_var = postprocess_precoding_convergence(data["raw_results"], data["simulation_params"], postprocess_params)

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

ax[:plot](xvals, results_mean["LogDetHeuristic"]["logdet_rates"][:,1], color="g", linestyle="-", label="FastDCP (rate)")
ax[:plot](xvals, results_mean["LogDetHeuristic"]["utilities"][:,1], color="g", linestyle="", marker=".", label="FastDCP (objective)")
ax[:plot](xvals, results_mean["Shi2011_WMMSE"]["logdet_rates"][:,1], color="b", linestyle="-", label="WMMSE [7]")
ax[:plot](xvals, results_mean["Du2013_ReweightedRCRM"]["logdet_rates"][:,1], color="r", linestyle="-", label="Reweighted RCRM [11]")
ax[:plot](xvals, results_mean["Eigenprecoding"]["intercell_tdma_logdet_rates"][:,1], color="c", linestyle="-", label="TDMA")
ax[:plot](xvals, results_mean["Eigenprecoding"]["uncoord_logdet_rates"][:,1], color="k", linestyle="-", label="Uncoord. transm.")

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
fig[:savefig]("convergence.eps")
