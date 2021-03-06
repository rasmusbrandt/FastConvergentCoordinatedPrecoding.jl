#!/usr/bin/env julia

##########################################################################
# plot_rho.jl
#
# Plots rho curves.
##########################################################################

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using LaTeXStrings

##########################################################################
# Load data
using JLD, Compat
data = load("rho-merged.jld")

##########################################################################
# Build figure
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
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

# Lines
ax[:plot](xvals, data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,3], color="g", linestyle="-")
ax[:plot](xvals, data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,2], color="g", linestyle="--")
ax[:plot](xvals, data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,1], color="g", linestyle=":")
# ax[:plot](xvals, data["results_mean"]["Eigenprecoding"]["intercell_tdma_logdet_rates"][:,1], color="c", linestyle="-")
ax[:plot](xvals, data["results_mean"]["Shi2011_WMMSE"]["logdet_rates"][:,1], color="b", linestyle="-")
ax[:plot](xvals, data["results_mean"]["Du2013_ReweightedRCRM"]["logdet_rates"][:,1], color="r", linestyle="-")
# ax[:plot](xvals, data["results_mean"]["Eigenprecoding"]["uncoord_logdet_rates"][:,1], color="k", linestyle="-")

# Markers
ax[:plot](xvals[1:3:end], data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,3][1:3:end], color="g", marker="o", markeredgecolor="g", linestyle="-", label=L"FastDCP ($L_\text{local} = 4$)")
ax[:plot](xvals[1:3:end], data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,2][1:3:end], color="g", marker="o", markeredgecolor="g", linestyle="--", label=L"FastDCP ($L_\text{local} = 2$)")
ax[:plot](xvals[1:3:end], data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,1][1:3:end], color="g", marker="o", markeredgecolor="g", linestyle=":", label=L"FastDCP ($L_\text{local} = 1$)")
# ax[:plot](xvals[1:3:end], data["results_mean"]["Eigenprecoding"]["intercell_tdma_logdet_rates"][:,1][1:3:end], color="c", marker="x", markeredgecolor="c", linestyle="-", label="TDMA")
ax[:plot](xvals[1:3:end], data["results_mean"]["Shi2011_WMMSE"]["logdet_rates"][:,1][1:3:end], color="b", marker="v", markeredgecolor="b", linestyle="-", label="WMMSE [7]")
ax[:plot](xvals[1:3:end], data["results_mean"]["Du2013_ReweightedRCRM"]["logdet_rates"][:,1][1:3:end], color="r", marker="^", markeredgecolor="r", linestyle="-", label="Reweighted RCRM [11]")
# ax[:plot](xvals[1:3:end], data["results_mean"]["Eigenprecoding"]["uncoord_logdet_rates"][:,1][1:3:end], color="k", marker=".", markeredgecolor="k", linestyle="-", label="Uncoordinated transmission")

ax[:set_xscale]("log")
ax[:set_ylim](0, 26)

ax[:set_xlabel](L"Regularization parameter $\rho$")
ax[:set_ylabel]("Sum rate [bits/s/Hz]")

legend = ax[:legend](loc="lower center")
# legend_lines = legend[:get_lines]()
legend_frame = legend[:get_frame]()
# PyPlot.setp(legend_lines, linewidth=0.5)
PyPlot.setp(legend_frame, linewidth=0.5)

##########################################################################
# Write file
fig[:savefig]("rho.eps")
