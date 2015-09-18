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
using JLD, Compat
data = load("SNR-max_iters-merged.jld")

##########################################################################
# Build figure
PyPlot.rc("lines", linewidth=1., markersize=3, markeredgewidth=0.5)
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

# Markers
l1 = ax[:plot](xvals, data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,3], color="g", markeredgecolor="g", marker="o", linestyle="-", label=L"FastDCP ($L_\text{OTA} = 3$)")
l2 = ax[:plot](xvals, data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,2], color="g", markeredgecolor="g", marker="o", linestyle="--", label=L"FastDCP ($L_\text{OTA} = 2$)")
l3 = ax[:plot](xvals, data["results_mean"]["LogDetHeuristic"]["logdet_rates"][:,1], color="g", markeredgecolor="g", marker="o", linestyle=":", label=L"FastDCP ($L_\text{OTA} = 1$)")
l4 = ax[:plot](xvals, data["results_mean"]["Shi2011_WMMSE"]["logdet_rates"][:,3], color="b", markeredgecolor="b", marker="v", linestyle="-", label=L"WMMSE ($L_\text{OTA} = 3$)")
l5 = ax[:plot](xvals, data["results_mean"]["Shi2011_WMMSE"]["logdet_rates"][:,2], color="b", markeredgecolor="b", marker="v", linestyle="--", label=L"WMMSE ($L_\text{OTA} = 2$)")
l6 = ax[:plot](xvals, data["results_mean"]["Shi2011_WMMSE"]["logdet_rates"][:,1], color="b", markeredgecolor="b", marker="v", linestyle=":", label=L"WMMSE ($L_\text{OTA} = 1$)")
l7 = ax[:plot](xvals, data["results_mean"]["Du2013_ReweightedRCRM"]["logdet_rates"][:,3], color="r", markeredgecolor="r", marker="^", linestyle="-", label=L"Reweighted RCRM ($L_\text{OTA} = 3$)")
l8 = ax[:plot](xvals, data["results_mean"]["Du2013_ReweightedRCRM"]["logdet_rates"][:,2], color="r", markeredgecolor="r", marker="^", linestyle="--", label=L"Reweighted RCRM ($L_\text{OTA} = 2$)")
l9 = ax[:plot](xvals, data["results_mean"]["Du2013_ReweightedRCRM"]["logdet_rates"][:,1], color="r", markeredgecolor="r", marker="^", linestyle=":", label=L"Reweighted RCRM ($L_\text{OTA} = 1$)")

ax[:set_ylim](0, 25)

ax[:set_xlabel]("Signal-to-noise ratio [dB]")
ax[:set_ylabel]("Sum rate [bits/s/Hz]")

legend1 = PyPlot.legend(handles = [l1, l2, l3, l4, l5, l6], loc="upper left")
legend1_frame = legend1[:get_frame]()
PyPlot.setp(legend1_frame, linewidth=0.5)
ax[:add_artist](legend1)

legend2 = ax[:legend](handles = [l7, l8, l9], loc="lower right")
legend2_frame = legend2[:get_frame]()
PyPlot.setp(legend2_frame, linewidth=0.5)

##########################################################################
# Write file
fig[:savefig]("SNR-max_iters.eps")
