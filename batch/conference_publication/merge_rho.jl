#!/usr/bin/env julia

require("../../src/MGRegularizedWSR.jl")
using MGRegularizedWSR, CoordinatedPrecoding
using JLD, Compat

##########################################################################
# Postprocessing parameters
postprocess_params = @Compat.Dict(
    "objective" => :sum,
    "methods" => @Compat.Dict(
        "LogDetHeuristic" => [
            ("logdet_rates",),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates",),
        ],

        "Du2013_ReweightedRCRM" => [
            ("logdet_rates",),
        ],

        "Eigenprecoding" => [
            ("intercell_tdma_logdet_rates",),
            ("intracell_tdma_logdet_rates",),
            ("uncoord_logdet_rates",),
        ],
    ),
)

##########################################################################
# Load data
sim_names = vcat([ "rho-$n.jld" for n = 0:9 ], [ "rho-W$n.jld" for n = 0:9 ])

# Load first
println("Loading from $(sim_names[1])")
data = load(sim_names[1])
simulation_params = data["simulation_params"]
raw_results = data["raw_results"]

for sim_name in sim_names[2:end]
    println("Loading from $(sim_name)")
    data = load(sim_name)
    raw_results.simulation_results = cat(1, raw_results.simulation_results, data["raw_results"].simulation_results)
end
data = []

##########################################################################
# Perform post processing
results, results_mean, results_var = postprocess(raw_results, simulation_params, postprocess_params)

println("-- Saving merged results")
save("rho-merged.jld",
     "simulation_params", simulation_params,
     "results", results,
     "results_mean", results_mean,
     "results_var", results_var)
