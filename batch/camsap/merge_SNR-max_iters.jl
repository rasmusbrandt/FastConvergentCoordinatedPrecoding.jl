#!/usr/bin/env julia

using FastConvergentCoordinatedPrecoding, CoordinatedPrecoding
using JLD, Compat

##########################################################################
# Postprocessing parameters
postprocess_params = @compat Dict(
    "objective" => :sum,
    "methods" => Dict(
        "LogDetHeuristic" => [
            ("logdet_rates",),
        ],

        "Shi2011_WMMSE" => [
            ("logdet_rates",),
        ],

        "Du2013_ReweightedRCRM" => [
            ("logdet_rates",),
        ],
    ),
)

##########################################################################
# Load data
sim_names = vcat([ "SNR-max_iters-$n.jld" for n = 0:9 ], [ "SNR-max_iters-W$n.jld" for n = 0:9 ])

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
save("SNR-max_iters-merged.jld",
     "simulation_params", simulation_params,
     "results", results,
     "results_mean", results_mean,
     "results_var", results_var)
