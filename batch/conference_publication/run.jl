#!/usr/bin/env julia

##########################################################################
# run.jl
#
# Run all simulations
##########################################################################\

include("run_convergence.jl")
include("run_rho.jl")
include("run_SNR-max_iters.jl")
include("run_SNR-turbo_iters.jl")
