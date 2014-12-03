##########################################################################
# Default settings
function check_and_defaultize_settings!(settings::Dict{ASCIIString, Any})
    if !haskey(settings, "user_priorities")
        error("Supply user_priorities.")
    end
    if !haskey(settings, "output_protocol")
        settings["output_protocol"] = 1
    end
    if !haskey(settings, "stop_crit")
        settings["stop_crit"] = 1e-3
    end
    if !haskey(settings, "max_iters")
        settings["max_iters"] = 500
    end
    if !haskey(settings, "initial_precoders")
        settings["initial_precoders"] = "dft"
    end
    if settings["output_protocol"] != 1 && settings["output_protocol"] != 2
        error("Unknown output protocol")
    end

    if !haskey(settings, "rho")
        settings["rho"] = 1.
    end
    if !haskey(settings, "delta")
        settings["delta"] = 1e-2
    end
end

##########################################################################
# Utility and rate calculation
HeuristicState = Union(LogDetHeuristicState, NuclearNormHeuristicState)

function calculate_logdet_rates(state::HeuristicState, settings)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    logdet_rates_objective = 0.
    logdet_rates = Array(Float64, K, max_d)

    for k = 1:K
        # W is p.d., so we should only get real eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less than 1, so we need to handle that to not get negative
        # rates.
        r = log2(max(1, real(eigvals(state.W[k]))))
        logdet_rates_objective += settings["user_priorities"][k]*sum(r)

        if ds[k] < max_d
            logdet_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            logdet_rates[k,:] = r
        end
    end

    return logdet_rates, logdet_rates_objective
end

function calculate_MMSE_rates(state::HeuristicState, settings)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    MMSE_rates_objective = 0.
    MMSE_rates = Array(Float64, K, max_d)

    for k = 1:K
        # Invert W to get MMSE matrix and pick out diagonal elements, to obtain
        # MSE performance of MMSE receiver. Then take log2 of the reciprocal
        # to get the rates for the streams.
        # N.B. it is a little wasteful to take the inverse here, but it makes
        # for cleaner code in the algorithms. Also, these matrices are typically
        # small (like 2-by-2 or 3-by-3), so the complexity isn't too bad.
        E = state.W[k]\eye(state.W[k])
        r = log2(max(1, real(1./diag(E))))
        MMSE_rates_objective += settings["user_priorities"][k]*sum(r)

        if ds[k] < max_d
            MMSE_rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            MMSE_rates[k,:] = r
        end
    end

    return MMSE_rates, MMSE_rates_objective
end

function calculate_allocated_power(state::HeuristicState)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    allocated_power = Array(Float64, K, max_d)

    for k = 1:K
        p = [ vecnorm(state.V[k][:,n])^2 for n = 1:ds[k] ]

        if ds[k] < max_d
            allocated_power[k,:] = cat(1, p, zeros(Float64, max_d - ds[k]))
        else
            allocated_power[k,:] = p
        end
    end

    return allocated_power
end