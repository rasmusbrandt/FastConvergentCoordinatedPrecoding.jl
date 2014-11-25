##########################################################################
# Global settings and consistency checks
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
end
