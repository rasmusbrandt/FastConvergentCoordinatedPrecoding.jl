immutable LogDetHeuristicState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # inverse of MMSE matrix
    Y::Array{Hermitian{Complex128},1} # inverse of received interference covariance
    Z::Array{Hermitian{Complex128},1} # inverse of MSE matrix
    V::Array{Matrix{Complex128},1}
end

function LogDetHeuristic(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    check_and_defaultize_settings!(settings, LogDetHeuristicState)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    state = LogDetHeuristicState(
        Array(Matrix{Complex128}, K),
        unity_MSE_weights(ds),
        unity_MSE_weights(ds),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings))
    objective = Float64[]
    logdet_rates = Array(Float64, K, maximum(ds), settings["max_iters"])
    MMSE_rates = Array(Float64, K, maximum(ds), settings["max_iters"])
    allocated_power = Array(Float64, K, maximum(ds), settings["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < settings["max_iters"]
        update_MSs!(state, channel, sigma2s, cell_assignment)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters], t = calculate_logdet_rates(state, settings)
        push!(objective, t)
        MMSE_rates[:,:,iters], _ = calculate_MMSE_rates(state, settings)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < settings["stop_crit"]
                Lumberjack.debug("LogDetHeuristic converged.",
                    { :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => settings["stop_crit"],
                      :max_iters => settings["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < settings["max_iters"]
            update_BSs!(state, channel, Ps, cell_assignment, settings)
        end
    end
    if iters == settings["max_iters"]
        Lumberjack.debug("LogDetHeuristic did NOT converge.",
            { :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => settings["stop_crit"],
              :max_iters => settings["max_iters"] })
    end

    results = Dict{ASCIIString, Any}()
    if settings["output_protocol"] == 1
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif settings["output_protocol"] == 2
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function check_and_defaultize_settings!(settings, ::Type{LogDetHeuristicState})
    # Global settings
    check_and_defaultize_settings!(settings)

    # Local settings
    if !haskey(settings, "LogDetHeuristic:bisection_Gamma_cond")
        settings["LogDetHeuristic:bisection_Gamma_cond"] = 1e10
    end
    if !haskey(settings, "LogDetHeuristic:bisection_singular_Gamma_mu_lower_bound")
        settings["LogDetHeuristic:bisection_singular_Gamma_mu_lower_bound"] = 1e-14
    end
    if !haskey(settings, "LogDetHeuristic:bisection_max_iters")
        settings["LogDetHeuristic:bisection_max_iters"] = 5e1
    end
    if !haskey(settings, "LogDetHeuristic:bisection_tolerance")
        settings["LogDetHeuristic:bisection_tolerance"] = 1e-3
    end
end

function update_MSs!(state::LogDetHeuristicState, channel::SinglecarrierChannel,
    sigma2s::Vector{Float64}, cell_assignment::CellAssignment)

    rho = 1.
    delta = 1e-2

    ds = [ size(state.W[k], 1) for k = 1:channel.K ]; dtot = sum(ds)

    for i = 1:channel.I
        for k in served_MS_ids(i, cell_assignment)
            # Received covariance
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            J = Array(Complex128, channel.Ns[k], dtot - ds[k]); J_offset = 0
            for j = 1:channel.I
                for l in served_MS_ids(j, cell_assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)

                    # Interferers
                    if l != k
                        J[:,(1 + J_offset):(J_offset+ds[k])] = channel.H[k,j]*state.V[l]
                        J_offset += ds[k]
                    end
                end
            end

            # MMSE stuff
            eff_chan = channel.H[k,i]*state.V[k]
            Ummse = Phi\eff_chan
            state.W[k] = Hermitian((eye(ds[k]) - Ummse'*eff_chan)\eye(ds[k]))

            # Receive filter
            JJh = J*J'
            state.U[k] = rho*reshape((kron(transpose(full(state.Y[k])), JJh) + rho*kron(transpose(full(state.Z[k])), full(Phi)))\vec(channel.H[k,i]*state.V[k]*state.Z[k]), channel.Ns[k], ds[k])

            # MSE
            E = eye(ds[k]) - state.U[k]'*eff_chan - eff_chan'*state.U[k] + state.U[k]'*Phi*state.U[k]
            state.Z[k] = Hermitian(inv(E))

            # Received interference
            F = state.U[k]'*JJh*state.U[k]
            state.Y[k] = Hermitian(inv(F + delta*eye(F)))
        end
    end
end

function update_BSs!(state::LogDetHeuristicState, channel::SinglecarrierChannel, 
    Ps::Vector{Float64}, cell_assignment::CellAssignment, settings)

    rho = 1.

    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        for j = 1:channel.I
            for l in served_MS_ids(j, cell_assignment)
                Gamma += rho*Hermitian(channel.H[l,i]'*(state.U[l]*state.Z[l]*state.U[l]')*channel.H[l,i])

                # Only works in IC....
                if l != i
                    Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.Y[l]*state.U[l]')*channel.H[l,i])
                end
            end
        end

        # Find optimal Lagrange multiplier
        mu_star, Gamma_eigen =
            optimal_mu(i, Gamma, state, channel, Ps, cell_assignment, settings)

        # Precoders (reuse EVD)
        for k in served_MS_ids(i, cell_assignment)
            state.V[k] = rho*Gamma_eigen.vectors*Diagonal(1./(Gamma_eigen.values .+ mu_star))*Gamma_eigen.vectors'*channel.H[k,i]'*state.U[k]*state.Z[k]
        end
    end
end

function optimal_mu(i::Int, Gamma::Hermitian{Complex128},
    state::LogDetHeuristicState, channel::SinglecarrierChannel,
    Ps::Vector{Float64}, cell_assignment::CellAssignment, settings)

    rho = 1.

    # Build bisector function
    bis_M = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
    for k in served_MS_ids(i, cell_assignment)
        #bis_M += Hermitian(channel.H[k,i]'*(state.U[k]*(state.Z[k]*state.Z[k])*state.U[k]')*channel.H[k,i])
        Base.LinAlg.BLAS.herk!(bis_M.uplo, 'N', complex(1.), rho*channel.H[k,i]'*state.U[k]*state.Z[k], complex(1.), bis_M.S)
    end
    Gamma_eigen = eigfact(Gamma)
    bis_JMJ_diag = real(diag(Gamma_eigen.vectors'*bis_M*Gamma_eigen.vectors))
    f(mu) = sum(bis_JMJ_diag./((Gamma_eigen.values .+ mu).*(Gamma_eigen.values .+ mu)))

    # mu lower bound
    if abs(maximum(Gamma_eigen.values))/abs(minimum(Gamma_eigen.values)) < settings["LogDetHeuristic:bisection_Gamma_cond"]
        # Gamma is invertible
        mu_lower = 0
    else
        mu_lower = settings["LogDetHeuristic:bisection_singular_Gamma_mu_lower_bound"]
    end

    if f(mu_lower) <= Ps[i]
        # No bisection needed
        return mu_lower, Gamma_eigen
    else
        # mu upper bound
        mu_upper = sqrt(channel.Ms[i]/Ps[i]*maximum(bis_JMJ_diag)) - abs(minimum(Gamma_eigen.values))
        if f(mu_upper) > Ps[i]
            error("Power bisection: infeasible mu upper bound.")
        end

        no_iters = 0
        while no_iters < settings["LogDetHeuristic:bisection_max_iters"]
            conv_crit = (Ps[i] - f(mu_upper))/Ps[i]

            if conv_crit < settings["LogDetHeuristic:bisection_tolerance"]
                break
            else
                mu = (1/2)*(mu_lower + mu_upper)

                if f(mu) < Ps[i]
                    # New point feasible, replace upper point
                    mu_upper = mu
                else
                    # New point not feasible, replace lower point
                    mu_lower = mu
                end
            end

            no_iters += 1
        end

        if no_iters == settings["LogDetHeuristic:bisection_max_iters"]
            println("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, Gamma_eigen
    end
end
