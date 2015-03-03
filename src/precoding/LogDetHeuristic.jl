immutable LogDetHeuristicState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # inverse of MMSE matrix
    Y::Array{Hermitian{Complex128},1} # inverse of MSE matrix
    Z::Array{Hermitian{Complex128},1} # inverse of leakage matrix
    V::Array{Matrix{Complex128},1}
end

function LogDetHeuristic(channel::SinglecarrierChannel, network::Network)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network); alphas_diagonal = Diagonal(alphas)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "LogDetHeuristic:bisection_matrix_cond" 1e10
    @defaultize_param! aux_params "LogDetHeuristic:bisection_singular_matrix_mu_lower_bound" 1e-14
    @defaultize_param! aux_params "LogDetHeuristic:bisection_max_iters" 1e2
    @defaultize_param! aux_params "LogDetHeuristic:bisection_tolerance" 1e-3

    state = LogDetHeuristicState(
        Array(Matrix{Complex128}, K),
        unity_MSE_weights(ds),
        unity_MSE_weights(ds),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    utilities = Array(Float64, K, max_d, aux_params["max_iters"])
    logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, sigma2s, assignment, aux_params)
        iters += 1

        # Results after this iteration
        utilities[:,:,iters] = calculate_utilities(state, aux_params)
        push!(objective, sum(alphas_diagonal*utilities[:,:,iters]))
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("LogDetHeuristic converged.",
                    { :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, alphas, assignment, aux_params)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("LogDetHeuristic did NOT converge.",
            { :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] })
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["utilities"] = utilities
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["utilities"] = utilities[:,:,iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::LogDetHeuristicState, channel::SinglecarrierChannel,
    sigma2s, assignment, aux_params)

    rho = aux_params["rho"]; delta = aux_params["delta"]
    ds = [ size(state.W[k], 1) for k = 1:channel.K ]

    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            # Received interference covariance
            Psi = Hermitian(complex(zeros(channel.Ns[k], channel.Ns[k])))
            for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
                #Psi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                Base.LinAlg.BLAS.herk!(Psi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Psi.S)
            end; end

            # Effective desired channel
            effective_channel = channel.H[k,i]*state.V[k]

            # Total received signal covariance
            Phi = Hermitian(Psi + effective_channel*effective_channel' + sigma2s[k]*eye(channel.Ns[k]))

            # MSE weight for rate calculation (w/ MMSE filter)
            state.W[k] = Hermitian(inv(eye(ds[k]) - effective_channel'*(Phi\effective_channel)))

            for turbo_iters = 1:aux_params["turbo_iters"]
                # Receive filter (N.B., not MMSE filter!)
                state.U[k] = reshape((kron(transpose(full(state.Y[k])), full(Phi)) + (1/rho)*kron(transpose(full(state.Z[k])), full(Psi)))\vec(effective_channel*state.Y[k]), channel.Ns[k], ds[k])

                # MSE
                E = eye(ds[k]) - state.U[k]'*effective_channel - effective_channel'*state.U[k] + state.U[k]'*Phi*state.U[k]
                state.Y[k] = Hermitian(inv(E))

                # Leakage
                F = state.U[k]'*Psi*state.U[k]
                state.Z[k] = Hermitian(inv(delta*eye(ds[k]) + F))
            end
        end
    end
end

function update_BSs!(state::LogDetHeuristicState, channel::SinglecarrierChannel, 
    Ps, alphas, assignment, aux_params)

    for i = 1:channel.I
        served = served_MS_ids(i, assignment); Kc = length(served)

        Gamma = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            Gamma += Hermitian(alphas[l]*channel.H[l,i]'*(state.U[l]*state.Y[l]*state.U[l]')*channel.H[l,i])
        end; end

        Lambdas = Array(Hermitian{Complex128}, Kc); k_idx = 1
        for k in served
            Lambdas[k_idx] = Hermitian(complex(zeros(channel.Ms[i], channel.Ms[i])))
            for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
                Lambdas[k_idx] += Hermitian(alphas[l]*channel.H[l,i]'*(state.U[l]*state.Z[l]*state.U[l]')*channel.H[l,i])
            end; end
            k_idx += 1
        end

        # Find optimal Lagrange multiplier
        mu_star, eigens =
            optimal_mu(i, Gamma, Lambdas, state, channel, Ps, alphas, assignment, aux_params)

        # Precoders (reuse EVDs)
        k_idx = 1
        for k in served
            state.V[k] = alphas[k]*eigens[k_idx].vectors*Diagonal(1./(abs(eigens[k_idx].values) .+ mu_star))*eigens[k_idx].vectors'*channel.H[k,i]'*state.U[k]*state.Y[k]
            k_idx += 1
        end
    end
end

function optimal_mu(i, Gamma, Lambdas, state::LogDetHeuristicState,
    channel::SinglecarrierChannel, Ps, alphas, assignment, aux_params)

    rho = aux_params["rho"]
    served = served_MS_ids(i, assignment); Kc = length(served)

    # Build bisector function
    eigens = Array(Factorization{Complex128}, Kc)
    bis = Array(Float64, channel.Ms[i], Kc)
    k_idx = 1
    for k in served
        eigens[k_idx] = eigfact(Gamma + (1/rho)*Lambdas[k_idx])

        effective_channel = alphas[k]*channel.H[k,i]'*state.U[k]*state.Y[k]
        bis[:,k_idx] = abs(diag(eigens[k_idx].vectors'*(effective_channel*effective_channel')*eigens[k_idx].vectors))
        k_idx += 1
    end

    f(mu) =
        begin
            a = 0.; k_idx = 1
            for k in served
                a += sum(bis[:,k_idx]./((eigens[k_idx].values .+ mu).*(eigens[k_idx].values .+ mu)))
                k_idx += 1
            end
            return a
        end

    # mu lower bound
    mu_lower = 0.; k_idx = 1
    for k in served
        if maximum(abs(eigens[k_idx].values))/minimum(abs(eigens[k_idx].values)) > aux_params["LogDetHeuristic:bisection_matrix_cond"]
            # Matrix not invertible, resort to non-zero default
            mu_lower = aux_params["LogDetHeuristic:bisection_singular_matrix_mu_lower_bound"]
            break
        end
        k_idx += 1
    end

    if f(mu_lower) <= Ps[i]
        # No bisection needed
        return mu_lower, eigens
    else
        # mu upper bound
        mu_upper = sqrt((1/Ps[i])*maximum(bis)*channel.Ms[i]*sum(alphas.^2)) - minimum([minimum(abs(eigens[k_idx].values)) for k_idx = 1:Kc])
        if f(mu_upper) > Ps[i]
            error("Power bisection: infeasible mu upper bound.")
        end

        no_iters = 0
        while no_iters < aux_params["LogDetHeuristic:bisection_max_iters"]
            conv_crit = (Ps[i] - f(mu_upper))/Ps[i]

            if conv_crit < aux_params["LogDetHeuristic:bisection_tolerance"]
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

        if no_iters == aux_params["LogDetHeuristic:bisection_max_iters"]
            println("Power bisection: reached max iterations.")
        end

        # The upper point is always feasible, therefore we use it
        return mu_upper, eigens
    end
end

function calculate_utilities(state::LogDetHeuristicState, aux_params)
    rho = aux_params["rho"]

    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    utilities = Array(Float64, K, max_d)

    for k = 1:K
        user_util = abs(log2(eigvals(state.Y[k]))) + (1/rho)*abs(log2(eigvals(state.Z[k])))

        if ds[k] < max_d
            utilities[k,:] = cat(1, user_util, zeros(Float64, max_d - ds[k]))
        else
            utilities[k,:] = user_util
        end
    end

    return utilities
end
