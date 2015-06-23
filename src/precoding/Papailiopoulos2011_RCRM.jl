immutable Papailiopoulos2011_RCRMState
    U_optim::Array{Matrix{Complex128},1} # normalized U
    W::Array{Hermitian{Complex128},1}
    V_optim::Array{Matrix{Complex128},1} # normalized V
    V::Array{Matrix{Complex128},1}
end

function Papailiopoulos2011_RCRM(channel::SinglecarrierChannel, network::Network)
    assignment = get_assignment(network)

    K = get_no_MSs(network); I = get_no_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network); alphas_diagonal = Diagonal(alphas)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "Papailiopoulos2011_RCRM:solver" Mosek.MosekSolver(LOG=0)
    @defaultize_param! aux_params "Papailiopoulos2011_RCRM:epsilon" 1e-1

    state = Papailiopoulos2011_RCRMState(
        Array(Matrix{Complex128}, K),
        initial_MSE_weights(channel, Ps, sigma2s, ds, assignment, aux_params),
        initial_precoders(channel, ones(I), ones(K), ds, assignment, aux_params),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    utilities = Array(Float64, K, max_d, aux_params["max_iters"])
    logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, ds, assignment, aux_params)
        orthogonalize!(state, channel, Ps, ds, sigma2s, assignment, aux_params)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        push!(objective, sum(alphas_diagonal*logdet_rates[:,:,iters]))
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("Papailiopoulos2011_RCRM converged.",
                    { :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, ds, assignment, aux_params)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("Papailiopoulos2011_RCRM did NOT converge.",
            { :no_iters => iters, :final_objective => objective[end],
              :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] })
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::Papailiopoulos2011_RCRMState,
    channel::SinglecarrierChannel, ds, assignment, aux_params)

    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            # Optimization variable
            Us = Convex.Variable(2*channel.Ns[k], ds[k])

            # Interference subspace basis and objective
            constraints = Convex.Constraint[]
            Q = Convex.Variable(ds[k], 2*(sum(ds) - ds[k])); q_offset = 0
            for j = 1:channel.I
                for l in served_MS_ids_except_me(k, j, assignment)
                    push!(constraints, Q[:, (1 + q_offset):(q_offset + 2*ds[l])] == Us'*cmat(channel.H[k,j]*state.V_optim[l]))
                    q_offset += ds[l]
                end
            end
            objective = Convex.nuclear_norm(Q)

            # Effective channel and constraint
            F = cvec(channel.H[k,i]*state.V_optim[k])
            push!(constraints, Convex.lambda_min(Us'*F) >= aux_params["Papailiopoulos2011_RCRM:epsilon"])

            # Solve local approximated RCRM optimization problem
            problem = Convex.minimize(objective, constraints)
            Convex.solve!(problem, aux_params["Papailiopoulos2011_RCRM:solver"])
            if problem.status == :Optimal
                Ur = Us.value[1:channel.Ns[k],:]; Ui = Us.value[channel.Ns[k]+1:end,:]
                state.U_optim[k] = Ur + im*Ui
            else
                Lumberjack.error("Optimization problem in update_MSs!")
            end
        end
    end
end

function orthogonalize!(state::Papailiopoulos2011_RCRMState,
    channel::SinglecarrierChannel, Ps, ds, sigma2s, assignment, aux_params)

    # Orthogonalize precoders
    for i = 1:channel.I
        served = served_MS_ids(i, assignment); Nserved = length(served)
        for k in served
            Q, _ = qr(state.V_optim[k], thin=false)
            state.V[k] = sqrt(Ps[i]/(ds[k]*Nserved))*Q[:,1:ds[k]]
        end
    end

    # Calculate MMSE weights (for rate calculation)
    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            # Received covariance and interference subspace basis
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            for j = 1:channel.I
                for l in served_MS_ids(j, assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)
                end
            end

            # MMSE weight
            F = channel.H[k,i]*state.V[k]
            state.W[k] = Hermitian(inv(eye(ds[k]) - F'*(Phi\F)))
        end
    end
end

function update_BSs!(state::Papailiopoulos2011_RCRMState,
    channel::SinglecarrierChannel, Ps, ds, assignment, aux_params)

    # Define optimization variables
    Vs = Array(Convex.Variable, channel.K)
    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            Vs[k] = Convex.Variable(2*channel.Ms[i], ds[k])
        end
    end

    # Build objective and constraints
    objective = Convex.Constant(0); constraints = Convex.Constraint[]
    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            # Interference subspace basis and user objective
            Q = Convex.Variable(2*ds[k], sum(ds) - ds[k]); q_offset = 0
            for j = 1:channel.I
                for l in served_MS_ids_except_me(k, j, assignment)
                    push!(constraints, Q[:, (1 + q_offset):(q_offset + ds[l])] == cmat(channel.H[k,j]'*state.U_optim[k])'*Vs[l])
                    q_offset += ds[l]
                end
            end
            objective += Convex.nuclear_norm(Q)

            # Effective channel and user constraint
            F = cvec(channel.H[k,i]'*state.U_optim[k])
            push!(constraints, Convex.lambda_min(F'*Vs[k]) >= aux_params["Papailiopoulos2011_RCRM:epsilon"])
        end
    end

    # Solve global approximated RCRM optimization problem
    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem, aux_params["Papailiopoulos2011_RCRM:solver"])
    if problem.status == :Optimal
        for i = 1:channel.I
            for k in served_MS_ids(i, assignment)
                Vr = Vs[k].value[1:channel.Ms[i],:]; Vi = Vs[k].value[channel.Ms[i]+1:end,:]
                state.V_optim[k] = Vr + im*Vi
            end
        end
    else
        Lumberjack.error("Optimization problem in update_BSs!")
    end
end
