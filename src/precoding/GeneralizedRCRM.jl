immutable GeneralizedRCRMState
    U_optim::Array{Matrix{Complex128},1} # normalized U
    W::Array{Hermitian{Complex128},1}    # MMSE weight, for rate calculation
    Y_ext::Array{Matrix{Float64},1}      # reweighting weights
    V_optim::Array{Matrix{Complex128},1} # normalized V
    V::Array{Matrix{Complex128},1}       # V with full power, for MMSE weight calculation
end

Papailiopoulos2011_RCRM(channel, network) =
    GeneralizedRCRM(channel, network, reweight=false, l2_reg=false)
Du2013_ReweightedRCRM(channel, network) =
    GeneralizedRCRM(channel, network, reweight=true, l2_reg=false)
Du2013_ReweightedRCRMl2Reg(channel, network) =
    GeneralizedRCRM(channel, network, reweight=true, l2_reg=true)

function GeneralizedRCRM(channel::SinglecarrierChannel, network; reweight=false, l2_reg=false)
    assignment = get_assignment(network)

    K = get_num_MSs(network); I = get_num_BSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_num_streams(network); max_d = maximum(ds)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "GeneralizedRCRM:solver" Mosek.MosekSolver(LOG=0, MAX_NUM_WARNINGS=0, NUM_THREADS=1)
    @defaultize_param! aux_params "GeneralizedRCRM:epsilon" 1e-1

    state = GeneralizedRCRMState(
        Array(Matrix{Complex128}, K),
        initial_MSE_weights(channel, Ps, sigma2s, ds, assignment, aux_params),
        [ eye(2ds[k]) for k = 1:K ],
        initial_precoders(channel, ones(I), ones(K), ds, assignment, aux_params),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    utilities = zeros(Float64, K, 1, aux_params["max_iters"])
    logdet_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = zeros(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        iter_utilities = update_MSs!(state, channel, ds, assignment, aux_params, reweight, l2_reg)
        orthogonalize!(state, channel, Ps, ds, sigma2s, assignment, aux_params)
        iters += 1

        # Results after this iteration
        utilities[:,:,iters] = iter_utilities
        push!(objective, sum(utilities[:,:,iters]))
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("GeneralizedRCRM converged.",
                    { :num_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, ds, assignment, aux_params, reweight, l2_reg)
            reweight && reweight!(state, channel, ds, assignment, aux_params)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("GeneralizedRCRM did NOT converge.",
            { :num_iters => iters, :final_objective => objective[end],
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

function update_MSs!(state::GeneralizedRCRMState, channel::SinglecarrierChannel,
    ds, assignment, aux_params, reweight, l2_reg)

    utilities = Array(Float64, channel.K)

    # Solve independently over MSs
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        # Define optimization variable
        Us = Convex.Variable(2*channel.Ns[k], ds[k]) # double first dimension for real and imag parts
        Us_r = Us[1:channel.Ns[k],:]; Us_i = Us[channel.Ns[k]+1:end,:]

        # Define user objective through auxiliary optimization variable Q_ext
        constraints = Convex.Constraint[]
        L = sum(ds) - ds[k]
        Q_ext = Convex.Variable(2ds[k], 2L); real_offset = 0; imag_offset = L
        for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
            F = channel.H[k,j]*state.V_optim[l]; F_r = real(F); F_i = imag(F)
            J_r = Us_r'*F_r + Us_i'*F_i
            J_i = Us_r'*F_i - Us_i'*F_r

            push!(constraints, Q_ext[1:ds[k],     (1 + real_offset):(ds[l] + real_offset)] ==  J_r)
            push!(constraints, Q_ext[1:ds[k],     (1 + imag_offset):(ds[l] + imag_offset)] == -J_i)
            push!(constraints, Q_ext[ds[k]+1:end, (1 + real_offset):(ds[l] + real_offset)] ==  J_i)
            push!(constraints, Q_ext[ds[k]+1:end, (1 + imag_offset):(ds[l] + imag_offset)] ==  J_r)
            real_offset += ds[l]; imag_offset += ds[l]
        end; end
        if reweight
            objective = 0.5*Convex.nuclear_norm(state.Y_ext[k]*Q_ext)
        else
            objective = 0.5*Convex.nuclear_norm(Q_ext)
        end
        if l2_reg
            objective += (1/aux_params["rho"])*Convex.sum_squares(Q_ext)
        end

        # Effective channel and user constraint
        F = channel.H[k,i]*state.V_optim[k]; F_r = real(F); F_i = imag(F)
        S_r = Us_r'*F_r + Us_i'*F_i
        S_i = Us_r'*F_i - Us_i'*F_r
        S_ext = hvcat(2, S_r, -S_i, S_i, S_r)
        push!(constraints, Convex.lambda_min(S_ext) >= aux_params["GeneralizedRCRM:epsilon"])

        # Solve local approximated RCRM optimization problem
        problem = Convex.minimize(objective, constraints)
        Convex.solve!(problem, aux_params["GeneralizedRCRM:solver"])
        if problem.status == :Optimal
            state.U_optim[k] = Convex.evaluate(Us_r) + im*Convex.evaluate(Us_i)
        else
            Lumberjack.error("Optimization problem in update_BSs!")
        end

        utilities[k] = problem.optval
    end; end

    return utilities
end

function update_BSs!(state::GeneralizedRCRMState, channel::SinglecarrierChannel,
    Ps, ds, assignment, aux_params, reweight, l2_reg)

    # Define optimization variables
    Vs = Array(Convex.Variable, channel.K)
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        Vs[k] = Convex.Variable(2*channel.Ms[i], ds[k]) # double first dimension for real and imag parts
    end; end

    # Build objective and constraints
    objective = Convex.Constant(0); constraints = Convex.Constraint[]
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        # Define user objective through auxiliary optimization variable Q_ext
        L = sum(ds) - ds[k]
        Q_ext = Convex.Variable(2ds[k], 2L); real_offset = 0; imag_offset = L
        for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
            Vs_r = Vs[l][1:channel.Ms[j],:]; Vs_i = Vs[l][channel.Ms[j]+1:end,:]
            G = channel.H[k,j]'*state.U_optim[k]; G_r = real(G); G_i = imag(G)
            J_r = G_r'*Vs_r + G_i'*Vs_i
            J_i = -G_i'*Vs_r + G_r'*Vs_i

            push!(constraints, Q_ext[1:ds[k],     (1 + real_offset):(ds[l] + real_offset)] ==  J_r)
            push!(constraints, Q_ext[1:ds[k],     (1 + imag_offset):(ds[l] + imag_offset)] == -J_i)
            push!(constraints, Q_ext[ds[k]+1:end, (1 + real_offset):(ds[l] + real_offset)] ==  J_i)
            push!(constraints, Q_ext[ds[k]+1:end, (1 + imag_offset):(ds[l] + imag_offset)] ==  J_r)
            real_offset += ds[l]; imag_offset += ds[l]
        end; end
        if reweight
            objective += 0.5*Convex.nuclear_norm(state.Y_ext[k]*Q_ext)
        else
            objective += 0.5*Convex.nuclear_norm(Q_ext)
        end
        if l2_reg
            objective += (1/aux_params["rho"])*Convex.sum_squares(Q_ext)
        end

        # Effective channel and user constraint
        Vs_r = Vs[k][1:channel.Ms[i],:]; Vs_i = Vs[k][channel.Ms[i]+1:end,:]
        G = channel.H[k,i]'*state.U_optim[k]; G_r = real(G); G_i = imag(G)
        S_r = G_r'*Vs_r + G_i'*Vs_i
        S_i = -G_i'*Vs_r + G_r'*Vs_i
        S_ext = hvcat(2, S_r, -S_i, S_i, S_r)
        push!(constraints, Convex.lambda_min(S_ext) >= aux_params["GeneralizedRCRM:epsilon"])
    end; end

    # Solve global approximated RCRM optimization problem
    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem, aux_params["GeneralizedRCRM:solver"])
    if problem.status == :Optimal
        for i = 1:channel.I; for k in served_MS_ids(i, assignment)
            Vs_r = Vs[k][1:channel.Ms[i],:]; Vs_i = Vs[k][channel.Ms[i]+1:end,:]
            state.V_optim[k] = Convex.evaluate(Vs_r) + im*Convex.evaluate(Vs_i)
        end; end
    else
        Lumberjack.error("Optimization problem in update_BSs!")
    end
end

function orthogonalize!(state::GeneralizedRCRMState,
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
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        # Received covariance and interference subspace basis
        Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
        for j = 1:channel.I; for l in served_MS_ids(j, assignment)
            #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
            Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)
        end; end

        # MMSE weight
        F = channel.H[k,i]*state.V[k]
        state.W[k] = Hermitian(inv(eye(ds[k]) - F'*(Phi\F)))
    end; end
end

function reweight!(state::GeneralizedRCRMState, channel::SinglecarrierChannel, ds, assignment, aux_params)
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        J = Array(Complex128, ds[k], sum(ds) - ds[k]); offset = 0
        for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
            Jt = state.U_optim[k]'*channel.H[k,j]*state.V_optim[l]
            J[:, (1 + offset):(ds[l] + offset)] = Jt
            offset += ds[l]
        end; end
        s = svdfact(J, thin=true)
        Y = s.U*Diagonal(1./(s.S + aux_params["delta"]))*s.U'
        Y_r = real(Y); Y_i = imag(Y)

        state.Y_ext[k] = hvcat(2, Y_r, -Y_i, Y_i, Y_r)
    end; end
end
