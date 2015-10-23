immutable NuclearNormHeuristicState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # inverse of MMSE matrix
    Y::Array{Hermitian{Complex128},1} # inverse of MSE matrix
    V::Array{Matrix{Complex128},1}
end

function NuclearNormHeuristic(channel::SinglecarrierChannel, network)
    assignment = get_assignment(network)

    K = get_num_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_num_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network); alphas_diagonal = Diagonal(alphas)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "NuclearNormHeuristic:solver" Mosek.MosekSolver(LOG=0, MAX_NUM_WARNINGS=0)

    state = NuclearNormHeuristicState(
        Array(Matrix{Complex128}, K),
        initial_MSE_weights(channel, Ps, sigma2s, ds, assignment, aux_params),
        initial_MSE_weights(channel, Ps, sigma2s, ds, assignment, aux_params),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    utilities = zeros(Float64, K, max_d, aux_params["max_iters"])
    logdet_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = zeros(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = zeros(Float64, K, max_d, aux_params["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        # Mosek fails hard when it doesn't find a solution, so we just wrap
        # the whole algorithm in this try/catch block.
        try
            update_MSs!(state, channel, ds, sigma2s, assignment, aux_params)
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
                    Lumberjack.debug("NuclearNormHeuristic converged.",
                        @compat Dict(
                            :num_iters => iters, :final_objective => objective[end],
                            :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                            :max_iters => aux_params["max_iters"]))
                    break
                end
            end

            # Begin next iteration, unless the loop will end
            if iters < aux_params["max_iters"]
                update_BSs!(state, channel, Ps, ds, assignment, aux_params)
            end
        catch e
            Lumberjack.warn("NuclearNormHeuristic bailing due to: $e.")
            break
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("NuclearNormHeuristic did NOT converge.",
            @compat Dict(
                :num_iters => iters, :final_objective => objective[end],
                :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                :max_iters => aux_params["max_iters"]))
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

function update_MSs!(state::NuclearNormHeuristicState,
    channel::SinglecarrierChannel, ds, sigma2s, assignment, aux_params)

    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            # Optimization variable
            Us = Convex.Variable(2*channel.Ns[k], ds[k])

            # Received covariance and interference subspace basis
            constraints = Convex.Constraint[]
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            Q = Convex.Variable(ds[k], 2*(sum(ds) - ds[k])); q_offset = 0
            for j = 1:channel.I
                for l in served_MS_ids(j, assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)

                    if l != k
                        push!(constraints, Q[:, (1 + q_offset):(q_offset + 2*ds[l])] == Us'*cmat(channel.H[k,j]*state.V[l]))
                        q_offset += ds[l]
                    end
                end
            end

            # Effective desired channel
            F = channel.H[k,i]*state.V[k]

            # Receive filter
            Us = Convex.Variable(2*channel.Ns[k], ds[k])
            MSE_term = Convex.sumsquares(cmat(sqrtm(full(Phi)))*Us) - 2*trace(Us'*cvec(F))
            IntfNN = Convex.nuclearnorm(Q)
            objective = MSE_term + aux_params["rho"]*IntfNN
            problem = Convex.minimize(objective, constraints)
            Convex.solve!(problem, aux_params["NuclearNormHeuristic:solver"])
            if problem.status == :Optimal
                Ur = Us.value[1:channel.Ns[k],:]; Ui = Us.value[channel.Ns[k]+1:end,:]
                state.U[k] = Ur + im*Ui
            else
                Lumberjack.error("Optimization problem in update_MSs!")
            end

            # MSE
            E = eye(ds[k]) - state.U[k]'*F - F'*state.U[k] + state.U[k]'*Phi*state.U[k]
            state.Y[k] = Hermitian(inv(E))

            # MSE weight for rate calculation (w/ MMSE filter)
            state.W[k] = Hermitian(inv(eye(ds[k]) - F'*(Phi\F)))
        end
    end
end

function update_BSs!(state::NuclearNormHeuristicState,
    channel::SinglecarrierChannel, Ps, ds, assignment, aux_params)

    Vs = Array(Convex.Variable, channel.K)
    objective = Convex.Constant(0)
    constraints = Convex.Constraint[]

    # Precoder definition, power constraints and weighted MMSE objective contribution
    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        
        for j = 1:channel.I
            for l in served_MS_ids(j, assignment)
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.Y[l]*state.U[l]')*channel.H[l,i])
            end
        end
        Gamma_ext_sqrtm = cmat(sqrtm(full(Gamma)))

        # Precoders and weighted MMSE objective contribution for served users
        used_power = Convex.Constant(0)
        for k in served_MS_ids(i, assignment)
            # Effective desired channel
            G = channel.H[k,i]'*state.U[k]*state.Y[k]

            Vs[k] = Convex.Variable(2*channel.Ms[i], ds[k])
            objective += Convex.sumsquares(Gamma_ext_sqrtm*Vs[k]) - 2*trace(cvec(G)'*Vs[k])
            used_power += Convex.sumsquares(Vs[k])
        end
        push!(constraints, used_power <= Ps[i])
    end

    # Interference subspace basis nuclear norm regularization
    for i = 1:channel.I
        for k in served_MS_ids(i, assignment)
            Q = Convex.Variable(2*ds[k], sum(ds) - ds[k]); q_offset = 0

            for j = 1:channel.I
                for l in served_MS_ids_except_me(k, j, assignment)
                    push!(constraints, Q[:, (1 + q_offset):(q_offset + ds[l])] == cmat(channel.H[k,j]'*state.U[k])'*Vs[l])
                    q_offset += ds[l]
                end
            end

            user_IntfNN = Convex.nuclearnorm(Q)

            objective += aux_params["rho"]*user_IntfNN
        end
    end

    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem, aux_params["NuclearNormHeuristic:solver"])
    if problem.status == :Optimal
        for i = 1:channel.I
            for k in served_MS_ids(i, assignment)
                Vr = Vs[k].value[1:channel.Ms[i],:]; Vi = Vs[k].value[channel.Ms[i]+1:end,:]
                state.V[k] = Vr + im*Vi
            end
        end
    else
        Lumberjack.error("Optimization problem in update_BSs!")
    end
end
