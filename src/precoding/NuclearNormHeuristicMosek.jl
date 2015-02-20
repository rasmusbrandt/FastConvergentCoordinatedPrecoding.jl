immutable NuclearNormHeuristicMosekState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1}
    V::Array{Matrix{Complex128},1}
end

function NuclearNormHeuristicMosek(channel::SinglecarrierChannel, network::Network)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network); alphas_diagonal = Diagonal(alphas)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "NuclearNormHeuristicMosek:perform_regularization" true
    @defaultize_param! aux_params "NuclearNormHeuristicMosek:solver" SCS.SCSSolver(verbose=0)

    state = NuclearNormHeuristicMosekState(
        Array(Matrix{Complex128}, K),
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
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        push!(objective, sum(alphas_diagonal*logdet_rates[:,:,iters]))
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("NuclearNormHeuristic converged.",
                    { :no_iters => iters, :final_objective => objective[end],
                      :conv_crit => conv_crit, :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] })
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, assignment, aux_params)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("NuclearNormHeuristic did NOT converge.",
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

function update_MSs!(state::NuclearNormHeuristicMosekState,
    channel::SinglecarrierChannel, sigma2s::Vector{Float64},
    assignment::Assignment, aux_params)

    ds = [ size(state.W[k], 1) for k = 1:length(state.W) ]

    for i = 1:channel.I
        for k = served_MS_ids(i, assignment)
            # Received covariance and interference subspace basis
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            J_ext = Array(Float64, 2*channel.Ns[k], sum(ds) - ds[k]); j_offset = 0
            for j = 1:channel.I
                Hr = real(channel.H[k,j]); Hi = imag(channel.H[k,j])
                H_ext = hvcat((2,2), Hr, -Hi, Hi, Hr)
                for l in served_MS_ids(j, assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    Base.LinAlg.BLAS.herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)

                    if l != k
                        Vr = real(state.V[l]); Vi = imag(state.V[l])
                        V_ext = vcat(Vr, Vi)

                        J_ext[:,(1 + j_offset):(j_offset + ds[l])] = H_ext*V_ext
                        j_offset += ds[l]
                    end
                end
            end

            # Find receiver
            F = channel.H[k,i]*state.V[k]
            Fr = real(F); Fi = imag(F)
            F_ext = vcat(Fr, Fi) # only retain real part in multiplication
            Phi_r = real(full(Phi)); Phi_i = imag(full(Phi))
            Phi_ext = Symmetric(hvcat((2,2), Phi_r, -Phi_i, Phi_i, Phi_r))

            # MSE term
            Us = Convex.Variable(2*channel.Ns[k], ds[k])
            MSE = ds[k] - 2*trace(Us'*F_ext) + Convex.sum_squares(sqrtm(Phi_ext)*Us)

            if aux_params["NuclearNormHeuristicMosek:perform_regularization"]
                IntfNN_obj = Convex.nuclear_norm(Us'*J_ext)

                problem = Convex.minimize(MSE + aux_params["rho"]*IntfNN_obj)
            else
                # Solve standard MSE problem
                problem = Convex.minimize(MSE)
            end

            Convex.solve!(problem, aux_params["NuclearNormHeuristicMosek:solver"])
            if problem.status == :Optimal
                Ur = Us.value[1:channel.Ns[k],:]; Ui = Us.value[channel.Ns[k]+1:end,:]
                state.U[k] = Ur + im*Ui
            else
                error("Problem with Convex.jl in update_MSs!")
            end

            # Find MSE weight
            state.W[k] = Hermitian((eye(ds[k]) - state.U[k]'*F - F'*state.U[k] + state.U[k]'*Phi*state.U[k])\eye(ds[k]))
        end
    end
end

function update_BSs!(state::NuclearNormHeuristicMosekState,
    channel::SinglecarrierChannel, Ps::Vector{Float64},
    assignment::Assignment, aux_params)

    ds = [ size(state.W[k], 1) for k = 1:length(state.W) ]

    Vs = Array(Convex.Variable, channel.K)
    objective = Convex.Constant(0)
    constraints = Convex.Constraint[]

    # Precoder definition, power constraints and weighted MMSE objective contribution
    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        
        for j = 1:channel.I
            for l in served_MS_ids(j, assignment)
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.W[l]*state.U[l]')*channel.H[l,i])
            end
        end
        Gamma_r = real(full(Gamma)); Gamma_i = imag(full(Gamma))
        Gamma_ext = Symmetric(hvcat((2,2), Gamma_r, -Gamma_i, Gamma_i, Gamma_r))
        Gamma_ext_sqrtm = sqrtm(Gamma_ext)

        # Precoders and weighted MMSE objective contribution for served users
        used_power = Convex.Constant(0)
        for k in served_MS_ids(i, assignment)
            # Desired channel
            G = channel.H[k,i]'*state.U[k]*state.W[k]
            Gr = real(G); Gi = imag(G)
            G_ext = vcat(Gr, Gi) # only retain real part in multiplication

            Vs[k] = Convex.Variable(2*channel.Ms[i], size(state.W[k], 1))
            objective += Convex.sum_squares(Gamma_ext_sqrtm*Vs[k]) - 2*trace(G_ext'*Vs[k])
            used_power += Convex.sum_squares(Vs[k])
        end
        push!(constraints, used_power <= Ps[i])
    end

    # Interference subspace basis nuclear norm regularization
    if aux_params["NuclearNormHeuristicMosek:perform_regularization"]
        for i = 1:channel.I
            for k = served_MS_ids(i, assignment)
                J_ext = Convex.Variable(2*channel.Ms[i], sum(ds) - ds[k]); j_offset = 0

                for j = 1:channel.I
                    Hr = real(channel.H[k,j]); Hi = imag(channel.H[k,j])
                    H_ext = hvcat((2,2), Hr, -Hi, Hi, Hr)

                    for l in served_MS_ids_except_me(k, j, assignment)
                        push!(constraints, J_ext[:, (1 + j_offset):(j_offset + ds[l])] == H_ext*Vs[j])
                        j_offset += ds[l]
                    end
                end

                Ur = real(state.U[k]); Ui = imag(state.U[k])
                U_ext = vcat(Ur, Ui)

                IntfNN_obj = Convex.nuclear_norm(U_ext'*J_ext)

                objective += aux_params["rho"]*IntfNN_obj
            end
        end
    end

    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem, aux_params["NuclearNormHeuristicMosek:solver"])
    if problem.status == :Optimal
        for i = 1:channel.I
            for k in served_MS_ids(i, assignment)
                state.V[k] = Vs[k].value[1:channel.Ms[i],:] + im*Vs[k].value[channel.Ms[i]+1:end,:]
            end
        end
        objective = problem.optval
    else
        warn("Problem with Convex.jl in update_BSs!")
        objective = 0
    end

    return objective
end
