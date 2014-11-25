immutable NuclearNormHeuristicState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1}
    V::Array{Matrix{Complex128},1}
end

function NuclearNormHeuristic(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    check_and_defaultize_settings!(settings, NuclearNormHeuristicState)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    state = NuclearNormHeuristicState(
        Array(Matrix{Complex128}, K),
        unity_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings))
    objective = Float64[]
    logdet_rates = Array(Float64, K, maximum(ds), settings["max_iters"])
    MMSE_rates = Array(Float64, K, maximum(ds), settings["max_iters"])
    allocated_power = Array(Float64, K, maximum(ds), settings["max_iters"])

    iters = 0; conv_crit = Inf
    while iters < settings["max_iters"]
        update_MSs!(state, channel, sigma2s, cell_assignment, settings)
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

function check_and_defaultize_settings!(settings, ::Type{NuclearNormHeuristicState})
    # Global settings
    check_and_defaultize_settings!(settings)

    # Local settings
    if !haskey(settings, "NuclearNormHeuristic:perform_regularization")
        settings["NuclearNormHeuristic:perform_regularization"] = true
    end
    if !haskey(settings, "NuclearNormHeuristic:solver")
        if settings["NuclearNormHeuristic:perform_regularization"] == false
            # ECOS gives sufficiently accurate results for the BCD to converge.
            settings["NuclearNormHeuristic:solver"] = ECOS.ECOSMathProgModel()
        else
            # Empirically, it seems that SCS does not give sufficiently accurate
            # results to reproduce the closed-form WMMSE solutions when the
            # regularization is turned off. I probably need another solver, and
            # thus I have to wait for support for this in Convex.jl.
            settings["NuclearNormHeuristic:solver"] = SCS.SCSMathProgModel()
        end
    end
    if !haskey(settings, "NuclearNormHeuristic:regularization_factor")
        settings["NuclearNormHeuristic:regularization_factor"] = 0
    end
end

function update_MSs!(state::NuclearNormHeuristicState,
    channel::SinglecarrierChannel, sigma2s::Vector{Float64},
    cell_assignment::CellAssignment, settings)

    ds = [ size(state.W[k], 1) for k = 1:length(state.W) ]

    for i = 1:channel.I
        for k = served_MS_ids(i, cell_assignment)
            # Received covariance and interference subspace basis
            Phi = Hermitian(complex(sigma2s[k]*eye(channel.Ns[k])))
            J_ext = Array(Float64, 2*channel.Ns[k], sum(ds) - ds[k]); j_offset = 0
            for j = 1:channel.I
                Hr = real(channel.H[k,j]); Hi = imag(channel.H[k,j])
                H_ext = hvcat((2,2), Hr, -Hi, Hi, Hr)
                for l in served_MS_ids(j, cell_assignment)
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

            if settings["NuclearNormHeuristic:perform_regularization"]
                # Nuclear norm term reformulated as inspired by Convex.jl
                A = Convex.Variable(sum(ds) - ds[k], sum(ds) - ds[k])
                B = Convex.Variable(ds[k], ds[k])
                IntfNN_obj = 0.5*(trace(A) + trace(B))
                IntfNN_constr = Convex.isposdef([A J_ext'*Us;Us'*J_ext B])

                problem = Convex.minimize(MSE + settings["NuclearNormHeuristic:regularization_factor"]*IntfNN_obj, IntfNN_constr)
            else
                # Solve standard MSE problem
                problem = Convex.minimize(MSE)
            end

            Convex.solve!(problem, settings["NuclearNormHeuristic:solver"])
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

function update_BSs!(state::NuclearNormHeuristicState,
    channel::SinglecarrierChannel, Ps::Vector{Float64},
    cell_assignment::CellAssignment, settings)

    ds = [ size(state.W[k], 1) for k = 1:length(state.W) ]

    Vs = Array(Convex.Variable, channel.K)
    objective = Convex.Constant(0)
    constraints = Convex.Constraint[]

    # Precoder definition, power constraints and weighted MMSE objective contribution
    for i = 1:channel.I
        # Virtual uplink covariance
        Gamma = Hermitian(complex(zeros(channel.Ms[i],channel.Ms[i])))
        
        for j = 1:channel.I
            for l in served_MS_ids(j, cell_assignment)
                Gamma += Hermitian(channel.H[l,i]'*(state.U[l]*state.W[l]*state.U[l]')*channel.H[l,i])
            end
        end
        Gamma_r = real(full(Gamma)); Gamma_i = imag(full(Gamma))
        Gamma_ext = Symmetric(hvcat((2,2), Gamma_r, -Gamma_i, Gamma_i, Gamma_r))
        Gamma_ext_sqrtm = sqrtm(Gamma_ext)

        # Precoders and weighted MMSE objective contribution for served users
        used_power = Convex.Constant(0)
        for k in served_MS_ids(i, cell_assignment)
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
    if settings["NuclearNormHeuristic:perform_regularization"]
        for i = 1:channel.I
            for k = served_MS_ids(i, cell_assignment)
                J_ext = Convex.Variable(2*channel.Ms[i], sum(ds) - ds[k]); j_offset = 0

                for j = 1:channel.I
                    Hr = real(channel.H[k,j]); Hi = imag(channel.H[k,j])
                    H_ext = hvcat((2,2), Hr, -Hi, Hi, Hr)

                    for l in served_MS_ids_except_me(k, j, cell_assignment)
                        push!(constraints, J_ext[:, (1 + j_offset):(j_offset + ds[l])] == H_ext*Vs[j])
                        j_offset += ds[l]
                    end
                end

                Ur = real(state.U[k]); Ui = imag(state.U[k])
                U_ext = vcat(Ur, Ui)

                # Nuclear norm term reformulated as inspired by Convex.jl
                A = Convex.Variable(sum(ds) - ds[k], sum(ds) - ds[k])
                B = Convex.Variable(ds[k], ds[k])
                IntfNN_obj = 0.5*(trace(A) + trace(B))
                IntfNN_constr = Convex.isposdef([A J_ext'*U_ext;U_ext'*J_ext B])

                objective += settings["NuclearNormHeuristic:regularization_factor"]*IntfNN_obj
                push!(constraints, IntfNN_constr)
            end
        end
    end

    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem, settings["NuclearNormHeuristic:solver"])
    if problem.status == :Optimal
        for i = 1:channel.I
            for k in served_MS_ids(i, cell_assignment)
                state.V[k] = Vs[k].value[1:channel.Ms[i],:] + im*Vs[k].value[channel.Ms[i]+1:end,:]
            end
        end
        objective = problem.optval
    else
        error("Problem with Convex.jl in update_BSs!")
        objective = 0
    end

    return objective
end
