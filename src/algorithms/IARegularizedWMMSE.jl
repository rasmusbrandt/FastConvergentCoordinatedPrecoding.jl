immutable IARegularizedWMMSEState
    U::Array{Matrix{Complex128},1} # receive filters
    W::Array{Hermitian{Complex128},1} # MSE weights
    V::Array{Matrix{Complex128},1} # precoders
end

function IARegularizedWMMSE(channel::SinglecarrierChannel, network::Network,
    cell_assignment::CellAssignment, settings=Dict())

    settings = check_and_defaultize_settings(IARegularizedWMMSEState,
                                             settings, channel, network)

    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network)

    state = IARegularizedWMMSEState(
        Array(Matrix{Complex128}, channel.K),
        initial_MSE_weights(ds),
        initial_precoders(channel, Ps, sigma2s, ds, cell_assignment, settings))
    rates = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])
    objective = Array(Float64, channel.K, maximum(ds), settings["stop_crit"])

    objective[:,:,1] = 0
    for iter = 1:(settings["stop_crit"]-1)
        update_MSs!(state, channel, sigma2s, cell_assignment, settings)
        rates[:,:,iter] = calculate_rates(state)
        objective[:,:,iter+1] = update_BSs!(state, channel, Ps, cell_assignment, settings)
    end
    update_MSs!(state, channel, sigma2s, cell_assignment, settings)
    rates[:,:,end] = calculate_rates(state)

    if settings["output_protocol"] == 1
        return [ "rates" => rates,
                 "objective" => objective ]
    elseif settings["output_protocol"] == 2
        return [ "rates" => rates[:,:,end],
                 "objective" => objective[:,:,end] ]
    end
end

function check_and_defaultize_settings(::Type{IARegularizedWMMSEState},
    settings, channel::SinglecarrierChannel, network::Network)

    settings = copy(settings)

    # Global settings and consistency checks
    if !haskey(settings, "user_priorities")
        settings["user_priorities"] = ones(channel.K)
    end
    if !haskey(settings, "output_protocol")
        settings["output_protocol"] = 1
    end
    if !haskey(settings, "stop_crit")
        settings["stop_crit"] = 20
    end
    if !haskey(settings, "initial_precoders")
        settings["initial_precoders"] = "dft"
    end
    if settings["output_protocol"] != 1 && settings["output_protocol"] != 2
        error("Unknown output protocol")
    end

    # Local settings
    if !haskey(settings, "IARegularizedWMMSE:perform_regularization")
        settings["IARegularizedWMMSE:perform_regularization"] = true
    end
    if !haskey(settings, "IARegularizedWMMSE:solver")
        if settings["IARegularizedWMMSE:perform_regularization"] == false
            # ECOS gives sufficiently accurate results for the BCD to converge.
            settings["IARegularizedWMMSE:solver"] = ECOS.ECOSMathProgModel()
        else
            # Empirically, it seems that SCS does not give sufficiently accurate
            # results to reproduce the closed-form WMMSE solutions when the
            # regularization is turned off. I probably need another solver, and
            # thus I have to wait for support for this in Convex.jl.
            settings["IARegularizedWMMSE:solver"] = SCS.SCSMathProgModel()
        end
    end
    if !haskey(settings, "IARegularizedWMMSE:regularization_factor")
        settings["IARegularizedWMMSE:regularization_factor"] = 0
    end

    return settings
end

function update_MSs!(state::IARegularizedWMMSEState,
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
                    herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)

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

            if settings["IARegularizedWMMSE:perform_regularization"]
                # Nuclear norm term reformulated as inspired by Convex.jl
                A = Convex.Variable(sum(ds) - ds[k], sum(ds) - ds[k])
                B = Convex.Variable(ds[k], ds[k])
                IntfNN_obj = 0.5*(trace(A) + trace(B))
                IntfNN_constr = Convex.isposdef([A J_ext'*Us;Us'*J_ext B])

                problem = Convex.minimize(MSE + settings["IARegularizedWMMSE:regularization_factor"]*IntfNN_obj, IntfNN_constr)
            else
                # Solve standard MSE problem
                problem = Convex.minimize(MSE)
            end

            Convex.solve!(problem, settings["IARegularizedWMMSE:solver"])
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

function update_BSs!(state::IARegularizedWMMSEState,
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
    if settings["IARegularizedWMMSE:perform_regularization"]
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

                objective += settings["IARegularizedWMMSE:regularization_factor"]*IntfNN_obj
                push!(constraints, IntfNN_constr)
            end
        end
    end

    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem, settings["IARegularizedWMMSE:solver"])
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

function calculate_rates(state::IARegularizedWMMSEState)
    K = length(state.W)
    ds = Int[ size(state.W[k], 1) for k = 1:K ]; max_d = maximum(ds)

    rates = Array(Float64, K, max_d)

    for k = 1:K
        # W is p.d., so we should only get real eigenvalues. Numerically we may
        # get some imaginary noise however. Also, numerically the eigenvalues
        # may be less that 1, so we need to handle that to not get negative
        # rates.
        r = log2(max(1, real(eigvals(state.W[k]))))

        if ds[k] < max_d
            rates[k,:] = cat(1, r, zeros(Float64, max_d - ds[k]))
        else
            rates[k,:] = r
        end
    end

    return rates
end
