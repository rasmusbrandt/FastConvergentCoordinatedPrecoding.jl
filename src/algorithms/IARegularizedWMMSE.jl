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

    if !all(ds .== 1)
        error("Currently, Convex.jl does not work with nuclear_norm. Therefore the optimization problems are written with 1-norm formulations, i.e. only working for d = 1. Fix this!")
    end

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
    if !haskey(settings, "IARegularizedWMMSE:regularization_factor")
        settings["IARegularizedWMMSE:regularization_factor"] = 3
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
            J = Array(Complex128, channel.Ns[k], sum(ds) - ds[k]); j_offset = 0
            for j = 1:channel.I
                for l in served_MS_ids(j, cell_assignment)
                    #Phi += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
                    herk!(Phi.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Phi.S)

                    if l != k
                        J[:,(1 + j_offset):(j_offset + ds[l])] = channel.H[k,j]*state.V[l]
                        j_offset += ds[l]
                    end
                end
            end

            # Find receiver
            F = channel.H[k,i]*state.V[k]; Fr = real(F); Fi = imag(F)
            Fext = vcat(Fr, Fi)
            Phi_r = real(full(Phi)); Phi_i = imag(full(Phi))
            Phi_ext = Symmetric(hvcat((2,2), Phi_r, -Phi_i, Phi_i, Phi_r))
            J_r = real(J); J_i = imag(J)
            J_ext = hvcat((2,2), J_r, -J_i, J_i, J_r)

            Us = Convex.Variable(2*channel.Ns[k], ds[k])
            MSE = ds[k] - 2*trace(Fext'*Us) + Convex.sum_squares(sqrtm(Phi_ext)*Us)
            INN = Convex.norm(J_ext'*Us, 1)
            problem = Convex.minimize(MSE + settings["IARegularizedWMMSE:regularization_factor"]*INN)
            Convex.solve!(problem)
            if problem.status == :Optimal
                Ur = Us.value[1:channel.Ns[k],:]; Ui = Us.value[channel.Ns[k]+1:end,:]
                state.U[k] = Ur + im*Ui
            else
                println("Problem with Convex.jl in update_MSs!")
            end

            # Find MSE weight
            Ummse = Phi\F
            state.W[k] = Hermitian((eye(ds[k]) - Ummse'*F)\eye(ds[k]))
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
            Gext = vcat(Gr, Gi)

            Vs[k] = Convex.Variable(2*channel.Ms[i], size(state.W[k], 1))
            objective += Convex.sum_squares(Gamma_ext_sqrtm*Vs[k]) - 2*trace(Gext'*Vs[k])
            used_power += Convex.sum_squares(Vs[k])
        end
        push!(constraints, used_power <= Ps[i])
    end

    # Interference subspace basis contribution
    for i = 1:channel.I
        for k = served_MS_ids(i, cell_assignment)
            J = Convex.Variable(2*channel.Ms[i], sum(ds) - ds[k]); j_offset = 0

            for j = 1:channel.I
                Hr = real(channel.H[k,j]); Hi = imag(channel.H[k,j])
                H_ext = hvcat((2,2), Hr, -Hi, Hi, Hr)

                for l in served_MS_ids_except_me(k, j, cell_assignment)
                    constraints += (J[:, (1 + j_offset):(j_offset + ds[l])] == H_ext*Vs[j])
                    j_offset += ds[l]
                end
            end

            Ur = real(state.U[k]); Ui = imag(state.U[k])
            U_ext = vcat(Ur, Ui)
            objective += settings["IARegularizedWMMSE:regularization_factor"]*Convex.norm(U_ext'*J, 1)
        end
    end

    problem = Convex.minimize(objective, constraints)
    Convex.solve!(problem)
    if problem.status == :Optimal
        for i = 1:channel.I
            for k in served_MS_ids(i, cell_assignment)
                state.V[k] = Vs[k].value[1:channel.Ms[i],:] + im*Vs[k].value[channel.Ms[i]+1:end,:]
            end
        end
        objective = problem.optval
    else
        println("Problem with Convex.jl in update_BSs!")
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
