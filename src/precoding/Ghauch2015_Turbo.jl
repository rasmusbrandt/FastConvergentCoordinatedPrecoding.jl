immutable Ghauch2015_TurboState
    U::Array{Matrix{Complex128},1}
    W::Array{Hermitian{Complex128},1} # these are only used for rate calculations
    V::Array{Matrix{Complex128},1}
end

function Ghauch2015_Turbo(channel, network)
    assignment = get_assignment(network)

    K = get_no_MSs(network)
    Ps = get_transmit_powers(network)
    sigma2s = get_receiver_noise_powers(network)
    ds = get_no_streams(network); max_d = maximum(ds)
    alphas = get_user_priorities(network)
    aux_params = get_aux_precoding_params(network)

    aux_params = get_aux_precoding_params(network)
    @defaultize_param! aux_params "Ghauch2015_Turbo:bisection_max_iters" 5e1
    @defaultize_param! aux_params "Ghauch2015_Turbo:bisection_tolerance" 1e-3

    state = Ghauch2015_TurboState(
        Array(Matrix{Complex128}, K),
        Array(Hermitian{Complex128}, K),
        initial_precoders(channel, Ps, sigma2s, ds, assignment, aux_params))
    objective = Float64[]
    logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_logdet_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    weighted_MMSE_rates = Array(Float64, K, max_d, aux_params["max_iters"])
    allocated_power = Array(Float64, K, max_d, aux_params["max_iters"])

    # Initial point for Delta, Phi, Lambda, Gamma
    Delta = Array(Matrix{Complex128}, channel.K); Phi = Array(Matrix{Complex128}, channel.K)
    Lambda = Array(Matrix{Complex128}, channel.K); Gamma = Array(Matrix{Complex128}, channel.K)
    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        q, _ = qr(randn(channel.Ns[k], channel.Ns[k]) + im*randn(channel.Ns[k], channel.Ns[k]), thin=false)
        Delta[k] = q[:, 1:ds[k]]; Phi[k] = q[:, ds[k]+1:end] # Delta'*Phi = 0

        q, _ = qr(randn(channel.Ms[i], channel.Ms[i]) + im*randn(channel.Ms[i], channel.Ms[i]), thin=false)
        Lambda[k] = q[:, 1:ds[k]]; Gamma[k] = q[:, ds[k]+1:end] # Lambda'*Gamma = 0
    end; end

    iters = 0; conv_crit = Inf
    while iters < aux_params["max_iters"]
        update_MSs!(state, channel, Ps, sigma2s, ds, assignment, aux_params, Delta, Phi)
        iters += 1

        # Results after this iteration
        logdet_rates[:,:,iters] = calculate_logdet_rates(state)
        push!(objective, sum(logdet_rates[:,:,iters]))
        MMSE_rates[:,:,iters] = calculate_MMSE_rates(state)
        weighted_logdet_rates[:,:,iters] = calculate_weighted_logdet_rates(state, alphas)
        weighted_MMSE_rates[:,:,iters] = calculate_weighted_MMSE_rates(state, alphas)
        allocated_power[:,:,iters] = calculate_allocated_power(state)

        # Check convergence
        if iters >= 2
            conv_crit = abs(objective[end] - objective[end-1])/abs(objective[end-1])
            if conv_crit < aux_params["stop_crit"]
                Lumberjack.debug("Ghauch2015_Turbo converged.",
                    [ :no_iters => iters,
                      :final_objective => objective[end],
                      :conv_crit => conv_crit,
                      :stop_crit => aux_params["stop_crit"],
                      :max_iters => aux_params["max_iters"] ]
                )
                break
            end
        end

        # Begin next iteration, unless the loop will end
        if iters < aux_params["max_iters"]
            update_BSs!(state, channel, Ps, sigma2s, ds, assignment, aux_params, Lambda, Gamma)
        end
    end
    if iters == aux_params["max_iters"]
        Lumberjack.debug("Ghauch2015_Turbo did NOT converge.",
            [ :no_iters => iters,
              :final_objective => objective[end],
              :conv_crit => conv_crit,
              :stop_crit => aux_params["stop_crit"],
              :max_iters => aux_params["max_iters"] ]
        )
    end

    results = PrecodingResults()
    if aux_params["output_protocol"] == :all_iterations
        results["objective"] = objective
        results["logdet_rates"] = logdet_rates
        results["MMSE_rates"] = MMSE_rates
        results["weighted_logdet_rates"] = weighted_logdet_rates
        results["weighted_MMSE_rates"] = weighted_MMSE_rates
        results["allocated_power"] = allocated_power
    elseif aux_params["output_protocol"] == :final_iteration
        results["objective"] = objective[iters]
        results["logdet_rates"] = logdet_rates[:,:,iters]
        results["MMSE_rates"] = MMSE_rates[:,:,iters]
        results["weighted_logdet_rates"] = weighted_logdet_rates[:,:,iters]
        results["weighted_MMSE_rates"] = weighted_MMSE_rates[:,:,iters]
        results["allocated_power"] = allocated_power[:,:,iters]
    end
    return results
end

function update_MSs!(state::Ghauch2015_TurboState,
    channel::SinglecarrierChannel, Ps, sigma2s, ds, assignment, aux_params, Delta, Phi)

    for i = 1:channel.I; for k in served_MS_ids(i, assignment)
        # Received interference covariance
        Q = Hermitian(complex(zeros(channel.Ns[k], channel.Ns[k])))
        for j = 1:channel.I; for l in served_MS_ids_except_me(k, j, assignment)
            #Q += Hermitian(channel.H[k,j]*(state.V[l]*state.V[l]')*channel.H[k,j]')
            Base.LinAlg.BLAS.herk!(Q.uplo, 'N', complex(1.), channel.H[k,j]*state.V[l], complex(1.), Q.S)
        end; end

        # Turbo iterations
        A = sqrt(Ps[i]/(2*ds[k]))*complex(eye(ds[k], ds[k]))
        B = sqrt(Ps[i]/(2*ds[k]))*complex(eye(channel.Ns[k] - ds[k], ds[k]))
        for turbo_iters = 1:aux_params["turbo_iters"]
            A = Ghauch2015_lemma1(Phi[k]*B, Delta[k], Q, Ps[i] - abs(trace(B'*B)), aux_params)
            B = Ghauch2015_lemma1(Delta[k]*A, Phi[k], Q, Ps[i] - abs(trace(A'*A)), aux_params)
        end

        # Rank reduction
        Ghauch2015_reduce_rank!(A, ds[k])
        Ghauch2015_reduce_rank!(B, ds[k])

        # Final receive filter
        state.U[k] = Delta[k]*A + Phi[k]*B

        # MSE weight for rate calculation (w/ MMSE filter)
        effective_channel = channel.H[k,i]*state.V[k]
        Qtot = Hermitian(Q + effective_channel*effective_channel' + sigma2s[k]*eye(channel.Ns[k]))
        state.W[k] = Hermitian(inv(eye(ds[k]) - effective_channel'*(Qtot\effective_channel)))
    end; end
end

function update_BSs!(state::Ghauch2015_TurboState,
    channel::SinglecarrierChannel, Ps, sigma2s, ds, assignment, aux_params, Lambda, Gamma)

    for i in active_BSs(assignment); for k in served_MS_ids(i, assignment)
        # Virtual uplink interference covariance
        Q = Hermitian(zeros(Complex128, channel.Ms[i], channel.Ms[i]))
        for j = 1:channel.I; for l = served_MS_ids_except_me(k, j, assignment)
            #Q += Hermitian(channel.H[k,i]'*(state.U[k]*state.U[k]')*channel.H[k,i])
            Base.LinAlg.BLAS.herk!(Q.uplo, 'N', complex(1.), channel.H[l,i]'*state.U[l], complex(1.), Q.S)
        end; end

        # Turbo iterations
        C = sqrt(Ps[i]/(2*ds[k]))*complex(eye(ds[k], ds[k]))
        D = sqrt(Ps[i]/(2*ds[k]))*complex(eye(channel.Ms[i] - ds[k], ds[k]))
        for turbo_iters = 1:aux_params["turbo_iters"]
            C = Ghauch2015_lemma1(Gamma[k]*D, Lambda[k], Q, Ps[i] - abs(trace(D'*D)), aux_params)
            D = Ghauch2015_lemma1(Lambda[k]*C, Gamma[k], Q, Ps[i] - abs(trace(C'*C)), aux_params)
        end

        # Rank reduction
        Ghauch2015_reduce_rank!(C, ds[k])
        Ghauch2015_reduce_rank!(D, ds[k])

        # Final transmit filter
        state.V[k] = Lambda[k]*C + Gamma[k]*D
    end; end
end

function Ghauch2015_lemma1(theta, T, Q, zeta, aux_params)
    # Build bisector function
    G = T'*Q*theta
    Z = G*G'
    Y = T'*Q*T
    Y_eigen = eigfact(Y); Y_eigen_values = abs(Y_eigen.values); Y_eigen_vectors = Y_eigen.vectors
    bis_YZY_diag = abs(diag(Y_eigen_vectors'*Z*Y_eigen_vectors))
    f(mu) = sum(bis_YZY_diag./((Y_eigen_values .+ mu).*(Y_eigen_values .+ mu)))

    # Bounds
    mu_lower = -minimum(Y_eigen_values) + 1e-10 # fudge a bit so we don't get indefinite matrix
    mu_upper = vecnorm(G)/sqrt(abs(zeta))

    # Bisect
    iters = 0
    while iters < aux_params["Ghauch2015_Turbo:bisection_max_iters"]
        conv_crit = (zeta - f(mu_upper))/zeta

        if conv_crit < aux_params["Ghauch2015_Turbo:bisection_tolerance"]
            break
        else
            mu = (1/2)*(mu_lower + mu_upper)

            if f(mu) < zeta
                # New point feasible, replace upper point
                mu_upper = mu
            else
                # New point not feasible, replace lower point
                mu_lower = mu
            end
        end

        iters += 1
    end
    if iters == aux_params["Ghauch2015_Turbo:bisection_max_iters"]
        Lumberjack.warn("Ghauch2015_lemma1 bisection: reached max iterations.")
    end

    X = -(Y + mu_upper*eye(Y))\G
end

function Ghauch2015_reduce_rank!(X, d)
    r = rank(X)
    if r < d
        p = trace(X'*X)
        q, _ = qr(X)
        X[:] = sqrt(p/r)*hcat(q[:,1:r], zeros(size(X, 1), d - r))
    end
end
