const OPTIM_ITERS = 50
const MIN_COST = 0.1
const EARLY_TERM = 1e-9
const STATE_DIM = 5
const Δs_deriv = [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]
# const W = [5.0, 5.0, 5.0, 1.0, 1.0] # [x, y, θ, v, β]

# TODO: Figure out the best cost function
function state_cost(s::State; p::Real=2)
	# return norm(state_vec(s), p)
	sv = state_vec(s)
    # return norm(sv[[1,2,4]], p) + (Wrap2Pi(sv[3]) + Wrap2Pi(sv[5])) * 0.5
	return norm(sv[1:3], p)
	# return sum((sv).^p .* W)
end

function traj_cost!(s::State, params::Vector{Float64}, target::State,
					hyperparams::Vector{Float64}; bypass::Bool=false)
	cost = 0.0
	if bypass
		s = generate_last_state!(s, params, hyperparams)
		Δs = state_diff(target, s)
		cost += state_cost(Δs)
		return cost
	end

	base, ord, n, timestep, _, log, _ = hyperparams
	times = GetΔtVec(Bool(log), base, ord, Int64(n), timestep)

	s, states, _, _ = generate_trajectory!(s, params, hyperparams)
	for i in 1:length(states)-1
		Δs = state_diff(states[i], states[i+1])
		Δt = times[i+1] - times[i]
		cost += state_cost(Δs) * Δt
	end

	Δs = state_diff(s, target)
	cost += state_cost(Δs)

	return cost
end

function state_diff(target::State, curr::State)
    res = State()
    res.x = target.x - curr.x
    res.y = target.y - curr.y
    res.θ = WrapPosNegPi(target.θ - curr.θ)
    res.v = target.v - curr.v
    res.β = WrapPosNegPi(target.β - curr.β)

    return res
end

function get_jacobian_column(target::State, params::Vector{Float64},
								col::Int64, hyperparams::Vector{Float64},
								initial::State)
	new_params = copy(params)

	s = State()
	set_initial_state!(s, initial)

    new_params[col] += Δs_deriv[col]
    s = generate_last_state!(s, new_params, hyperparams)
    Δs_pos = state_diff(target, s)
    Δs_pos = state_vec(Δs_pos)

	set_initial_state!(s, initial)

    new_params[col] -= 2*Δs_deriv[col]
    s = generate_last_state!(s, new_params, hyperparams)
    Δs_neg = state_diff(target, s)
    Δs_neg = state_vec(Δs_neg)

    Δs_col = (Δs_pos - Δs_neg) ./ (2.0 * Δs_deriv[col])
    return Δs_col
end

function calc_jacobian(target::State, params::Vector{Float64},
						hyperparams::Vector{Float64},
						initial::State)
    J = zeros(STATE_DIM, length(params))
    for col in 1:length(params)
        J[:, col] = get_jacobian_column(target, params,
											col, hyperparams, initial)
    end
    return J
end

function α_line_search(Δp::Vector{Float64}, params::Vector{Float64},
                                target::State, hyperparams::Vector{Float64},
								initial::State)
    mincost = Inf
    s = State()
	set_initial_state!(s, initial)

    best_α = nothing
    for α in 0.1:0.1:1.0
		if α != 0
	        test_params = params .+ (α .* Δp)
	    	c = traj_cost!(s, test_params, target, hyperparams, bypass=true)

	        if c < mincost
	            mincost = c
	            best_α = α
	        end
		end

		set_initial_state!(s, initial)
    end

    return best_α
end

function set_initial_state!(s::State, initial::State)
	if !isnothing(initial)
		s.x = initial.x
		s.y = initial.y
		s.θ = initial.θ
		s.v = initial.v
		s.β = initial.β
	else
		s.x = 0.0
		s.y = 0.0
		s.θ = 0.0
		s.v = 0.0
		s.β = 0.0
	end
end

function optimise_trajectory(target::State, params::Vector{Float64},
                            	hyperparams::Vector{Float64};
								initial::State=nothing,
								verbose::Bool=true)
	states = nothing
	old_cost = Inf
    for i in 1:OPTIM_ITERS
		s = State()
		set_initial_state!(s, initial)

        s, states, _, _ = generate_trajectory!(s, params, hyperparams)
        Δs = state_diff(target, s)
		set_initial_state!(s, initial)
		c = traj_cost!(s, params, target, hyperparams, bypass=true)

		if verbose
			println("Iteration: ", i, " | Parameters: ", params, " | Cost: ", c)
		end
        if c <=  MIN_COST
			if verbose
				println("Optimisation finished!")
			end
            break
        end

		if abs(old_cost - c) < EARLY_TERM
			if verbose
				println("Optimisation in local minima. Terminating!")
			end
			break
		end
		old_cost = c

        J = calc_jacobian(target, params, hyperparams, initial)
        Δp = params
        try
            Δp = -pinv(J) * state_vec(Δs)
        catch e
			if verbose
	            println(J)
	            println("Jacobian is singular. Optimisation failed!")
	            println(e)
			end
            states = nothing
            break
        end

        α = α_line_search(Δp, params, target, hyperparams, initial)
		params += α .* Δp
		params[1:3] .= clamp!(params[1:3], -4.0, 3.5)
		params[4:6] .= clamp!(params[4:6], -0.6, 0.6)
    end

    return states, params
end
