const L_R = 2.0
const L_F = 2.0

@with_kw mutable struct State
	x::Float64 = 0.0
	y::Float64 = 0.0
	θ::Float64 = 0.0 # in radians
	v::Float64 = 0.0
	β::Float64 = 0.0 # in radians
end
Base.copy(s::State) = State(x=s.x, y=s.y, θ=s.θ, v=s.v, β=s.β)

function reset_state!(s::State)
	s.x = 0.0
	s.y = 0.0
	s.θ = 0.0
	s.v = 0.0
	s.β = 0.0
end

function state_vec(s::State)
    v = [s.x, s.y, s.θ, s.v, s.β]
    return v
end

function update!(s::State, a::Float64, δf::Float64, Δt::Float64)
	s.x += s.v * cos(s.θ + s.β) * Δt
	s.y += s.v * sin(s.θ + s.β) * Δt
	s.θ += (s.v / L_R) * sin(s.β) * Δt
	s.θ = WrapPosNegPi(s.θ)

	s.v += a * Δt
	s.v = max(s.v, 0.0)
	s.v = min(s.v, 15.0)
	s.β = atan( (L_R / (L_F + L_R)) * tan(δf * Δt))
	s.β = WrapPosNegPi(s.β)
	s.β = clamp(s.β, -π/4, π/4)
end

function get_knots(params::Vector{Float64}, knot::Float64, init_knot::Bool)
	if length(params) == 4
		if init_knot
			a_mid, a_fin, δ_mid, δ_fin = params
			a_knots = [knot, a_mid, a_fin]
			δ_knots = [knot, δ_mid, δ_fin]
		else
			a_init, a_fin, δ_init, δ_fin = params
			a_knots = [a_init, knot, a_fin]
			δ_knots = [δ_init, knot, δ_fin]
		end

		return a_knots, δ_knots
	elseif length(params) == 5
		a_init, a_fin, δ_init, δ_mid, δ_fin = params
		a_knots = [a_init, knot, a_fin]
		δ_knots = [δ_init, δ_mid, δ_fin]
		return a_knots, δ_knots
	else
		a_knots = params[1:3]
		δ_knots = params[4:6]
		return a_knots, δ_knots
	end
end

function generate_trajectory!(s::State, params::Vector{Float64},
								hyperparams::Vector{Float64};
								return_last::Bool=false,
								init_knot::Bool=true,
								noisy::Bool=false)
	base, ord, n, timestep, interp, log, knot = hyperparams
	# t_knots = [0.0, time/2.0, time]
	times = GetΔtVec(Bool(log), base, ord, Int64(n), timestep)
	t_knots = [times[1], times[div(Int64(n), 2)], times[end]]
	a_knots, δ_knots = get_knots(params, knot, init_knot)

	a_spline = Spline1D(t_knots, a_knots; k=Int64(interp))
	δ_spline = Spline1D(t_knots, δ_knots; k=Int64(interp))

	a_interm = a_spline(times)
	δ_interm = δ_spline(times)

	noise_a = OrnsteinUhlenbeckNoise([0.0], 0.2, θ=1.0)
	noise_δ = OrnsteinUhlenbeckNoise([0.0], 0.2, θ=1.0)

	states = [copy(s)]
	for i in 1:length(times)-1
		Δt = times[i+1] - times[i]
		if noisy
			a_noise = OrnsteinUhlenbeckNoise!(noise_a, Δt)
			δ_noise = OrnsteinUhlenbeckNoise!(noise_δ, Δt)
			update!(s, a_interm[i] + a_noise[1], δ_interm[i] + δ_noise[1], Δt)
		else
			update!(s, a_interm[i], δ_interm[i], Δt)
		end
		push!(states, copy(s))
	end

	if return_last
		return s
	end

	return s, states, a_interm[1], δ_interm[1]
end

function generate_last_state!(s::State, params::Vector{Float64},
								hyperparams::Vector{Float64})
	return generate_trajectory!(s, params, hyperparams::Vector{Float64},
									return_last=true)
end
