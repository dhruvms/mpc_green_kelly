function get_coord_trajectories(states::Vector{State})
    x, y, θ, v, β = [], [], [], [], []
    for s in states
        push!(x, s.x)
        push!(y, s.y)
        push!(θ, s.θ)
        push!(v, s.v)
        push!(β, s.β)
    end
    return x, y, θ, v, β
end

function plot_traj(states::Vector{State})
    start_state = states[1]
    x_i, y_i, θ_i, vel_i, β_i = state_vec(start_state)

    x, y, θ, vel, β = get_coord_trajectories(states)
    u, v = vel .* cos.(θ), vel .* sin.(θ)
    quiver!(x, y, quiver=(u,v), color=:salmon)
    u, v = (vel ./ 2.0) .* cos.(θ .+ β), (vel ./ 2.0) .* sin.(θ .+ β)
    quiver!(x, y, quiver=(u,v), color=:indianred)
    plot!(x, y, color=:red, legend=false)
    scatter!(x, y, color=:firebrick, legend=false)
    u, v = vel_i * cos(θ_i), vel_i * sin(θ_i)
    quiver!([x_i], [y_i], quiver=([u],[v]), color=:forestgreen)
    u, v = (vel_i / 2.0) * cos(θ_i + β_i), (vel_i / 2.0) * sin(θ_i + β_i)
    quiver!([x_i], [y_i], quiver=([u],[v]), color=:olive)
end

function plot_traj(x::Array{Any,1}, y::Array{Any,1}, θ::Array{Any,1},
                    v::Array{Any,1}, β::Array{Any,1})
    x_i, y_i, θ_i, v_i, β_i = x[1], y[1], θ[1], v[1], β[1]

    plot!(x, y, color=:red, legend=false)
    scatter!(x, y, color=:firebrick, legend=false)

    u, v = v .* cos.(θ), v .* sin.(θ)
    quiver!([x[end]], [y[end]], quiver=([u[end]],[v[end]]), color=:salmon)
    u, v = (v ./ 2.0) .* cos.(θ .+ β), (v ./ 2.0) .* sin.(θ .+ β)
    quiver!([x[end]], [y[end]], quiver=([u[end]],[v[end]]), color=:indianred)

    u, v = v_i * cos(θ_i), v_i * sin(θ_i)
    quiver!([x_i], [y_i], quiver=([u],[v]), color=:forestgreen)
    u, v = (v_i / 2.0) * cos(θ_i + β_i), (v_i / 2.0) * sin(θ_i + β_i)
    quiver!([x_i], [y_i], quiver=([u],[v]), color=:olive)
end

function calc_and_plot_traj!(s::State, params::Vector{Float64},
                                hyperparams::Vector{Float64}; plot::Bool=true,
                                action_noise::Bool=false)
    s, states, a1, δ1 = generate_trajectory!(s, params, hyperparams,
                                        noisy=action_noise)

    if plot
        plot_traj(states)
    end
    return states, a1, δ1
end
