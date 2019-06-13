include("../src/MPCTrajGen.jl")

# using .MPCTrajGen
using Statistics
using Plots
gr()
using Printf
using JLD2, FileIO

# Hyperparamaters
base = 2.0
ord = 10.0
n = 26
timestep = 0.2
time = 5.0
interp = 2
log = false
knot = 0.0
init_knot = false
param_noise = false
action_noise = false
hyperparams = [base, ord, n, timestep, interp, log, knot]

include("test_helpers.jl")

function state_from_dict_idx(traj::Dict{Any, Any}, p::Int64)
    return [traj["x"][p], traj["y"][p], traj["θ"][p], traj["v"][p], traj["β"][p]]
end

function get_traj_distance(traj1::Dict{Any, Any}, traj2::Dict{Any, Any})
    pts = length(traj1["x"])

    dist = 0.0
    Δβ = 0.0
    for p in 1:pts
        s1 = state_from_dict_idx(traj1, p)
        s2 = state_from_dict_idx(traj2, p)
        s1_euc = [s1[1], s1[2], s1[4]*cos(s1[3]), s1[4]*sin(s1[3])]
        s2_euc = [s2[1], s2[2], s2[4]*cos(s2[3]), s2[4]*sin(s2[3])]
        dist += norm(s1_euc - s2_euc)
        Δβ += abs(WrapPosNeg180(rad2deg(s1[5] - s2[5])))
    end
    return dist
end

function green_kelly(init_set::Array{Any,1}, N::Int64)
    seed = init_set[findall(traj->traj["params"] == [3.5, 3.5, 3.5, 0.0, 0.0, 0.0], init_set)[1]]
    final_set = []
    push!(final_set, seed)

    while length(final_set) < N
        min_distances = []
        for i in 1:length(init_set)
            distances = []
            for j in 1:length(final_set)
                d = get_traj_distance(init_set[i], final_set[j])
                push!(distances, d)
            end
            min_d = minimum(distances)
            push!(min_distances, min_d)
        end
        add_id = argmax(min_distances)

        push!(final_set, init_set[add_id])
    end
    return final_set
end

function rollout(N::Int64; saveimg::Bool=false, M::Int64=n)
    x_i = 0.0
    y_i = 0.0
    θ_i = 0.0

    for vel_i in range(15.0, stop=0.0, length=16)
        for β_i in range(-0.6, stop=0.6, length=11)
            all_trajs = []
            for a_init in range(-4.0, stop=3.5, length=5)
                for a_mid in range(-4.0, stop=3.5, length=5)
                    for a_fin in range(-4.0, stop=3.5, length=5)
                        for δ_init in range(-0.6, stop=0.6, length=5)
                            for δ_mid in range(-0.6, stop=0.6, length=5)
                                for δ_fin in range(-0.6, stop=0.6, length=5)
                                    traj = Dict()
                                    s = State(x=x_i, y=y_i, θ=θ_i, v=vel_i, β=β_i)
                                    params = [a_init, a_mid, a_fin, δ_init, δ_mid, δ_fin]

                                    states, a1, δ1 = calc_and_plot_traj!(s, params, hyperparams, plot=saveimg)
                                    x, y, θ, v, β = get_coord_trajectories(states)

                                    traj["x"] = x
                                    traj["y"] = y
                                    traj["θ"] = θ
                                    traj["v"] = v
                                    traj["β"] = β
                                    traj["params"] = params
                                    traj["cmd"] = [a1, δ1]
                                    push!(all_trajs, traj)
                                end # δ_fin
                            end # δ_mid
                        end # δ_init
                    end # a_fin
                end # a_mid
            end # a_init

            trajset = green_kelly(all_trajs, N)

            datadir = @sprintf "./data/gk%d/" N
            mkpath(datadir * "images/")
            mkpath(datadir * "mats/")
            if saveimg
                filename = @sprintf "images/(%d,%d).png" vel_i rad2deg(β_i)
                fig = plot()
                for i in 1:length(trajset)
                    plot_traj(trajset[i]["x"], trajset[i]["y"], trajset[i]["θ"], trajset[i]["v"], trajset[i]["β"])
                end
                savefig(fig, datadir * filename)
            end

            TRAJSMAT = zeros(N, M*5 + 2)
            for i in 1:length(trajset)
                for j in 1:M
                    TRAJSMAT[i, collect(1:5) .+ (j-1)*5] = [trajset[i]["x"][j], trajset[i]["y"][j], trajset[i]["θ"][j], trajset[i]["v"][j], trajset[i]["β"][j]]
                end
                TRAJSMAT[i, end-1:end] = trajset[i]["cmd"]
            end
            filename = @sprintf "mats/(%d,%d).jld2" vel_i rad2deg(β_i)
            filename = datadir * filename
            @save filename TRAJSMAT
        end # β_i
    end # vel_i
    # end # θ_i
end

rollout(30, saveimg=true, M=5)
