include("../src/MPCTrajGen.jl")

# using .MPCTrajGen
using Statistics
using Plots
gr()
using Printf
using Profile

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

function test_mpc_trajectory(save::Bool=false)
    initial = State(x=171.04, y=38.98, θ=deg2rad(84.79), v=10.13, β=deg2rad(-5.21))
    target = State(x=150.96, y=81.01, θ=deg2rad(172.18), v=10.0, β=0.0)

	params = zeros(6)

	states, params = optimise_trajectory(target, params, hyperparams,
											initial=initial, verbose=false)
    # if !isnothing(states)
    #     fig = plot()
	#
    #     tu, tv = target.v * cos(target.θ), target.v * sin(target.θ)
    #     quiver!([target.x], [target.y], quiver=([tu],[tv]), color=:forestgreen)
    #     scatter!([target.x], [target.y], markershape=:star5,
    #                                         markersize=10,
    #                                         markercolor=:green,
    #                                         legend=false)
	#
    #     plot_traj(states)
    #     display(fig)
	#
    #     if save
    #         mkpath("./images/")
    #         filename = @sprintf "res_%05d.png" length(readdir("./images/"))
    #         savefig(fig, "./images/" * filename)
    #     end
    # end
end

function test_initial_state(save::Bool=false)
    noise = nothing

    x_i = 0.0
    y_i = 0.0
    θ_i = 0.0
    # for θ_i in -π/4 : π/12 : π/4
    for vel_i in range(0.0, stop=15.0, length=16)
        for β_i in range(-0.6, stop=0.6, length=11)
            fig = plot()

            a_allowed = [((-10.0)/(time)), ((10.0)/(time))]
            δ_allowed = [((-π/3)/(time)), ((π/3)/(time))]
            # a_allowed = [((0.0-vel_i)/(timestep*n)), ((10.0-vel_i)/(timestep*n))]
            # δ_allowed = [(((-π/3)-β_i)/(timestep*n)), (((π/3)-β_i)/(timestep*n))]

            if param_noise
                noise = OrnsteinUhlenbeckNoise(zeros(5), 1.0, θ=1.0)
            end

			for a_init in a_allowed
            	for a_mid in a_allowed
	                for a_fin in a_allowed
	                    for δ_init in δ_allowed
	                        for δ_mid in δ_allowed
	                            for δ_fin in δ_allowed
	                                s = State(x=x_i, y=y_i, θ=θ_i, v=vel_i, β=β_i)
	                                params = [a_init, a_mid, a_fin, δ_init, δ_mid, δ_fin]
	                                if param_noise
	                                    Δt = 1.0/(length(a_allowed)^2 * length(δ_allowed)^3)
	                                    p_noise = OrnsteinUhlenbeckNoise!(noise, Δt)
	                                    params += p_noise
	                                end

	                                calc_and_plot_traj!(s, params, hyperparams, plot=save)
									# display(fig)
	                            end # δ_fin
	                        end # δ_mid
	                    end # δ_init
					end # a_fin
				end # a_mid
            end # a_init

            if save
                mkpath("./images/")
                filename = @sprintf "(%d,%d).png" vel_i rad2deg(β_i)
				@printf("Save: %s\n", "./images/" * filename)
                savefig(fig, "./images/" * filename)
            end
        end # β_i
    end # vel_i
    # end # θ_i
end

function δ_tests()
    vel_i = 5.0
    β_i = -π/4
    run = 0
    for δ_init in -π/2 : π/4 : π/2
        for δ_mid in -π/2 : π/4 : π/2
            for δ_fin in -π/2 : π/4 : π/2
                fig = plot()

                s = State(x=0.0, y=0.0, θ=0.0, v=vel_i, β=β_i)
                params = [2.0, 0.5, -1.0, δ_init, δ_mid, δ_fin]
                calc_and_plot_traj!(s, params, hyperparams, plot=true)

                run += 1

                mkpath("./images/")
                filename = @sprintf "test%03d_(%2.2d,%2.2d,%2.2d).png" run rad2deg(δ_init) rad2deg(δ_mid) rad2deg(δ_fin)
                savefig(fig, "./images/" * filename)
            end
        end
    end
end

function a_tests()
    vel_i = 5.0
    β_i = -π/4
    run = 0
    for a_init in -1.0 : 0.3 : 0.5
        for a_mid in -1.0 : 0.3 : 0.5
            for a_fin in -1.0 : 0.3 : 0.5
                fig = plot()

                s = State(x=0.0, y=0.0, θ=0.0, v=vel_i, β=β_i)
                params = [a_init, a_mid, a_fin, π/4, 0.0, -π/4]
                calc_and_plot_traj!(s, params, hyperparams, plot=true)

                run += 1

                mkpath("./images/")
                filename = @sprintf "test%03d_(%1.1f,%1.1f,%1.1f).png" run a_init a_mid a_fin
                savefig(fig, "./images/" * filename)
            end
        end
    end
end

test_mpc_trajectory(false)
# @profile (for i = 1:1000; test_mpc_trajectory(false); end)
# Profile.print()
@time begin
  for i = 1:1000
	  test_mpc_trajectory(false)
  end
end
# test_initial_state(true)

# δ_tests()
# a_tests()
