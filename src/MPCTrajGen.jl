# module MPCTrajGen

using Parameters
using Dierckx
using LinearAlgebra

# export  WrapPosNegPi,
#         Wrap2Pi,
#         WrapPosNeg180,
#         Wrap360,
#         GetΔtVec
include("helpers.jl")

# export  State,
#         state_vec,
#         generate_trajectory!
include("motion_model.jl")

# export  optimise_trajectory
include("mpc_traj_gen.jl")

# end
