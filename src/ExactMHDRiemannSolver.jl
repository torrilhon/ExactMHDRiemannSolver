"""
    ExactMHDRiemannSolver

Exact solver for Riemann problems of 1D ideal magnetohydrodynamics (γ-law gas),
following M. Torrilhon, "Exact Solver and Uniqueness Conditions for Riemann
Problems of Ideal Magnetohydrodynamics", SAM Research Report 2002-06, ETH Zürich.

Version 0.1 computes regular solutions only (fast/slow Lax shocks, fast/slow fans,
rotational discontinuities, contact). Every returned solution is checked
independently; inputs outside the supported domain are refused with a reason.
"""
module ExactMHDRiemannSolver

using LinearAlgebra
using Logging
using Printf
using StaticArrays
using ForwardDiff
using Roots
using NonlinearSolve
using OrdinaryDiffEqVerner
using SciMLBase

export RiemannProblem, SolverOptions, RiemannSolution, RetCode
export Success, InvalidInput, Unsupported, RegularLimit, NoConvergence, CheckFailed
export solve, sample, wavetable, write_csv, check, successful

include("types.jl")
include("eos.jl")
include("canonical.jl")
include("shocks.jl")
include("fans.jl")
include("side.jl")
include("check.jl")
include("solve.jl")
include("output.jl")

end
