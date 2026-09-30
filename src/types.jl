# Core types: options, return codes, internal half-states, waves, solutions.

"""
    RetCode

Outcome of [`solve`](@ref). Only `Success` carries a checked solution.

- `Success`        converged and every independent check passed
- `InvalidInput`   non-finite values, ρ ≤ 0, p ≤ 0, γ ≤ 1 or Bx differing between L and R (`:bn_jump`)
- `Unsupported`    input outside the supported domain (see `sol.reason`)
- `RegularLimit`   converged, but the solution sits at a regular-wave limit
- `NoConvergence`  the nonlinear solve failed from every start
- `CheckFailed`    solver converged but the independent checks failed (a bug)
"""
@enum RetCode Success InvalidInput Unsupported RegularLimit NoConvergence CheckFailed

"""
    SolverOptions(; kwargs...)

All thresholds of the solver in one place. Quantities are in canonical units
(ρ_L = p_L = 1, velocities in √(p_L/ρ_L), B in √p_L).
"""
Base.@kwdef struct SolverOptions
    bn_min::Float64 = 1e-3        # min Bn/√p of both input states
    bt_min_larger::Float64 = 1e-4   # min Bt/√p of the input state with the larger value
    bt_min_smaller::Float64 = 1e-6  # min Bt/√p of the input state with the smaller value
    bt_floor::Float64 = 1e-8        # min Bt/√p of all middle states (else RegularLimit)
    ratio_max::Float64 = 1e6      # max ρ_R/ρ_L, p_R/p_L (and inverse)
    δ::Float64 = 1e-6             # safety margin below the slow-shock limit
    s_vac::Float64 = log(1e8)     # cap of slow-fan strength s = log(ρ0/ρ)
    sat::Float64 = 5.0            # tanh argument above which a map counts as saturated
    rho_floor::Float64 = 1e-8     # min ρ, p of middle states (canonical)
    abstol::Float64 = 1e-12       # residual tolerance (∞-norm, canonical)
    maxiters::Int = 60            # nonlinear iterations per start
    homotopy::Bool = true
    homotopy_maxsteps::Int = 400
    time_limit::Float64 = 5.0     # wall-clock budget of the homotopy fallback [s]
    ode_tol::Float64 = 1e-12
    check_tol::Float64 = 1e-8     # relative tolerance of the post-solve checks
    verbose::Bool = false
end

# Internal state along one side of the solution, canonical frame.
# ρ, u (normal velocity), p, bt = |Bt| ≥ 0, φ = direction of Bt, vt = (vy, vz).
struct HState{T}
    ρ::T
    u::T
    p::T
    bt::T
    φ::T
    vt::SVector{2,T}
end
HState(ρ, u, p, bt, φ, vt) = (T = promote_type(typeof(ρ), typeof(u), typeof(p), typeof(bt), typeof(φ), eltype(vt));
    HState{T}(T(ρ), T(u), T(p), T(bt), T(φ), SVector{2,T}(vt)))

fvalue(x::Real) = ForwardDiff.value(x)
fvalue(h::HState) = HState(fvalue(h.ρ), fvalue(h.u), fvalue(h.p), fvalue(h.bt), fvalue(h.φ), fvalue.(h.vt))

struct Ctx
    γ::Float64
    κ::Float64      # (γ+1)/(γ-1)
    Bn::Float64     # canonical, > 0
    opts::SolverOptions
    scale::SVector{3,Float64}   # residual units: pressure, velocity, magnetic field
end
Ctx(γ, Bn, opts, scale = SVector(1.0, 1.0, 1.0)) = Ctx(γ, (γ + 1) / (γ - 1), Bn, opts, scale)

"""
Data of a rarefaction fan, enough to evaluate its interior exactly
(dense ODE output in the fan parameter τ ∈ [0, 1]).
"""
struct FanData
    family::Symbol        # :fast or :slow
    σ::Int
    up::HState{Float64}   # state at τ = 0 (outer end)
    par::Float64          # ψ (fast: ℓ = ψτ) or s_end (slow: s = s_end τ)
    sol::Any              # ODE solution, dense
end

"""
    Wave

One elementary wave, canonical frame, stored in physical left-to-right order:
`left` is the state on its left, `right` on its right. For a fan `s_left < s_right`.
`kind` ∈ (:fast_shock, :fast_fan, :rotation, :slow_shock, :slow_fan, :contact).
"""
struct Wave
    kind::Symbol
    side::Int
    s_left::Float64
    s_right::Float64
    left::HState{Float64}
    right::HState{Float64}
    fan::Union{Nothing,FanData}
end

"""
Symmetry map between the user frame and the canonical frame.
"""
struct Frame
    mirror::Bool            # x → -x combined with B → -B (sides swapped, Bn keeps its sign)
    sB::Int                 # sign flip of B
    θ::Float64              # rotation angle of the transverse plane
    vtR::SVector{2,Float64} # transverse velocity of R after flip/rotation (unscaled)
    ρ0::Float64
    p0::Float64
end

struct RiemannProblem
    L::SVector{8,Float64}   # (ρ, vx, vy, vz, Bx, By, Bz, p)
    R::SVector{8,Float64}
    γ::Float64
end
RiemannProblem(L, R; γ = 5 / 3) = RiemannProblem(SVector{8,Float64}(L), SVector{8,Float64}(R), Float64(γ))

"""
    RiemannSolution

Result of [`solve`](@ref). Fields: `retcode`, `reason` (Symbol), `prob`, `Ψ` (path
variables), `residual`, `waves` (canonical frame), `frame`, `ctx`, `check`
(report of the independent checks), `method` (how the solution was found).
"""
struct RiemannSolution
    retcode::RetCode
    reason::Symbol
    prob::RiemannProblem
    Ψ::SVector{5,Float64}
    residual::Float64
    waves::Vector{Wave}
    frame::Union{Nothing,Frame}
    ctx::Union{Nothing,Ctx}
    check::Any
    method::Symbol
end

function Base.show(io::IO, sol::RiemannSolution)
    print(io, "RiemannSolution(", sol.retcode, sol.reason === :none ? "" : " :" * String(sol.reason),
        ", ", length(sol.waves), " waves, residual = ", sol.residual, ")")
end

successful(sol::RiemannSolution) = sol.retcode == Success
