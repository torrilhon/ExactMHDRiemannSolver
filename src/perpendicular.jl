# Quasi-Euler solver for a vanishing normal field (Bn = 0, used for Bn/√p ≤ bn_euler).
#
# With Bn = 0 the slow and Alfvén waves merge into the contact, which becomes a
# tangential discontinuity: u and the total pressure P = p + Bt²/2 are continuous,
# while ρ, p, |Bt|, the direction of Bt and vt may jump. Across the fast waves Bt/ρ,
# the direction of Bt and vt stay constant, and the fast speed is cf² = (γp + Bt²)/ρ.
# This is gas dynamics with an extra "magnetic pressure" ∝ ρ²; it is strictly
# hyperbolic and genuinely nonlinear for every Bt ≥ 0, so no Bt thresholds apply.
# One unknown, P*, is found from u*_L(P*) = u*_R(P*) by bracketing.

ptot(ρ, p, bt) = p + bt^2 / 2

# Fast fan on U down to total pressure P < P0, parametrized by s = log(ρ0/ρ) = s* τ.
function perp_fan_rhs(y, q, τ)
    ρ0, p0, bt0, send, γ, σ = q
    s = send * τ
    ρ = ρ0 * exp(-s); p = p0 * exp(-γ * s); bt = bt0 * exp(-s)
    cf = sqrt((γ * p + bt^2) / ρ)
    return SVector(send * (-σ * cf))
end

function perp_fan(U::HState, P, σ::Int, ctx::Ctx; record::Bool = false)
    γ = ctx.γ
    g(ρ) = ptot(ρ, U.p * (ρ / U.ρ)^γ, U.bt * ρ / U.ρ) - P
    ρs = find_zero(g, (1e-300, U.ρ), Roots.Brent(); xatol = 0.0, xrtol = 4eps())
    send = log(U.ρ / ρs)
    prob = ODEProblem{false}(perp_fan_rhs, SVector(U.u), (0.0, 1.0), (U.ρ, U.p, U.bt, send, γ, σ))
    sol = OrdinaryDiffEqVerner.solve(prob, Vern9(); abstol = ctx.opts.ode_tol, reltol = ctx.opts.ode_tol,
        dense = record, save_everystep = record)
    SciMLBase.successful_retcode(sol) || throw(DomainError(P, "quasi-Euler fan ODE failed"))
    D = HState(ρs, sol.u[end][1], U.p * (ρs / U.ρ)^γ, U.bt * ρs / U.ρ, U.φ, U.vt)
    return D, (record ? FanData(:fast0, σ, U, send, sol) : nothing)
end

# Fast shock on U up to total pressure P > P0: compression r = ρ1/ρ0 ∈ (1, κ) from the
# energy jump condition with Bt ∝ ρ, then the mass flux from the momentum condition.
function perp_shock(U::HState, P, σ::Int, ctx::Ctx)
    γ, κ = ctx.γ, ctx.κ
    A = U.bt / sqrt(U.p)
    function ph(r)
        v = 1 / r
        return ((1 - v) * (A * r - A)^2 / 4 + 1 / (γ - 1) - (v - 1) / 2) / (v / (γ - 1) + (v - 1) / 2)
    end
    f(r) = U.p * ph(r) + (U.bt * r)^2 / 2 - P
    rmax = κ * (1 - 1e-14)
    # beyond P/P0 ~ 1e14 the compression is indistinguishable from κ in Float64
    f(rmax) > 0 || throw(DomainError(P, "quasi-Euler shock too strong"))
    r = find_zero(f, (1.0, rmax), Roots.Brent(); xatol = 0.0, xrtol = 4eps())
    ρ1 = U.ρ * r
    P0 = ptot(U.ρ, U.p, U.bt)
    m = sqrt((P - P0) / (1 / U.ρ - 1 / ρ1))            # mass flux
    D = HState(ρ1, U.u + σ * (P - P0) / m, U.p * ph(r), U.bt * r, U.φ, U.vt)
    return D, U.u + σ * m / U.ρ
end

function perp_side(U::HState, P, σ::Int, ctx::Ctx; record::Bool = false)
    P0 = ptot(U.ρ, U.p, U.bt)
    if P > P0 * (1 + 1e-14)
        D, s = perp_shock(U, P, σ, ctx)
        w = record ? (σ < 0 ? Wave(:fast_shock, σ, s, s, U, D, nothing) : Wave(:fast_shock, σ, s, s, D, U, nothing)) : nothing
        return D, w
    elseif P < P0 * (1 - 1e-14)
        D, fan = perp_fan(U, P, σ, ctx; record)
        if record
            c0 = sqrt((ctx.γ * U.p + U.bt^2) / U.ρ); c1 = sqrt((ctx.γ * D.p + D.bt^2) / D.ρ)
            w = σ < 0 ? Wave(:fast_fan, σ, U.u - c0, D.u - c1, U, D, fan) : Wave(:fast_fan, σ, D.u + c1, U.u + c0, D, U, fan)
            return D, w
        end
        return D, nothing
    else
        return U, nothing        # no fast wave on this side
    end
end

"""
    solve_perpendicular(prob, opts) -> RiemannSolution

Quasi-Euler solution for Bn ≈ 0 (fast waves and a tangential discontinuity), in the
user's units. Refuses only vacuum generation.
"""
function solve_perpendicular(prob::RiemannProblem, opts::SolverOptions)
    L, R, γ = prob.L, prob.R, prob.γ
    F = Frame(false, 1, 0.0, SVector(0.0, 0.0), 1.0, 1.0)           # identity map
    ctx = Ctx(γ, L[5], opts)
    UL, UR = to_hstate(L), to_hstate(R)
    PL, PR = ptot(UL.ρ, UL.p, UL.bt), ptot(UR.ρ, UR.p, UR.bt)
    g(P) = perp_side(UL, P, -1, ctx)[1].u - perp_side(UR, P, +1, ctx)[1].u
    lo, hi = 1e-12 * min(PL, PR), 2 * max(PL, PR)
    g(lo) < 0 && return failed(prob, RegularLimit, :vacuum; frame = F, ctx, method = :quasi_euler)
    k = 0
    while g(hi) > 0
        hi *= 2; k += 1
        k > 200 && return failed(prob, NoConvergence, :no_convergence; frame = F, ctx, method = :quasi_euler)
    end
    Ps = find_zero(g, (lo, hi), Roots.Brent(); xatol = 0.0, xrtol = 4eps())
    SL, wl = perp_side(UL, Ps, -1, ctx; record = true)
    SR, wr = perp_side(UR, Ps, +1, ctx; record = true)
    us = (SL.u + SR.u) / 2
    waves = Wave[]
    wl === nothing || push!(waves, wl)
    push!(waves, Wave(:tangential, 0, us, us, SL, SR, nothing))
    wr === nothing || push!(waves, wr)
    r = abs(SL.u - SR.u) / max(sqrt((γ * max(L[8], R[8]) + max(UL.bt, UR.bt)^2) / min(L[1], R[1])), eps())
    chk = check_waves(waves, ctx, F, γ)
    lim = min(SL.ρ, SR.ρ, SL.p, SR.p) < opts.rho_floor * min(L[1], R[1], L[8], R[8]) ? :vacuum : :none
    Ψ = SVector(Ps, 0.0, 0.0, 0.0, 0.0)
    lim === :none || return RiemannSolution(RegularLimit, lim, prob, Ψ, r, waves, F, ctx, chk, :quasi_euler)
    chk.ok || return RiemannSolution(CheckFailed, :check, prob, Ψ, r, waves, F, ctx, chk, :quasi_euler)
    return RiemannSolution(Success, :none, prob, Ψ, r, waves, F, ctx, chk, :quasi_euler)
end
