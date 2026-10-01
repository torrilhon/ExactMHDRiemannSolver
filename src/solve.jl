# Driver: input checks, canonical frame, nonlinear solve with a small fixed set of
# starts and a homotopy fallback, limit detection and independent checks.

import SciMLBase: solve

const BIG = 1e6

# Residual that never throws: domain errors map to a large residual, so the trust
# region simply rejects the step.
function safe_residual(Ψ, UL, UR, ctx)
    try
        r = residual(Ψ, UL, UR, ctx)
        all(x -> isfinite(fvalue(x)), r) && return r
    catch err
        err isa DomainError || rethrow()
    end
    T = eltype(Ψ)
    return SVector{5,T}(BIG, BIG, BIG, BIG, BIG)
end

"""Natural scales (total pressure, fast speed, field strength) of the two input states."""
function problem_scale(UL, UR, Bn, γ)
    P = max(UL.p + (Bn^2 + UL.bt^2) / 2, UR.p + (Bn^2 + UR.bt^2) / 2)
    C = max(speeds(UL.ρ, UL.p, UL.bt, Bn, γ)[1], speeds(UR.ρ, UR.p, UR.bt, Bn, γ)[1])
    B = sqrt(Bn^2 + max(UL.bt, UR.bt)^2)
    return SVector(P, C, B)
end

resnorm(Ψ, UL, UR, ctx) = maximum(abs, safe_residual(SVector{5,Float64}(Ψ), UL, UR, ctx))

function nl_solve(Ψ0, UL, UR, ctx)
    f(u, p) = safe_residual(u, UL, UR, ctx)
    nlp = NonlinearProblem{false}(f, SVector{5,Float64}(Ψ0))
    s = try
        solve(nlp, TrustRegion(; autodiff = AutoForwardDiff()); abstol = ctx.opts.abstol,
            maxiters = ctx.opts.maxiters)
    catch err
        # numerical breakdown only; programming errors must surface
        err isa Union{InterruptException, OutOfMemoryError, StackOverflowError, ArgumentError, MethodError,
            UndefVarError, BoundsError, TypeError, UndefKeywordError} && rethrow()
        ctx.opts.verbose && @warn "nonlinear solve threw" err
        return SVector{5,Float64}(Ψ0), Inf
    end
    Ψ = SVector{5,Float64}(s.u)
    return Ψ, resnorm(Ψ, UL, UR, ctx)
end

accept(r, ctx) = r <= max(1e3 * ctx.opts.abstol, 1e-10)

"""Starting guesses in the canonical frame (α = twist angle of the right state)."""
function starts(α)
    g = [SVector(0.01, 0.01, α / 2, 0.01, 0.01)]
    if abs(α) > π - 1e-6          # coplanar: α/2 = ±π/2 is a poor start
        push!(g, SVector(0.01, 0.01, 0.0, 0.01, 0.01), SVector(0.01, 0.01, π, 0.01, 0.01))
    end
    push!(g, SVector(0.01, 0.01, α / 2 + π, 0.01, 0.01))
    # slow fans as start: on a side with tiny Bt the slow-shock branch is capped at
    # ΔB̂t ≤ A and its path variable saturates at once, so Newton needs ψs ≤ 0 there
    push!(g, SVector(0.01, -0.01, α / 2, -0.01, 0.01), SVector(0.01, 0.01, α / 2, -0.01, 0.01))
    return g
end

"""
Homotopy from the trivial problem (right state = left state) to the real one.
Returns (Ψ, residual, λ, Ψλ): on failure λ < 1 is the last value reached and Ψλ its
solution, which is used to diagnose why the continuation stalled.
"""
function homotopy(UL, UR, ctx)
    Rλ(λ) = HState(exp((1 - λ) * log(UL.ρ) + λ * log(UR.ρ)), (1 - λ) * UL.u + λ * UR.u,
        exp((1 - λ) * log(UL.p) + λ * log(UR.p)), (1 - λ) * UL.bt + λ * UR.bt,
        λ * UR.φ, (1 - λ) * UL.vt + λ * UR.vt)
    λ, h = 0.0, 0.05
    Ψ, Ψold, λold = zero(SVector{5,Float64}), zero(SVector{5,Float64}), 0.0
    t0 = time()
    for _ in 1:ctx.opts.homotopy_maxsteps
        time() - t0 > ctx.opts.time_limit && break
        λn = min(λ + h, 1.0)
        # linear predictor from the last two points
        pred = λ > 0 ? Ψ + (Ψ - Ψold) * ((λn - λ) / (λ - λold)) : Ψ + SVector(1e-4, 1e-4, 0.0, 1e-4, 1e-4)
        Ψn, r = nl_solve(pred, UL, Rλ(λn), ctx)
        if accept(r, ctx)
            Ψold, λold, Ψ, λ = Ψ, λ, Ψn, λn
            λ >= 1.0 && return Ψ, r, 1.0, Rλ(1.0)
            h = min(2h, 0.25)
        else
            h /= 2
            h < 1e-4 && break
        end
    end
    return Ψ, Inf, λ, Rλ(λ)
end

failed(prob, rc, reason; frame = nothing, ctx = nothing, Ψ = zero(SVector{5,Float64}), r = Inf, waves = Wave[], method = :none, chk = nothing) =
    RiemannSolution(rc, reason, prob, Ψ, r, waves, frame, ctx, chk, method)

# Regular-limit diagnostics of a converged Ψ.
function limit_reason(Ψ, waves, UL, UR, ctx)
    o = ctx.opts
    for w in waves
        for h in (w.left, w.right)
            (h.ρ < o.rho_floor || h.p < o.rho_floor) && return :vacuum
            if h.bt / sqrt(h.p) < o.bt_floor
                return any(x -> x.kind === :fast_fan, waves) ? :fast_switch_off : :bt_small
            end
        end
    end
    for (ψs, σ) in ((Ψ[2], -1), (Ψ[4], +1))
        if ψs > SLOW_EPS
            idx = findfirst(x -> x.kind === :rotation && x.side == σ, waves)
            Rt = σ < 0 ? waves[idx].right : waves[idx].left
            Δmax = (1 - o.δ) * slow_limit(Rt, ctx)
            ψs / Δmax > o.sat && return :slow_limit
        elseif -ψs / o.s_vac > o.sat
            return :vacuum
        end
    end
    return :none
end

"""
    solve(prob::RiemannProblem; opts = SolverOptions(), guess = nothing) -> RiemannSolution

Solve the Riemann problem. `guess` (5 path variables in the canonical frame) is
tried first if given. Never throws for physical reasons; inspect `sol.retcode`.
"""
function solve(prob::RiemannProblem; opts::SolverOptions = SolverOptions(), guess = nothing)
    opts.verbose && return _solve(prob, opts, guess)
    return Logging.with_logger(Logging.NullLogger()) do   # silence solver-internal warnings
        _solve(prob, opts, guess)
    end
end

function _solve(prob::RiemannProblem, opts::SolverOptions, guess)
    rc, reason = validate_input(prob, opts)
    rc == Success || return failed(prob, rc, reason)
    if reason === :quasi_euler
        return try
            solve_perpendicular(prob, opts)
        catch err
            err isa DomainError || rethrow()
            failed(prob, NoConvergence, :no_convergence; method = :quasi_euler)
        end
    end
    F = make_frame(prob.L, prob.R)
    Lc, Rc = canonical_states(prob.L, prob.R, F)
    UL, UR = to_hstate(Lc), to_hstate(Rc)
    ctx = Ctx(prob.γ, Lc[5], opts, problem_scale(UL, UR, Lc[5], prob.γ))
    if maximum(abs, Lc - Rc) <= 1e-14
        w = [Wave(:contact, 0, UL.u, UL.u, UL, UR, nothing)]
        return RiemannSolution(Success, :trivial, prob, zero(SVector{5,Float64}), 0.0, w, F, ctx,
            check_waves(w, ctx, unmirrored(F), prob.γ), :trivial)
    end
    α = UR.φ
    cand = guess === nothing ? starts(α) : vcat([SVector{5,Float64}(guess)], starts(α))
    Ψ, r, method = zero(SVector{5,Float64}), Inf, :none
    for Ψ0 in cand
        Ψ, r = nl_solve(Ψ0, UL, UR, ctx)
        if accept(r, ctx)
            method = :direct
            break
        end
    end
    if !accept(r, ctx) && opts.homotopy
        Ψ, r, λ, URλ = homotopy(UL, UR, ctx)
        method = :homotopy
        if !accept(r, ctx) && λ > 0
            # why did the continuation stall? A limit reached on the way is reported
            # as such (e.g. vacuum generation), otherwise plain non-convergence.
            lim = try
                limit_reason(Ψ, assemble(Ψ, UL, URλ, ctx), UL, URλ, ctx)
            catch err
                err isa DomainError || rethrow()
                :none
            end
            lim === :none || return failed(prob, RegularLimit, lim; frame = F, ctx, Ψ, r, method = :homotopy_stalled)
        end
    end
    accept(r, ctx) || return failed(prob, NoConvergence, :no_convergence; frame = F, ctx, Ψ, r, method)
    # recording the waves (dense fan output) and the checks can still break down
    # numerically in degenerate cases; that must never escape as an exception
    waves, chk, lim = try
        w = assemble(Ψ, UL, UR, ctx)
        w, check_waves(w, ctx, unmirrored(F), prob.γ), limit_reason(Ψ, w, UL, UR, ctx)
    catch err
        err isa DomainError || rethrow()
        return failed(prob, NoConvergence, :assembly_failed; frame = F, ctx, Ψ, r, method)
    end
    lim === :none || return RiemannSolution(RegularLimit, lim, prob, Ψ, r, waves, F, ctx, chk, method)
    chk.ok || return RiemannSolution(CheckFailed, :check, prob, Ψ, r, waves, F, ctx, chk, method)
    return RiemannSolution(Success, :none, prob, Ψ, r, waves, F, ctx, chk, method)
end
