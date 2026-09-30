# Fast and slow rarefaction fans: integral curves of the fast/slow eigenvectors.
# Along a fan ρ = ρ0 e^{-s}, p = p0 e^{-γ s} (paper eq. 53); u, Bt and the transverse
# velocity follow ODEs (paper eq. 54, second line corrected: see docs).
#
# Fast fan: parametrized by ℓ = log(Bt/Bt0) = ψ τ, τ ∈ [0,1], ψ ≤ 0. Then
#   ds/dℓ = cA²/cf² - 1 ∈ (-1, 0),  du/dℓ = -σ cf ds/dℓ,  dw/dℓ = -σ Bt Bn/(ρ cf),
# which is smooth for every ψ ∈ (-∞, 0]: Bt → 0 (switch-off) is reached only as ψ → -∞.
# Slow fan: in terms of s,
#   du/ds = -σ cs,  dBt/ds = Bt cs²/(cA² - cs²),  dw/ds = -σ cs Bt cA²/(Bn (cA² - cs²)).
# When the upstream Bt is tiny and a > cA, dBt/ds ~ 1/Bt (switch-on start, Bt ~ √s), so
# s is a poor parameter there. The fan is therefore parametrized by the arc length
# ς of the curve (s, Bt/√p0):  dς² = ds² + (dBt/√p0)²,  ς = ς_end τ with
# ς_end = s_vac tanh(-ψ/s_vac). This is s for ordinary fans, Bt/√p0 near switch-on,
# and every right-hand side stays bounded.
# w is the transverse velocity change along the (fixed) direction of Bt.

function fast_fan_rhs(y, q, τ)
    ρ0, p0, bt0, ψ, γ, Bn, σ = q
    s, _, _ = y
    ρ = ρ0 * exp(-s)
    p = p0 * exp(-γ * s)
    bt = bt0 * exp(ψ * τ)
    cf, cA, _ = speeds(ρ, p, bt, Bn, γ)
    dsdl = cA^2 / cf^2 - 1
    return ψ * SVector(dsdl, -σ * cf * dsdl, -σ * bt * Bn / (ρ * cf))
end

function slow_fan_rhs(y, q, τ)
    ρ0, p0, send, γ, Bn, σ = q
    s, _, bt, _ = y
    ρ = ρ0 * exp(-s)
    p = p0 * exp(-γ * s)
    cf, cA, cs = speeds(ρ, p, bt, Bn, γ)
    if cA^2 < γ * p / ρ
        # a > cA: cs → cA as Bt → 0; cA² - cs² = cA² bt²/(ρ (cf² - cA²)) avoids the cancellation
        g = ρ * (cf^2 - cA^2)
        Bp = cs^2 * g / (cA^2 * bt)                  # dBt/ds
        Wp = -σ * cs * g / (Bn * bt)                 # dw/ds
    else
        den = cA^2 - cs^2                            # cs < a ≤ cA: well separated
        Bp = bt * cs^2 / den
        Wp = -σ * cs * bt * cA^2 / (Bn * den)
    end
    N = sqrt(1 + (Bp / sqrt(p0))^2)                  # dς/ds
    return send * SVector(1 / N, -σ * cs / N, Bp / N, Wp / N)
end

function _fan_solve(rhs, y0, q, ctx; dense::Bool)
    prob = ODEProblem{false}(rhs, y0, (0.0, 1.0), q)
    sol = OrdinaryDiffEqVerner.solve(prob, Vern9(); abstol = ctx.opts.ode_tol, reltol = ctx.opts.ode_tol,
        dense = dense, save_everystep = dense, maxiters = 100_000)
    SciMLBase.successful_retcode(sol) || throw(DomainError(fvalue(q[4]), "fan ODE failed: $(sol.retcode)"))
    return sol
end

"""
    fast_fan(U, ψ, σ, ctx; record=false) -> (downstream, FanData or nothing)
"""
function fast_fan(U::HState, ψ, σ::Int, ctx::Ctx; record::Bool = false)
    T = promote_type(typeof(U.ρ), typeof(ψ))
    q = (U.ρ, U.p, U.bt, ψ, ctx.γ, ctx.Bn, σ)
    y0 = SVector{3,T}(0, U.u, 0)
    sol = _fan_solve(fast_fan_rhs, y0, q, ctx; dense = record)
    s, u, w = sol.u[end]
    e = SVector(cos(U.φ), sin(U.φ))
    D = HState(U.ρ * exp(-s), u, U.p * exp(-ctx.γ * s), U.bt * exp(ψ), U.φ, U.vt + w * e)
    fan = record ? FanData(:fast, σ, fvalue(U), fvalue(ψ), sol) : nothing
    return D, fan
end

slow_fan_send(ψ, ctx) = ctx.opts.s_vac * tanh(-ψ / ctx.opts.s_vac)

"""
    slow_fan(U, ψ, σ, ctx; record=false) -> (downstream, FanData or nothing)
"""
function slow_fan(U::HState, ψ, σ::Int, ctx::Ctx; record::Bool = false)
    send = slow_fan_send(ψ, ctx)
    T = promote_type(typeof(U.ρ), typeof(send))
    q = (U.ρ, U.p, send, ctx.γ, ctx.Bn, σ)
    y0 = SVector{4,T}(0, U.u, U.bt, 0)
    sol = _fan_solve(slow_fan_rhs, y0, q, ctx; dense = record)
    s, u, bt, w = sol.u[end]
    e = SVector(cos(U.φ), sin(U.φ))
    D = HState(U.ρ * exp(-s), u, U.p * exp(-ctx.γ * s), bt, U.φ, U.vt + w * e)
    fan = record ? FanData(:slow, σ, fvalue(U), fvalue(send), sol) : nothing
    return D, fan
end

"""State inside a fan at parameter τ ∈ [0, 1]."""
function fan_state(f::FanData, τ, ctx::Ctx)
    U = f.up
    y = f.sol(τ)
    e = SVector(cos(U.φ), sin(U.φ))
    if f.family === :fast
        s, u, w = y
        return HState(U.ρ * exp(-s), u, U.p * exp(-ctx.γ * s), U.bt * exp(f.par * τ), U.φ, U.vt + w * e)
    else
        s, u, bt, w = y
        return HState(U.ρ * exp(-s), u, U.p * exp(-ctx.γ * s), bt, U.φ, U.vt + w * e)
    end
end

"""Characteristic speed of the fan family at τ."""
function fan_speed(f::FanData, τ, ctx::Ctx)
    h = fan_state(f, τ, ctx)
    cf, _, cs = speeds(h, ctx)
    return h.u + f.σ * (f.family === :fast ? cf : cs)
end
