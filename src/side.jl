# The single wave-sequence builder. Residual, output and checks all use it, so what
# Newton solves and what is reported cannot differ (the root cause of the false
# intermediate-shock solutions of the C code).

"""Alfvén (rotational) discontinuity turning Bt from U.φ to α (requires Bn > 0)."""
function rotation(U::HState, α, σ::Int, ctx::Ctx)
    e0 = SVector(cos(U.φ), sin(U.φ))
    e1 = SVector(cos(α), sin(α))
    dvt = -σ * U.bt * (e1 - e0) / sqrt(U.ρ)
    D = HState(U.ρ, U.u, U.p, U.bt, α, U.vt + dvt)
    _, cA, _ = speeds(U, ctx)
    return D, fvalue(U.u + σ * cA)
end

const SLOW_EPS = 1e-9   # slow ψ below this is treated as a (vanishing) fan

# Store a wave in physical left-right order.
function _push!(waves, kind, σ, up, down, s_up, s_down, fan = nothing)
    up, down = fvalue(up), fvalue(down)
    if σ < 0
        push!(waves, Wave(kind, σ, s_up, s_down, up, down, fan))
    else
        push!(waves, Wave(kind, σ, s_down, s_up, down, up, fan))
    end
end

"""
    build_side(U, ψf, ψs, αR, σ, ctx; record=false) -> (inner state, waves)

Apply fast wave, rotation to angle αR and slow wave to the outer state `U` of side
σ (-1 left, +1 right). With `record=true` the waves are returned in outer-to-inner
order, otherwise `waves` is `nothing`.
"""
function build_side(U::HState, ψf, ψs, αR, σ::Int, ctx::Ctx; record::Bool = false)
    waves = record ? Wave[] : nothing
    # fast wave
    if ψf > 0
        F, s = fast_shock(U, ψf, σ, ctx)
        record && _push!(waves, :fast_shock, σ, U, F, s, s)
    else
        F, fan = fast_fan(U, ψf, σ, ctx; record)
        if record
            cf0, _, _ = speeds(fvalue(U), ctx); cf1, _, _ = speeds(fvalue(F), ctx)
            _push!(waves, :fast_fan, σ, U, F, fvalue(U.u) + σ * cf0, fvalue(F.u) + σ * cf1, fan)
        end
    end
    # rotation
    Rt, s = rotation(F, αR, σ, ctx)
    record && _push!(waves, :rotation, σ, F, Rt, s, s)
    # slow wave
    if ψs > SLOW_EPS
        Δmax = (1 - ctx.opts.δ) * slow_limit(Rt, ctx)
        Δ = Δmax * tanh(ψs / Δmax)
        S, s, _, _ = slow_shock_state(Rt, Δ, σ, ctx)
        record && _push!(waves, :slow_shock, σ, Rt, S, s, s)
    else
        S, fan = slow_fan(Rt, ψs, σ, ctx; record)
        if record
            _, _, cs0 = speeds(fvalue(Rt), ctx); _, _, cs1 = speeds(fvalue(S), ctx)
            _push!(waves, :slow_fan, σ, Rt, S, fvalue(Rt.u) + σ * cs0, fvalue(S.u) + σ * cs1, fan)
        end
    end
    return S, waves
end

"""
    residual(Ψ, UL, UR, α, ctx)

Mismatch of (p, u, |Bt|, vt) between the inner states of the two sides, in units of
the problem's own scales (`ctx.scale`: total pressure, fast speed, field strength),
so that the tolerance means the same relative accuracy for every problem.
Ψ = (ψf⁻, ψs⁻, αR, ψs⁺, ψf⁺); in the canonical frame the left rotation turns 0 → αR
and the right one α → αR.
"""
function residual(Ψ, UL::HState, UR::HState, ctx::Ctx)
    SL, _ = build_side(UL, Ψ[1], Ψ[2], Ψ[3], -1, ctx)
    SR, _ = build_side(UR, Ψ[5], Ψ[4], Ψ[3], +1, ctx)
    P, C, B = ctx.scale
    return SVector((SL.p - SR.p) / P, (SL.u - SR.u) / C, (SL.bt - SR.bt) / B,
        (SL.vt[1] - SR.vt[1]) / C, (SL.vt[2] - SR.vt[2]) / C)
end

"""Assemble all waves of a solution, left to right, including the contact."""
function assemble(Ψ, UL::HState, UR::HState, ctx::Ctx)
    SL, wl = build_side(UL, Ψ[1], Ψ[2], Ψ[3], -1, ctx; record = true)
    SR, wr = build_side(UR, Ψ[5], Ψ[4], Ψ[3], +1, ctx; record = true)
    uc = (SL.u + SR.u) / 2
    waves = vcat(wl, [Wave(:contact, 0, uc, uc, SL, SR, nothing)], reverse(wr))
    return waves
end
