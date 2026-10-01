# Canonical frame: larger Bt/√p on the left, Bn > 0, left Bt along +y, right vt = 0,
# units ρ_L = p_L = 1. Every step is an exact symmetry of ideal MHD, so the map is
# invertible.
#
# Side swap: reflection x → -x (swap L and R, vx → -vx, Bx → -Bx) combined with field
# reversal B → -B gives W'(ξ) = T(W(-ξ)) with T(ρ, vx, vy, vz, Bx, By, Bz, p) =
# (ρ, -vx, vy, vz, Bx, -By, -Bz, p): Bn keeps its sign and the sides trade places. It is
# used to put the state with the smaller Bt/√p on the right, where the homotopy start
# (right state = left state) stays non-degenerate.

rot2(θ) = SMatrix{2,2}(cos(θ), sin(θ), -sin(θ), cos(θ))

mirrorT(W) = SVector(W[1], -W[2], W[3], W[4], W[5], -W[6], -W[7], W[8])

btp(W) = hypot(W[6], W[7]) / sqrt(W[8])

function make_frame(L0, R0)
    mirror = btp(L0) < btp(R0)
    L, R = mirror ? (mirrorT(R0), mirrorT(L0)) : (L0, R0)
    sB = L[5] < 0 ? -1 : 1
    θ = atan(sB * L[7], sB * L[6])
    Rm = rot2(-θ)
    vtR = Rm * SVector(R[3], R[4])
    return Frame(mirror, sB, θ, vtR, L[1], L[8])
end

"""Canonical 8-vectors (left, right) of a problem in the frame F."""
function canonical_states(L, R, F::Frame)
    F.mirror && ((L, R) = (mirrorT(R), mirrorT(L)))
    return _to_canonical(L, F), _to_canonical(R, F)
end

"""The same frame without the side swap: the physical frame of the (possibly mirrored)
problem, used by the independent checks."""
unmirrored(F::Frame) = Frame(false, F.sB, F.θ, F.vtR, F.ρ0, F.p0)

"""Primitive 8-vector of the (already side-swapped) problem → canonical 8-vector."""
function _to_canonical(W, F::Frame)
    Rm = rot2(-F.θ)
    vt = Rm * SVector(W[3], W[4]) - F.vtR
    bt = Rm * (F.sB * SVector(W[6], W[7]))
    c0 = sqrt(F.p0 / F.ρ0); b0 = sqrt(F.p0)
    return SVector(W[1] / F.ρ0, W[2] / c0, vt[1] / c0, vt[2] / c0, F.sB * W[5] / b0,
        bt[1] / b0, bt[2] / b0, W[8] / F.p0)
end

"""Canonical 8-vector → user-frame primitive 8-vector."""
function from_canonical(W, F::Frame)
    Rm = rot2(F.θ)
    c0 = sqrt(F.p0 / F.ρ0); b0 = sqrt(F.p0)
    vt = Rm * (SVector(W[3], W[4]) * c0 + F.vtR)
    bt = F.sB * (Rm * (SVector(W[6], W[7]) * b0))
    U = SVector(W[1] * F.ρ0, W[2] * c0, vt[1], vt[2], F.sB * W[5] * b0, bt[1], bt[2], W[8] * F.p0)
    return F.mirror ? mirrorT(U) : U
end

speed_to_user(s, F::Frame) = (F.mirror ? -s : s) * sqrt(F.p0 / F.ρ0)

"""HState (canonical) → canonical primitive 8-vector."""
to_prim(h::HState, ctx::Ctx) = SVector(h.ρ, h.u, h.vt[1], h.vt[2], ctx.Bn, h.bt * cos(h.φ), h.bt * sin(h.φ), h.p)

"""Canonical primitive 8-vector → HState (angle in (-π, π])."""
to_hstate(W) = HState(W[1], W[2], W[8], hypot(W[6], W[7]), atan(W[7], W[6]), SVector(W[3], W[4]))

"""
    validate_input(prob, opts) -> (retcode, reason)

Checks of the v1 input domain, applied before any solve.
"""
function validate_input(prob::RiemannProblem, opts::SolverOptions)
    L, R, γ = prob.L, prob.R, prob.γ
    (all(isfinite, L) && all(isfinite, R) && isfinite(γ)) || return (InvalidInput, :nonfinite)
    γ > 1 || return (InvalidInput, :gamma)
    (L[1] > 0 && R[1] > 0) || return (InvalidInput, :density)
    (L[8] > 0 && R[8] > 0) || return (InvalidInput, :pressure)
    abs(L[5] - R[5]) <= 1e-12 * max(abs(L[5]), 1.0) || return (InvalidInput, :bn_jump)
    rr, pr = R[1] / L[1], R[8] / L[8]
    (1 / opts.ratio_max <= rr <= opts.ratio_max && 1 / opts.ratio_max <= pr <= opts.ratio_max) ||
        return (Unsupported, :extreme_ratio)
    # vanishing normal field: quasi-Euler solver, any Bt (including 0) is fine
    max(abs(L[5]) / sqrt(L[8]), abs(R[5]) / sqrt(R[8])) <= opts.bn_euler && return (Success, :quasi_euler)
    # transverse field: one side may be very small as long as the other is not
    btL, btR = btp(L), btp(R)
    tolr = 1 - 1e-12        # a value set exactly at a threshold must not be refused by round-off
    (min(btL, btR) >= tolr * opts.bt_min_smaller && max(btL, btR) >= tolr * opts.bt_min_larger) ||
        return (Unsupported, :switch_on_off)
    return (Success, :none)
end
