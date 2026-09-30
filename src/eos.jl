# Characteristic speeds, conserved variables and fluxes of 1D ideal MHD, γ-law gas.

"""
    speeds(ρ, p, bt, Bn, γ) -> (cf, cA, cs)

Fast, Alfvén and slow speeds. `cs` is evaluated as a·cA/cf, which avoids the
cancellation of the textbook formula.
"""
@inline function speeds(ρ, p, bt, Bn, γ)
    a2 = γ * p / ρ
    ca2 = Bn^2 / ρ
    b2 = ca2 + bt^2 / ρ
    d = sqrt(max((a2 + b2)^2 - 4 * a2 * ca2, zero(a2)))
    cf2 = (a2 + b2 + d) / 2
    cf = sqrt(cf2)
    cs = sqrt(a2 * ca2 / cf2)
    return cf, sqrt(ca2), cs
end
speeds(h::HState, ctx::Ctx) = speeds(h.ρ, h.p, h.bt, ctx.Bn, ctx.γ)

"""Speeds for a full primitive state W = (ρ, vx, vy, vz, Bx, By, Bz, p)."""
speeds(W::AbstractVector, γ) = speeds(W[1], W[8], hypot(W[6], W[7]), W[5], γ)

"""Conserved variables (ρ, ρvx, ρvy, ρvz, By, Bz, E) of a primitive state."""
function conserved(W, γ)
    ρ, u, v, w, bx, by, bz, p = W
    E = p / (γ - 1) + ρ * (u^2 + v^2 + w^2) / 2 + (bx^2 + by^2 + bz^2) / 2
    return SVector(ρ, ρ * u, ρ * v, ρ * w, by, bz, E)
end

"""Physical flux in x of a primitive state (Bx constant, its equation omitted)."""
function flux(W, γ)
    ρ, u, v, w, bx, by, bz, p = W
    B2 = bx^2 + by^2 + bz^2
    pt = p + B2 / 2
    E = p / (γ - 1) + ρ * (u^2 + v^2 + w^2) / 2 + B2 / 2
    vB = u * bx + v * by + w * bz
    return SVector(ρ * u, ρ * u^2 + pt - bx^2, ρ * u * v - bx * by, ρ * u * w - bx * bz,
        by * u - bx * v, bz * u - bx * w, (E + pt) * u - bx * vB)
end

"""Quasi-linear matrix of the primitive system q = (ρ, u, v, w, By, Bz, p)."""
function primitive_jacobian(W, γ)
    ρ, u, v, w, bx, by, bz, p = W
    A = zeros(7, 7)
    A[1, 1] = u; A[1, 2] = ρ
    A[2, 2] = u; A[2, 5] = by / ρ; A[2, 6] = bz / ρ; A[2, 7] = 1 / ρ
    A[3, 3] = u; A[3, 5] = -bx / ρ
    A[4, 4] = u; A[4, 6] = -bx / ρ
    A[5, 2] = by; A[5, 3] = -bx; A[5, 5] = u
    A[6, 2] = bz; A[6, 4] = -bx; A[6, 6] = u
    A[7, 2] = γ * p; A[7, 7] = u
    return A
end
