# Fast and slow shocks (paper Sec. 4.2 and Appendices A/B), in local dimensionless
# variables v̂ = v1/v0, p̂ = p1/p0, B̂t = Bt/√p0, A = B̂t of the upstream state,
# B = Bn/√p0, X = γM² (M = shock speed relative to the gas, over the sound speed).

"""
    ift_polish(f, fv, x0)

Given a root `x0` of the value-level function `fv`, return x0 - f(x0)/fv'(x0), where
`f` is the same function with possibly dual-valued parameters. For Float64 this is
one Newton polish step; for ForwardDiff duals it yields the exact first derivative
of the root by the implicit function theorem.
"""
function ift_polish(f, fv, x0::Float64)
    d = ForwardDiff.derivative(fv, x0)
    return x0 - f(x0) / d
end

# Fast shock, formulated in the downstream transverse field y = B̂t1 (derived
# symbolically from paper eqs. 43). With D = X - B² ≥ 0 the compression is exactly
#   1 - v̂ = D (y - A) / (X y),
# and after dividing out the trivial root y = A the Hugoniot is the cubic
#   c3 y³ + c2 y² + c1 y + c0 = 0,
#   c3 = -B²(k-1),  c2 = A(2D - X(k-1)),
#   c1 = -D(A²(k+1) + 2Dk - 2X(k-1) + 2(k+1)),  c0 = 2 A D² k,
# whose largest real root in the physical range 0 < y ≤ √(A² + 2(1+X)) (p̂ ≥ 0) is the
# fast shock. Unlike the cubic in v̂ (paper eq. 73) this
# stays well conditioned near switch-on (A → 0 with c_A > a): there the root tends to
# the switch-on value of eq. (75) instead of approaching the pole v̂ = B²/X.
@inline fast_cubic_bt(A, B2, X, D, k) =
    (2 * A * D^2 * k, -D * (A^2 * (k + 1) + 2 * D * k - 2 * X * (k - 1) + 2 * (k + 1)),
     A * (2 * D - X * (k - 1)), -B2 * (k - 1))

# Largest real root of a real cubic in the interval (lo, hi], by bracketing on its
# monotone pieces (no complex arithmetic). The interval matters when Bn → 0: then the
# leading coefficient c3 ∝ B² vanishes and a spurious root of size ~1/B² appears far
# outside the physical range.
function largest_root_in(c::NTuple{4,Float64}, lo::Float64, hi::Float64)
    c0, c1, c2, c3 = c
    f(x) = evalpoly(x, c)
    disc = c2^2 - 3 * c3 * c1                      # critical points of the cubic
    crit = disc > 0 ? [(-c2 - sqrt(disc)) / (3c3), (-c2 + sqrt(disc)) / (3c3)] : Float64[]
    pts = sort([lo; filter(x -> lo < x < hi, crit); hi])
    for i in length(pts)-1:-1:1
        a, b = pts[i], pts[i+1]
        fa, fb = f(a), f(b)
        fb == 0 && return b
        fa * fb < 0 && return find_zero(f, (a, b), Roots.Brent(); xatol = 0.0, xrtol = 4eps())
    end
    throw(DomainError(c, "fast shock: no root of the Hugoniot cubic in the physical range"))
end

"""
    fast_shock(U, ψ, σ, ctx) -> (downstream, speed)

Fast shock on the state `U` travelling in direction σ (-1 left, +1 right), with path
variable ψ = √(γ(M² - ĉA²)) - √(γ(ĉf² - ĉA²)) ≥ 0 (M: shock Mach number).
"""
function fast_shock(U::HState, ψ, σ::Int, ctx::Ctx)
    γ, κ, Bn = ctx.γ, ctx.κ, ctx.Bn
    sp = sqrt(U.p)
    A, B = U.bt / sp, Bn / sp
    a0 = sqrt(γ * U.p / U.ρ)
    cf, cA, cs = speeds(U, ctx)
    # c_f² - c_A² without cancellation: (cf² - cA²)(cA² - cs²) = cA² bt²/ρ
    df2 = cA > sqrt(γ * U.p / U.ρ) ? cA^2 * U.bt^2 / (U.ρ * (cA^2 - cs^2)) : cf^2 - cA^2
    # Path variable ψ = √D - √D_f with D = X - B² = γ(M² - ĉA²) and D_f its value at
    # M = ĉ_f. For ordinary states this is ~ linear in M; near switch-on (tiny A, cA > a)
    # the downstream field grows like √D, so it grows ~ linearly in ψ instead of like √ψ.
    Df = γ * df2 / a0^2
    D = (sqrt(Df) + ψ)^2
    X = B^2 + D
    M = sqrt(X / γ)
    c = fast_cubic_bt(A, B^2, X, D, κ)
    # physical range of the downstream field: 0 < y and p̂ ≥ 0, i.e. y² ≤ A² + 2(1 + X)
    yhi = sqrt(fvalue(A)^2 + 2 * (1 + fvalue(X))) * (1 + 1e-12)
    y0 = largest_root_in(fvalue.(c), 1e-300, yhi)
    y = ift_polish(x -> evalpoly(x, c), x -> evalpoly(x, fvalue.(c)), y0)
    t = D * (y - A) / (X * y)                    # 1 - v̂
    ph = 1 + X * t - (y^2 - A^2) / 2
    bt1 = y * sp
    Cc = -σ * Bn / (U.ρ * a0 * M)
    e = SVector(cos(U.φ), sin(U.φ))
    Dn = HState(U.ρ / (1 - t), U.u + σ * a0 * M * t, U.p * ph, bt1, U.φ, U.vt + Cc * (bt1 - U.bt) * e)
    return Dn, fvalue(U.u + σ * a0 * M)
end

# Slow shock lowering B̂t from A to A - Δ (paper eqs. 78-83). The quadratic (79) has
# the double root v̂ = 1 at Δ = 0, so it is rewritten for t = 1 - v̂ = Δ z:
#   a z² - e1 z + e0 = 0,  e1 = (2a+b)/Δ,  e0 = (a+b+c)/Δ²  (derived symbolically),
# which is well conditioned for weak shocks. The compressive root is z > 0 (branch v̂⁺).
# Returns (v̂, X = γM²) with X = B²/(1 + (A-Δ) z), also free of cancellation.
@inline function slow_volume(Δ, A, B2, k)
    D = Δ
    a = A^2 * D * k - 3 * A * D^2 * k / 2 + A * D^2 / 2 + A * k + A + B2 * D * k + D^3 * k / 2 - D^3 / 2 - D * k - D
    e1 = (2 * A^2 * k - 2 * A^2 - 5 * A * D * k + 3 * A * D + 2 * B2 * k - 2 * B2 + 2 * D^2 * k - 2 * D^2 - 2 * k - 2) / 2
    e0 = -(2 * A - D) * (k - 1) / 2
    disc = e1^2 - 4 * a * e0
    (fvalue(a) > 0 && fvalue(e0) < 0) || throw(DomainError(fvalue(Δ), "slow shock outside the regular branch"))
    z = fvalue(e1) > 0 ? (e1 + sqrt(disc)) / (2 * a) : -2 * e0 / (sqrt(disc) - e1)
    return 1 - D * z, B2 / (1 + (A - D) * z)
end

@inline slow_pH(v, Bth, A, κ) = (v - κ + (Bth - A)^2 * (v - 1) / 2) / (1 - κ * v)

"""
    slow_shock_state(U, Δ, σ, ctx) -> (downstream, speed, M, v̂)

Slow shock that lowers B̂t by Δ > 0.
"""
function slow_shock_state(U::HState, Δ, σ::Int, ctx::Ctx)
    γ, κ, Bn = ctx.γ, ctx.κ, ctx.Bn
    sp = sqrt(U.p)
    A, B2 = U.bt / sp, (Bn / sp)^2
    a0 = sqrt(γ * U.p / U.ρ)
    Bth = A - Δ
    v, X = slow_volume(Δ, A, B2, κ)
    M = sqrt(X / γ)
    # momentum balance (Rayleigh line) instead of the Hugoniot form slow_pH: the latter
    # divides by 1 - κ v̂, which vanishes when Bn → 0 makes the slow shock maximally compressive
    ph = 1 - X * (v - 1) - (Bth^2 - A^2) / 2
    bt1 = Bth * sp
    Cc = -σ * Bn / (U.ρ * a0 * M)
    e = SVector(cos(U.φ), sin(U.φ))
    D = HState(U.ρ / v, U.u + σ * a0 * M * (1 - v), U.p * ph, bt1, U.φ, U.vt + Cc * (bt1 - U.bt) * e)
    return D, fvalue(U.u + σ * a0 * M), M, v
end

# Lax defect of a slow shock: (downstream relative speed)² - (downstream slow speed)².
# Negative for a regular 3→4 shock, zero at the maximal Mach number.
function slow_defect(U::HState, Δ, ctx::Ctx)
    D, _, M, v = slow_shock_state(U, Δ, 1, ctx)
    w1 = sqrt(ctx.γ * U.p / U.ρ) * M * v
    _, _, cs1 = speeds(D, ctx)
    return (w1^2 - cs1^2) / (ctx.γ * U.p / U.ρ)
end

"""
    slow_limit(U, ctx) -> Δmax

Largest drop Δ of B̂t for which the slow shock on `U` stays a regular 3→4 shock with
B̂t > 0: either the point of maximal Mach number (paper eqs. 86-88) or Δ = A.
Found by a coarse scan plus a bracketing root, replacing the magic start value of
the C code.
"""
function slow_limit(U::HState, ctx::Ctx; N::Int = 64)
    Uv = fvalue(U)
    A = Uv.bt / sqrt(Uv.p)
    fv(Δ) = slow_defect(Uv, Δ, ctx)
    Δprev = 1e-4 * A
    fprev = fv(Δprev)
    fprev >= 0 && return oftype(U.bt / sqrt(U.p), Δprev)   # degenerate: no regular slow shock
    for i in 1:N
        Δi = i == N ? A * (1 - 1e-9) : A * i / N
        fi = try
            fv(Δi)
        catch err
            err isa DomainError || rethrow()
            NaN
        end
        if isnan(fi)             # quadratic ceased to exist: end of the branch
            return oftype(U.bt / sqrt(U.p), Δprev)
        elseif fi >= 0
            Δ0 = find_zero(fv, (Δprev, Δi), Roots.Brent(); xatol = 0.0, xrtol = 4eps())
            return ift_polish(Δ -> slow_defect(U, Δ, ctx), fv, Δ0)
        end
        Δprev, fprev = Δi, fi
    end
    return U.bt / sqrt(U.p)      # regular up to B̂t = 0
end
