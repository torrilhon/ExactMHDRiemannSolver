# Independent checks of an assembled solution, in the user frame, using the full
# conservative flux and a separate eigenvector integration for fans. They share no
# code with the wave kernels except the characteristic speeds.

user_state(h::HState, ctx::Ctx, F::Frame) = from_canonical(to_prim(h, ctx), F)

# Region index of a relative flow speed w: 1 super-fast … 4 sub-slow.
function lax_region(w, cf, cA, cs; tol = 1e-9)
    w > cf * (1 + tol) && return 1
    w > cA * (1 + tol) && return 2
    w > cs * (1 + tol) && return 3
    return 4
end

function fan_eigen_rhs(q, par, ρ)
    bx, γ, fam = par
    W = SVector(q[1], q[2], q[3], q[4], bx, q[5], q[6], q[7])
    E = eigen(primitive_jacobian(W, γ))
    idx = sortperm(real.(E.values))[fam]
    r = real.(E.vectors[:, idx])
    return SVector{7}(r ./ r[1])
end

"""Integrate the fan family `fam` (1, 3, 5 or 7) from state W0 to density ρ1."""
function fan_eigen_integrate(W0, ρ1, fam, γ)
    q0 = SVector(W0[1], W0[2], W0[3], W0[4], W0[6], W0[7], W0[8])
    prob = ODEProblem{false}(fan_eigen_rhs, q0, (W0[1], ρ1), (W0[5], γ, fam))
    s = OrdinaryDiffEqVerner.solve(prob, Vern7(); abstol = 1e-11, reltol = 1e-11, save_everystep = false)
    q = s.u[end]
    return SVector(q[1], q[2], q[3], q[4], W0[5], q[5], q[6], q[7])
end

# distance of a fan's family speed to the nearest other characteristic speed, relative
# to the fast speed (the scale of the Jacobian): numerical eigenvectors are only
# reliable when this is not small. For the slow family both the Alfvén speed and the
# contact (cs → 0 as Bn → 0) count.
function speed_gap(W, γ, kind)
    cf, cA, cs = speeds(W, γ)
    return kind === :fast_fan ? (cf - cA) / cf : min(cA - cs, cs) / cf
end

prim7(W) = SVector(W[1], W[2], W[3], W[4], W[6], W[7], W[8])

"""Largest relative defect of the eigen relation (A - λI) dq/dτ = 0 inside a fan."""
function fan_pointwise_defect(f::FanData, ctx::Ctx, F::Frame, γ)
    d = 0.0
    h = 1e-3
    q(t) = prim7(user_state(fan_state(f, t, ctx), ctx, F))
    for τ in (0.02, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.98)
        W = user_state(fan_state(f, τ, ctx), ctx, F)
        dq = (-q(τ + 2h) + 8q(τ + h) - 8q(τ - h) + q(τ - 2h)) / 12h      # 4th-order difference
        λ = speed_to_user(fan_speed(f, τ, ctx), F)
        A = primitive_jacobian(W, γ)
        nq = maximum(abs, dq)
        nq > 0 || continue
        # dq of a vanishingly weak fan is as small as the noise of the dense output divided
        # by h; measuring against a floor of 1e-6 of the state size keeps that noise from
        # being reported as a mismatch (the defect of such a fan is bounded by its strength)
        nfloor = 1e-6 * max(maximum(abs, prim7(W)), 1.0)
        d = max(d, maximum(abs, (A - λ * I) * dq) / (opnorm(A, Inf) * max(nq, nfloor)))
    end
    return d
end

relerr(a, b) = maximum(abs.(a .- b) ./ max.(abs.(a), abs.(b), 1.0))

"""
    check_waves(waves, ctx, F, γ) -> NamedTuple

Returns `(ok, maxerr, messages)`.
"""
function check_waves(waves::Vector{Wave}, ctx::Ctx, F::Frame, γ)
    tol = ctx.opts.check_tol
    msgs = String[]
    maxerr = 0.0
    note(ok, m) = ok || push!(msgs, m)
    for (i, wv) in enumerate(waves)
        Wl, Wr = user_state(wv.left, ctx, F), user_state(wv.right, ctx, F)
        sl, sr = speed_to_user(wv.s_left, F), speed_to_user(wv.s_right, F)
        σ = wv.side
        Wup, Wdn = σ < 0 ? (Wl, Wr) : (Wr, Wl)
        if wv.kind in (:fast_shock, :slow_shock, :rotation, :contact, :tangential)
            s = sl
            rh = s * (conserved(Wr, γ) - conserved(Wl, γ)) - (flux(Wr, γ) - flux(Wl, γ))
            scale = max(maximum(abs, flux(Wl, γ)), maximum(abs, flux(Wr, γ)), abs(s) * maximum(abs, conserved(Wl, γ)), 1.0)
            e = maximum(abs, rh) / scale
            maxerr = max(maxerr, e)
            note(e <= tol, "wave $i ($(wv.kind)): Rankine-Hugoniot error $e")
        end
        weak = abs(Wdn[1] - Wup[1]) <= 1e-7 * Wup[1]   # vanishing wave: only RH applies
        if wv.kind in (:fast_shock, :slow_shock) && !weak
            s = sl
            cfu, cAu, csu = speeds(Wup, γ); cfd, cAd, csd = speeds(Wdn, γ)
            wu = σ * (s - Wup[2]); wd = σ * (s - Wdn[2])
            ru, rd = lax_region(wu, cfu, cAu, csu), lax_region(wd, cfd, cAd, csd)
            want = wv.kind === :fast_shock ? (1, 2) : (3, 4)
            note((ru, rd) == want, "wave $i ($(wv.kind)): shock type $ru→$rd, expected $(want[1])→$(want[2])")
            note(Wdn[1] > Wup[1], "wave $i ($(wv.kind)): density does not increase (entropy)")
        elseif wv.kind === :rotation
            e = relerr(SVector(Wl[1], Wl[2], Wl[8], hypot(Wl[6], Wl[7])), SVector(Wr[1], Wr[2], Wr[8], hypot(Wr[6], Wr[7])))
            _, cA, _ = speeds(Wl, γ)
            e = max(e, abs(sl - (Wl[2] + σ * cA)) / max(abs(sl), 1.0))
            maxerr = max(maxerr, e)
            note(e <= tol, "wave $i (rotation): not a rotational discontinuity ($e)")
        elseif wv.kind === :tangential
            # Bn = 0: u and total pressure continuous, everything else may jump
            Pl = Wl[8] + (Wl[6]^2 + Wl[7]^2) / 2; Pr = Wr[8] + (Wr[6]^2 + Wr[7]^2) / 2
            e = max(abs(Wl[2] - Wr[2]) / max(1.0, abs(Wl[2])), abs(Pl - Pr) / max(Pl, Pr))
            maxerr = max(maxerr, e)
            note(e <= tol, "tangential discontinuity: jump in u or total pressure of $e")
            note(abs(sl - Wl[2]) <= tol * max(1.0, abs(sl)), "tangential discontinuity: speed differs from flow speed")
        elseif wv.kind === :contact
            e = relerr(Wl[2:8], Wr[2:8])
            maxerr = max(maxerr, e)
            note(e <= tol, "contact: jump in (v, B, p) of $e")
            note(abs(sl - Wl[2]) <= tol * max(1.0, abs(sl)), "contact: speed differs from flow speed")
        elseif wv.kind in (:fast_fan, :slow_fan)
            fam = wv.kind === :fast_fan ? (σ < 0 ? 1 : 7) : (σ < 0 ? 3 : 5)
            cfu, _, csu = speeds(Wup, γ); cfd, _, csd = speeds(Wdn, γ)
            cu, cd = wv.kind === :fast_fan ? (cfu, cfd) : (csu, csd)
            shead, stail = σ < 0 ? (sl, sr) : (sr, sl)
            e = max(abs(shead - (Wup[2] + σ * cu)), abs(stail - (Wdn[2] + σ * cd))) / max(abs(shead), abs(stail), 1.0)
            note(weak || sl <= sr + tol * max(1.0, abs(sr)), "wave $i ($(wv.kind)): fan speeds not increasing")
            if abs(Wdn[1] - Wup[1]) > 1e-12 * Wup[1] && wv.fan !== nothing
                # (a) pointwise: along the fan, (A(q) - λ I) dq/dτ = 0 with λ the family speed.
                #     Needs no eigen-decomposition, so it stays well conditioned when the
                #     slow and Alfvén (or fast and Alfvén) speeds nearly coincide.
                e = max(e, fan_pointwise_defect(wv.fan, ctx, F, γ))
                # (b) where the family is well separated: independent integration along
                #     eigenvectors of the numerical Jacobian, compared at the fan's end.
                if speed_gap(Wup, γ, wv.kind) > 1e-3 && speed_gap(Wdn, γ, wv.kind) > 1e-3
                    e = max(e, relerr(fan_eigen_integrate(Wup, Wdn[1], fam, γ), Wdn))
                end
            end
            maxerr = max(maxerr, e)
            note(e <= 100 * tol, "wave $i ($(wv.kind)): fan mismatch $e")
        end
        if i > 1
            prev = waves[i-1]
            ok = speed_to_user(prev.s_right, F) <= sl + tol * max(1.0, abs(sl))
            note(ok, "waves $(i-1) and $i overlap")
        end
    end
    return (ok = isempty(msgs), maxerr = maxerr, messages = msgs)
end

"""
    check(sol::RiemannSolution)

Re-run the independent checks on a solution; returns `(ok, maxerr, messages)`.
"""
check(sol::RiemannSolution) = check_waves(sol.waves, sol.ctx, unmirrored(sol.frame), sol.prob.γ)
