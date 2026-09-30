# Sampling of the self-similar solution and tabular output, all in the user frame.

"""
    sample(sol, ξ) -> SVector{8}

Primitive state (ρ, vx, vy, vz, Bx, By, Bz, p) at similarity coordinate ξ = x/t.
Fan interiors are evaluated exactly from the dense ODE output.
"""
function sample(sol::RiemannSolution, ξ::Real)
    isfinite(ξ) || throw(ArgumentError("sample needs a finite ξ = x/t"))
    isempty(sol.waves) && throw(ArgumentError("solution has no waves (retcode $(sol.retcode))"))
    F, ctx = sol.frame, sol.ctx
    ξc = (F.mirror ? -ξ : ξ) / sqrt(F.p0 / F.ρ0)   # canonical speed
    for w in sol.waves
        ξc < w.s_left && return user_state(w.left, ctx, F)
        if w.fan !== nothing && ξc <= w.s_right
            f = w.fan
            g(τ) = fan_speed(f, τ, ctx) - ξc
            τ = find_zero(g, (0.0, 1.0), Roots.Brent(); xatol = 1e-14)
            return user_state(fan_state(f, τ, ctx), ctx, F)
        end
    end
    return user_state(sol.waves[end].right, ctx, F)
end

"""
    sample(sol, x, t)

State at position x and time t > 0 (x a number: SVector{8}; x a vector: matrix
length(x) × 8).
"""
function sample(sol::RiemannSolution, x::Real, t::Real)
    t > 0 || throw(ArgumentError("sample needs t > 0"))
    return sample(sol, x / t)
end

function sample(sol::RiemannSolution, x::AbstractVector, t::Real)
    t > 0 || throw(ArgumentError("sample needs t > 0"))
    M = zeros(length(x), 8)
    for (i, xi) in enumerate(x)
        M[i, :] = sample(sol, xi / t)
    end
    return M
end

"""
    wavetable(sol) -> Vector{NamedTuple}

One row per wave, left to right: kind, side, speeds (user frame) and the state
right of the wave.
"""
function wavetable(sol::RiemannSolution)
    F, ctx = sol.frame, sol.ctx
    rows = NamedTuple[]
    # in a mirrored frame the canonical waves run right to left in the user frame
    ws = F.mirror ? reverse(sol.waves) : sol.waves
    for w in ws
        W = user_state(F.mirror ? w.left : w.right, ctx, F)
        a, b = speed_to_user(w.s_left, F), speed_to_user(w.s_right, F)
        push!(rows, (kind = w.kind, s_left = min(a, b), s_right = max(a, b),
            ρ = W[1], vx = W[2], vy = W[3], vz = W[4], Bx = W[5], By = W[6], Bz = W[7], p = W[8]))
    end
    return rows
end

function Base.show(io::IO, ::MIME"text/plain", sol::RiemannSolution)
    show(io, sol)
    isempty(sol.waves) && return
    println(io)
    @printf(io, "%-11s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "wave", "s_left", "s_right",
        "rho", "vx", "vy", "vz", "By", "Bz", "p")
    W = sol.prob.L
    @printf(io, "%-11s %10s %10s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", "left state", "", "",
        W[1], W[2], W[3], W[4], W[6], W[7], W[8])
    for r in wavetable(sol)
        @printf(io, "%-11s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", r.kind,
            r.s_left, r.s_right, r.ρ, r.vx, r.vy, r.vz, r.By, r.Bz, r.p)
    end
end

"""
    write_csv(path, sol; t = 1.0, x = range(-1, 1; length = 2001))

Write the sampled solution at time t as CSV with columns x, rho, vx, vy, vz, Bx, By, Bz, p.
"""
function write_csv(path::AbstractString, sol::RiemannSolution; t::Real = 1.0, x = range(-1, 1; length = 2001))
    M = sample(sol, collect(x), t)
    open(path, "w") do io
        println(io, "x,rho,vx,vy,vz,Bx,By,Bz,p")
        for (i, xi) in enumerate(x)
            println(io, join((repr(Float64(xi)), (repr(M[i, j]) for j in 1:8)...), ","))
        end
    end
    return path
end
