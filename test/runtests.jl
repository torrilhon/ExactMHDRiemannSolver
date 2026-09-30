using ExactMHDRiemannSolver
using Test, Random, DelimitedFiles, LinearAlgebra
using StaticArrays, ForwardDiff
const E = ExactMHDRiemannSolver

const DATA = joinpath(@__DIR__, "data")
const s4 = sqrt(4π)

# Rankine-Hugoniot defect of a jump W0 → W1 at speed s, relative to the flux size.
function rh_defect(W0, W1, s, γ)
    d = s * (E.conserved(W1, γ) - E.conserved(W0, γ)) - (E.flux(W1, γ) - E.flux(W0, γ))
    return maximum(abs, d) / max(maximum(abs, E.flux(W0, γ)), maximum(abs, E.flux(W1, γ)), 1.0)
end

randstate(rng; Bn = 0.3 + 2rand(rng)) = E.HState(0.2 + 3rand(rng), 2rand(rng) - 1, 0.2 + 3rand(rng),
    0.1 + 2rand(rng), 2π * rand(rng), SVector(rand(rng) - 0.5, rand(rng) - 0.5)), Bn

@testset "ExactMHDRiemannSolver" begin

@testset "kernels: Rankine-Hugoniot" begin
    rng = MersenneTwister(1)
    for γ in (5 / 3, 1.4, 2.0), _ in 1:40
        U, Bn = randstate(rng)
        ctx = E.Ctx(γ, Bn, SolverOptions())
        σ = rand(rng, (-1, 1))
        W0 = E.to_prim(U, ctx)
        # fast shock over a range of strengths, including very weak ones
        for ψ in (1e-10, 1e-4, 0.05, 0.5, 3.0)
            D, s = E.fast_shock(U, ψ, σ, ctx)
            @test rh_defect(W0, E.to_prim(D, ctx), s, γ) < 1e-12
            ψ > 1e-6 && @test D.ρ > U.ρ
        end
        # slow shock up to its regular limit
        Δmax = E.slow_limit(U, ctx)
        A = U.bt / sqrt(U.p)
        @test 0 < Δmax <= A
        for f in (1e-7, 1e-3, 0.3, 0.9, 0.999)
            D, s, _, _ = E.slow_shock_state(U, f * Δmax, σ, ctx)
            @test rh_defect(W0, E.to_prim(D, ctx), s, γ) < 1e-12
            @test D.bt > 0 && D.ρ > U.ρ
        end
        # rotation
        D, s = E.rotation(U, U.φ + 1.3, σ, ctx)
        @test rh_defect(W0, E.to_prim(D, ctx), s, γ) < 1e-12
    end
end

@testset "kernels: fans vs eigenvector integration" begin
    rng = MersenneTwister(2)
    for γ in (5 / 3, 2.0), _ in 1:15
        U, Bn = randstate(rng)
        ctx = E.Ctx(γ, Bn, SolverOptions())
        W0 = E.to_prim(U, ctx)
        for σ in (-1, 1)
            D, _ = E.fast_fan(U, -0.7, σ, ctx)
            Wt = E.fan_eigen_integrate(W0, D.ρ, σ < 0 ? 1 : 7, γ)
            @test E.relerr(Wt, E.to_prim(D, ctx)) < 1e-8
            D, _ = E.slow_fan(U, -0.4, σ, ctx)
            Wt = E.fan_eigen_integrate(W0, D.ρ, σ < 0 ? 3 : 5, γ)
            @test E.relerr(Wt, E.to_prim(D, ctx)) < 1e-8
        end
    end
end

@testset "Jacobian: ForwardDiff vs central differences" begin
    L = [3.0, 0, 0, 0, 1.5, 1.0, 0, 3.0]; R = [1.0, 0, 0, 0, 1.5, cos(1.5), sin(1.5), 1.0]
    F = E.make_frame(SVector{8}(L), SVector{8}(R))
    UL, UR = E.to_hstate.(E.canonical_states(L, R, F))
    ctx = E.Ctx(5 / 3, 1.5 / sqrt(3), SolverOptions())
    for Ψ in (SVector(-0.3, -0.05, 0.9, 0.3, 0.08), SVector(0.2, 0.1, 0.4, -0.2, -0.1))
        f(x) = E.residual(x, UL, UR, ctx)
        J = ForwardDiff.jacobian(f, Ψ)
        h = 1e-6
        Jfd = hcat([(f(Ψ + h * SVector(ntuple(k -> k == j ? 1.0 : 0.0, 5))) -
                     f(Ψ - h * SVector(ntuple(k -> k == j ? 1.0 : 0.0, 5)))) / 2h for j in 1:5]...)
        @test maximum(abs, J - Jfd) < 1e-6 * max(1, maximum(abs, J))
    end
end

@testset "paper example, Tables 1-2 (α = 1.5)" begin
    L = [3.0, 0, 0, 0, 1.5, 1.0, 0, 3.0]
    R = [1.0, 0, 0, 0, 1.5, cos(1.5), sin(1.5), 1.0]
    sol = solve(RiemannProblem(L, R; γ = 5 / 3))
    @test sol.retcode == Success
    @test [r.kind for r in wavetable(sol)] ==
          [:fast_fan, :rotation, :slow_fan, :contact, :slow_shock, :rotation, :fast_shock]
    T = wavetable(sol)
    # (s_left, ρ, vx, vy, vz, By, Bz, p) right of each wave; paper values, 6 decimals
    ref = [(-1.474922, 2.340949, 0.348797, -0.144157, 0.0, 0.642777, 0.0, 1.984139),
        (-0.631585, 2.340949, 0.348797, -0.339270, 0.354780, 0.344252, 0.542820, 1.984139),
        (-0.521395, 2.200167, 0.402052, -0.286284, 0.438321, 0.413199, 0.651535, 1.789281),
        (0.402052, 1.408739, 0.402052, -0.286284, 0.438321, 0.413199, 0.651535, 1.789281),
        (1.279598, 1.054703, 0.107484, -0.514217, 0.078923, 0.601050, 0.947741, 1.093004),
        (1.568067, 1.054703, 0.107484, -0.006260, -0.088275, 0.079386, 1.119452, 1.093004),
        (2.072332, 1.0, 0.0, 0.0, 0.0, cos(1.5), sin(1.5), 1.0)]
    for (r, t) in zip(T, ref)
        got = (r.s_left, r.ρ, r.vx, r.vy, r.vz, r.By, r.Bz, r.p)
        # vz of the slow fan tail is printed as 0.438321 in the paper, 0.438329 here and in the C code
        @test maximum(abs.(got .- t)) < 1e-5
    end
    @test T[1].s_right ≈ -0.990247 atol = 2e-6
    @test T[3].s_right ≈ -0.445268 atol = 2e-6
end

@testset "benchmarks against the C code" begin
    cases = [
        (2.0, [1, 0, 0, 0, 0.75, 1, 0, 1], [0.125, 0, 0, 0, 0.75, -1, 0, 0.1],
            [-1.792285, -0.280742, -0.117758, 1.419872, 3.315279], 3.683667),
        (5 / 3, [1.08, 1.2, 0.01, 0.5, 2 / s4, 3.6 / s4, 2 / s4, 0.95], [1, 0, 0, 0, 2 / s4, 4 / s4, 2 / s4, 1],
            [-0.957840, 0.143728, 0.259691, 0.902124, 1.027452, 2.263788], nothing),
        (5 / 3, [1, 10, 0, 0, 5 / s4, 5 / s4, 0, 20], [1, -10, 0, 0, 5 / s4, 5 / s4, 0, 1],
            [-4.802843, -0.116150, 0.723762, 1.406909, 4.600452], nothing)]
    for (γ, L, R, speeds, extra) in cases
        sol = solve(RiemannProblem(L, R; γ))
        @test sol.retcode == Success
        s = [w.s_left for w in wavetable(sol)]
        for v in speeds
            @test minimum(abs.(s .- v)) < 2e-6
        end
        extra === nothing || @test wavetable(sol)[end].s_right ≈ extra atol = 2e-6
    end
end

@testset "random problems: cross-check with the C code" begin
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    REF = readdlm(joinpath(DATA, "c_reference.csv"), ',')
    ref = Dict(Int(REF[i, 1]) => REF[i, :] for i in 1:size(REF, 1))
    nfail, maxd = 0, 0.0
    for i in 1:300
        seed = Int(P[i, 1]); L = P[i, 2:9]; R = P[i, 10:17]
        sol = solve(RiemannProblem(L, R; γ = 5 / 3))
        sol.retcode == Success || (nfail += 1; continue)
        if haskey(ref, seed)
            sc = only(r.s_left for r in wavetable(sol) if r.kind === :contact)
            d = max(maximum(abs.(sample(sol, sc - 1e-9) .- ref[seed][18:25])),
                maximum(abs.(sample(sol, sc + 1e-9) .- ref[seed][26:33])))
            maxd = max(maxd, d)
        end
    end
    @test nfail == 0
    @test maxd < 1e-7          # the C output carries 8 decimals
end

@testset "symmetries" begin
    rng = MersenneTwister(3)
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    rot(W, t) = (c = cos(t); s = sin(t);
        [W[1], W[2], c * W[3] - s * W[4], s * W[3] + c * W[4], W[5], c * W[6] - s * W[7], s * W[6] + c * W[7], W[8]])
    negB(W) = [W[1:4]; -W[5:7]; W[8]]
    shift(W, u, v, w) = [W[1], W[2] + u, W[3] + v, W[4] + w, W[5:8]...]
    mirror(W) = [W[1], -W[2], W[3], W[4], -W[5], W[6], W[7], W[8]]
    ξs = range(-3, 3; length = 61)
    for i in 1:25
        L = P[i, 2:9]; R = P[i, 10:17]
        base = solve(RiemannProblem(L, R))
        base.retcode == Success || continue
        S(sol, f) = [sample(sol, ξ) for ξ in ξs]
        B = S(base, identity)
        t = 2π * rand(rng)
        s1 = solve(RiemannProblem(rot(L, t), rot(R, t)))
        @test maximum(maximum(abs.(sample(s1, ξ) - rot(b, t))) for (ξ, b) in zip(ξs, B)) < 1e-8
        s2 = solve(RiemannProblem(negB(L), negB(R)))
        @test maximum(maximum(abs.(sample(s2, ξ) - negB(b))) for (ξ, b) in zip(ξs, B)) < 1e-8
        u, v, w = 0.7, -0.3, 0.2
        s3 = solve(RiemannProblem(shift(L, u, v, w), shift(R, u, v, w)))
        @test maximum(maximum(abs.(sample(s3, ξ + u) - shift(b, u, v, w))) for (ξ, b) in zip(ξs, B)) < 1e-8
        s4m = solve(RiemannProblem(mirror(R), mirror(L)))
        @test maximum(maximum(abs.(sample(s4m, -ξ) - mirror(b))) for (ξ, b) in zip(ξs, B) if abs(ξ) > 0) < 1e-8
    end
end

@testset "single-wave problems (vanishing waves)" begin
    γ = 5 / 3
    R = [1.0, 0, 0, 0, 0.8, 0.6, 0, 1.0]
    ctx = E.Ctx(γ, 0.8, SolverOptions())
    U = E.to_hstate(R)
    for (kind, make) in ((:fast_shock, u -> E.fast_shock(u, 0.3, 1, ctx)[1]),
                         (:slow_shock, u -> E.slow_shock_state(u, 0.2, 1, ctx)[1]),
                         (:fast_fan, u -> E.fast_fan(u, -0.4, 1, ctx)[1]))
        # left state = state behind one right-moving wave on R: the solution is that wave alone
        L = Vector(E.to_prim(make(U), ctx))
        sol = solve(RiemannProblem(L, R; γ))
        @test sol.retcode == Success
        @test maximum(abs.(sample(sol, 5.0) .- R)) < 1e-10
        @test maximum(abs.(sample(sol, -5.0) .- L)) < 1e-10
    end
end

@testset "small transverse field" begin
    # fast shock near switch-on (cA > a, tiny A): machine-precision Rankine-Hugoniot
    ctx = E.Ctx(5 / 3, 2.0, SolverOptions())
    W0(U) = E.to_prim(U, ctx)
    for A in (1e-3, 1e-5, 1e-7), ψ in (1e-8, 1e-4, 0.1, 1.0)
        U = E.HState(1.0, 0.0, 1.0, A, 0.0, SVector(0.0, 0.0))
        D, s = E.fast_shock(U, ψ, 1, ctx)
        @test rh_defect(W0(U), W0(D), s, 5 / 3) < 1e-13
    end
    # one side with Bt/√p = 1e-5, both regimes; placing it left or right gives mirror images
    P = readdlm(joinpath(DATA, "random_problems.csv"), ',')
    setbt(W, b) = (W = copy(W); f = b * sqrt(W[8]) / hypot(W[6], W[7]); W[6] *= f; W[7] *= f; W)
    for i in 1:6, lowβ in (false, true)
        L, R = copy(P[i, 2:9]), copy(P[i, 10:17])
        if lowβ      # cA = 2a on the small-Bt side
            f = (R[5]^2 / (4 * 5 / 3)) / R[8]; L[8] *= f; R[8] *= f
        end
        sR = solve(RiemannProblem(L, setbt(R, 1e-5)))
        sL = solve(RiemannProblem(setbt(R, 1e-5), L))
        sR.retcode == RegularLimit && continue       # genuine vacuum limit
        @test sR.retcode == Success && sL.retcode == Success
        @test sR.frame.mirror == false && sL.frame.mirror == true
    end
end

@testset "rejection policy" begin
    good = [1.0, 0, 0, 0, 1.0, 1.0, 0, 1.0]
    r(L, R; kw...) = solve(RiemannProblem(L, R; kw...))
    @test r([NaN; good[2:8]], good).retcode == InvalidInput
    @test r([-1.0; good[2:8]], good).retcode == InvalidInput
    @test r(good, [good[1:7]; 0.0]).retcode == InvalidInput
    @test r(good, good; γ = 1.0).retcode == InvalidInput
    @test r(good, [1.0, 0, 0, 0, 2.0, 1.0, 0, 1.0]).reason === :bn_jump
    @test r([1.0, 0, 0, 0, 0.7, 0, 0, 1.0], [0.3, 0, 0, 1, 0.7, 1, 0, 0.2]).reason === :switch_on_off   # RJ4d type
    @test r([1.0, 0, 0, 0, 0, 1, 0, 1.0], [0.125, 0, 0, 0, 0, 0.5, 0, 0.1]).reason === :perpendicular
    @test r(good, [1e-7, 0, 0, 0, 1.0, 1.0, 0, 1.0]).reason === :extreme_ratio
    # transverse-field thresholds: smaller side ≥ 1e-6, larger side ≥ 1e-4 (units of √p)
    bt(b) = [1.0, 0, 0, 0, 1.0, b, 0, 1.0]
    @test r(good, bt(5e-7)).reason === :switch_on_off
    @test r(bt(5e-5), bt(5e-5)).reason === :switch_on_off
    @test r(bt(2e-6), good).retcode != Unsupported
    @test r(bt(2e-4), bt(2e-6)).retcode != Unsupported
    t = r(good, good)
    @test t.retcode == Success && t.reason === :trivial && sample(t, 0.3) ≈ good
    # strongly diverging flow: vacuum is approached, the solver refuses
    v = r([1.0, -8, 0, 0, 1.0, 1.0, 0, 1.0], [1.0, 8, 0, 0, 1.0, 1.0, 0, 1.0])
    @test v.retcode == RegularLimit && v.reason === :vacuum
end

@testset "output" begin
    L = [3.0, 0, 0, 0, 1.5, 1.0, 0, 3.0]; R = [1.0, 0, 0, 0, 1.5, cos(1.5), sin(1.5), 1.0]
    sol = solve(RiemannProblem(L, R))
    @test sample(sol, -10.0) ≈ L
    @test sample(sol, 0.1, 0.4) == sample(sol, 0.25)
    @test_throws ArgumentError sample(sol, 0.1, 0.0)
    @test sample(sol, 10.0) ≈ R
    # inside the left fast fan the speed at the sampled point equals ξ
    W = sample(sol, -1.2)
    @test W[2] - E.speeds(W, 5 / 3)[1] ≈ -1.2 atol = 1e-10
    path = write_csv(tempname() * ".csv", sol; t = 0.4, x = range(-1, 1; length = 101))
    M = readdlm(path, ','; skipstart = 1)
    @test size(M) == (101, 9)
    @test M[1, 2] ≈ 3.0 && M[end, 2] ≈ 1.0
    @test check(sol).ok
end

end
