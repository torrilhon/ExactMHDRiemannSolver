# Stress scan: harder random problems than the test suite. Run: julia --project=. benchmark/stress.jl 1000
using ExactMHDRiemannSolver, Random
function side(rng, Bn; bmax=3.0, ratio=100.0)
    bt = exp(rand(rng)*log(bmax/0.01))*0.01; th = 2π*rand(rng)
    ρ = exp((2rand(rng)-1)*log(ratio)/2); p = exp((2rand(rng)-1)*log(ratio)/2)
    [ρ, 3*(2rand(rng)-1), 2rand(rng)-1, 2rand(rng)-1, Bn, bt*cos(th), bt*sin(th), p]
end
function gen(rng, i)
    Bn = exp((2rand(rng)-1)*log(30))*0.3
    L = side(rng, Bn); R = side(rng, Bn)
    if i % 5 == 0
        R[6:7] = -hypot(R[6],R[7]) .* L[6:7] ./ hypot(L[6],L[7])
    end
    return L, R, rand(rng, (5/3, 1.4, 2.0))
end
rng = MersenneTwister(7)
counts = Dict{String,Int}(); worst = 0.0; tmax = 0.0; bad = Int[]
for i in 1:(isempty(ARGS) ? 1000 : parse(Int, ARGS[1]))
    L, R, γ = gen(rng, i)            # every 5th problem is coplanar
    t = @elapsed sol = solve(RiemannProblem(L, R; γ))
    i > 1 && (global tmax = max(tmax, t))
    k = string(sol.retcode) * (sol.reason === :none ? "" : ":" * string(sol.reason))
    counts[k] = get(counts, k, 0) + 1
    sol.retcode == Success && (global worst = max(worst, sol.check.maxerr))
    sol.retcode == CheckFailed && push!(bad, i)
end
println(sort(collect(counts), by = x -> -x[2]))
println("worst check err (Success) = $worst, max time = $(round(tmax, digits=2)) s, CheckFailed cases: $bad")
