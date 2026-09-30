#!/usr/bin/env julia
# Command-line front end:  julia --project=<pkg> bin/riemann.jl problem.toml [out.csv]
#
# problem.toml:
#   gamma = 1.6666666666666667
#   left  = [rho, vx, vy, vz, Bx, By, Bz, p]
#   right = [rho, vx, vy, vz, Bx, By, Bz, p]
#   [output]            # optional
#   t = 0.4
#   x = [-1.0, 1.0]
#   n = 2001
using ExactMHDRiemann, TOML

function main(args)
    isempty(args) && (println(stderr, "usage: riemann.jl problem.toml [out.csv]"); return 2)
    cfg = TOML.parsefile(args[1])
    γ = Float64(get(cfg, "gamma", 5 / 3))
    prob = RiemannProblem(Float64.(cfg["left"]), Float64.(cfg["right"]); γ)
    sol = solve(prob)
    show(stdout, MIME"text/plain"(), sol); println()
    if sol.retcode != Success
        println(stderr, "no solution: ", sol.retcode, " (", sol.reason, ")")
        return 1
    end
    out = get(cfg, "output", Dict{String,Any}())
    t = Float64(get(out, "t", 1.0)); xr = Float64.(get(out, "x", [-1.0, 1.0])); n = Int(get(out, "n", 2001))
    path = length(args) >= 2 ? args[2] : splitext(args[1])[1] * ".csv"
    write_csv(path, sol; t, x = range(xr[1], xr[2]; length = n))
    println("wrote ", path)
    return 0
end

exit(main(ARGS))
