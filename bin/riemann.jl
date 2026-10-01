#!/usr/bin/env julia
# Command-line front end:  julia --project=<pkg> bin/riemann.jl problem.toml [out.csv]
# File format: see `problem_load` (src/problemfile.jl) and examples/.
using ExactMHDRiemannSolver

function main(args)
    isempty(args) && (println(stderr, "usage: riemann.jl problem.toml [out.csv]"); return 2)
    L, R, γ, t, x = try
        problem_load(args[1])
    catch err
        err isa ArgumentError || rethrow()
        println(stderr, err.msg); return 2
    end
    sol = solve(RiemannProblem(L, R; γ))
    show(stdout, MIME"text/plain"(), sol); println()
    if sol.retcode != Success
        println(stderr, "no solution: ", sol.retcode, " (", sol.reason, ")")
        return 1
    end
    path = length(args) >= 2 ? args[2] : splitext(args[1])[1] * ".csv"
    write_csv(path, sol; t, x)
    println("wrote ", path)
    return 0
end

exit(main(ARGS))
