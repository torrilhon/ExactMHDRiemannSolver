# Problem files (TOML), shared by the REPL and the command-line front end.
#
#   gamma = 1.6666666666666667            # optional, default 5/3
#   left  = [rho, vx, vy, vz, Bx, By, Bz, p]
#   right = [rho, vx, vy, vz, Bx, By, Bz, p]
#   [output]                              # optional
#   t = 0.4                               # default 1.0
#   x = [-1.0, 1.0]                       # default [-1, 1]
#   n = 2001                              # default 2001

"""
    problem_load(path) -> (; L, R, γ, t, x)

Read a problem file (TOML, see `examples/`). Returns the primitive states `L`, `R`
(ρ, vx, vy, vz, Bx, By, Bz, p), the adiabatic exponent `γ` and the output time `t`
and grid `x` from the optional `[output]` table. A relative `path` that does not exist
in the current directory is also looked up in the package directory, so
`problem_load("examples/paper.toml")` works for an installed package too.

```julia
L, R, γ, t, x = problem_load("examples/paper.toml")
sol = solve(RiemannProblem(L, R; γ))
write_csv("paper.csv", sol; t, x)
```
"""
function problem_load(path::AbstractString)
    file = path
    if !isfile(file) && !isabspath(path)
        alt = joinpath(pkgdir(@__MODULE__), path)
        isfile(alt) && (file = alt)
    end
    isfile(file) || throw(ArgumentError("problem file not found: $path"))
    cfg = TOML.parsefile(file)
    state(key) = begin
        haskey(cfg, key) || throw(ArgumentError("$path: missing `$key`"))
        v = cfg[key]
        (v isa AbstractVector && length(v) == 8 && all(x -> x isa Real, v)) ||
            throw(ArgumentError("$path: `$key` must be 8 numbers (ρ, vx, vy, vz, Bx, By, Bz, p)"))
        Float64.(v)
    end
    L, R = state("left"), state("right")
    γ = Float64(get(cfg, "gamma", 5 / 3))
    out = get(cfg, "output", Dict{String,Any}())
    t = Float64(get(out, "t", 1.0))
    xr = Float64.(get(out, "x", [-1.0, 1.0]))
    length(xr) == 2 || throw(ArgumentError("$path: `output.x` must be [xmin, xmax]"))
    x = range(xr[1], xr[2]; length = Int(get(out, "n", 2001)))
    return (; L, R, γ, t, x)
end
