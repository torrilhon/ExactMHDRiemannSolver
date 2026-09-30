# ExactMHDRiemannSolver.jl

Exact solutions of Riemann problems of one-dimensional ideal magnetohydrodynamics
(γ-law gas), for benchmarking numerical schemes. The algorithm follows

> M. Torrilhon, *Exact Solver and Uniqueness Conditions for Riemann Problems of Ideal
> Magnetohydrodynamics*, Research Report 2002-06, Seminar für Angewandte Mathematik,
> ETH Zürich (2002).

Version 0.1 computes **regular solutions**: fast and slow Lax shocks, fast and slow
rarefaction fans, rotational (Alfvén) discontinuities and the contact. Every solution
it returns has been checked independently; inputs it cannot handle safely are refused
with a stated reason.

## How to cite

If you use this code, please cite the report above. A `CITATION.cff` file is included,
so GitHub shows a "Cite this repository" button; releases will carry a Zenodo DOI.

## Installation

The package is not registered. Install it from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/torrilhon/ExactMHDRiemannSolver")
```

Julia 1.10 or newer. The first call compiles the nonlinear and ODE solvers (about
20 s); afterwards a solve takes a few milliseconds.

## Quick start

```julia
using ExactMHDRiemannSolver

#     (ρ,   vx,  vy,  vz,  Bx,  By,       Bz,       p)
L = [3.0, 0.0, 0.0, 0.0, 1.5, 1.0,      0.0,      3.0]
R = [1.0, 0.0, 0.0, 0.0, 1.5, cos(1.5), sin(1.5), 1.0]

sol = solve(RiemannProblem(L, R; γ = 5/3))
sol.retcode            # Success
display(sol)           # wave table: kind, speeds, state right of each wave
sample(sol, 0.25)      # state at x/t = 0.25
write_csv("paper.csv", sol; t = 0.4, x = range(-1, 1; length = 2001))
```

From the command line:

```
julia --project=. bin/riemann.jl examples/paper.toml out.csv
```

## Interface

| Function | Purpose |
| --- | --- |
| `RiemannProblem(L, R; γ = 5/3)` | problem; states are primitive 8-vectors (ρ, vx, vy, vz, Bx, By, Bz, p), Bx equal on both sides |
| `solve(prob; opts = SolverOptions(), guess = nothing)` | solve; never throws for physical reasons |
| `sol.retcode`, `sol.reason` | outcome, see below |
| `sample(sol, ξ)`, `sample(sol, x, t)` | exact state at ξ = x/t, fans included |
| `wavetable(sol)` | one row per wave: kind, speeds, state right of it |
| `write_csv(path, sol; t, x)` | sampled solution as CSV |
| `check(sol)` | re-run the independent checks |

Return codes:

| `retcode` | Meaning |
| --- | --- |
| `Success` | converged, all independent checks passed |
| `InvalidInput` | non-finite values, ρ ≤ 0, p ≤ 0, γ ≤ 1, or Bx differs between the states |
| `Unsupported` | outside the v0.1 domain: `:perpendicular` (Bn ≈ 0), `:switch_on_off` (Bt/√p below the thresholds below), `:extreme_ratio` |
| `RegularLimit` | the solution needs a limit of the regular waves: `:vacuum`, `:fast_switch_off`, `:slow_limit` (a compound or intermediate wave would be needed), `:bt_small` |
| `NoConvergence` | no solution found from any start, including the homotopy |
| `CheckFailed` | converged but failed the independent checks; please report it |

All thresholds are fields of `SolverOptions`, in units of √p of the state concerned:

| Option | Default | Meaning |
| --- | --- | --- |
| `bn_min` | 1e-3 | minimum Bn/√p of both input states |
| `bt_min_larger` | 1e-4 | minimum Bt/√p of the input state with the larger value |
| `bt_min_smaller` | 1e-6 | minimum Bt/√p of the input state with the smaller value |
| `bt_floor` | 1e-8 | minimum Bt/√p of every middle state; below it the result is `RegularLimit` |
| `time_limit` | 5 s | wall-clock budget of the homotopy fallback |

So one side may carry an almost vanishing transverse field as long as the other one
does not; with both sides small the solver gives up more often (as `NoConvergence`,
never with a wrong answer).

## How it works

1. **Canonical frame.** Every input is mapped by exact symmetries of ideal MHD to a
   frame with the larger Bt/√p on the left, Bn > 0, the left Bt along +y, zero
   transverse velocity on the right and units ρ_L = p_L = 1. Results are mapped back.
   The side swap combines x → −x with B → −B, so Bn keeps its sign; it keeps the
   homotopy start (right state = left state) away from a degenerate state.
2. **Unknowns.** Five path variables Ψ = (ψf⁻, ψs⁻, α_R, ψs⁺, ψf⁺) as in the report,
   with maps chosen so that every Ψ ∈ ℝ⁵ gives a valid regular wave sequence:
   - fast shock: ψ = √D − √D_f with D = γ(M² − ĉ_A²) and D_f its value at M = ĉ_f
     (ψ > 0); near switch-on (tiny Bt, c_A > a) the downstream field then grows linearly
     in ψ, not like √ψ;
   - fast fan: ψ = log(Bt/Bt₀) (ψ ≤ 0), which makes the fan ODE smooth up to switch-off;
   - slow shock: drop of B̂t equal to Δmax·tanh(ψ/Δmax), capped below the regular limit Δmax;
   - slow fan: arc length ς of the curve (s, Bt/√p₀), s = log(ρ₀/ρ), with
     ς = s_vac·tanh(−ψ/s_vac); this is s for ordinary fans and Bt near switch-on, where
     Bt grows like √s.
3. **Kernels.** The fast-shock Hugoniot is written as a cubic in the downstream
   transverse field (after dividing out the trivial root); its largest real root is the
   fast shock, found by bracketing. Unlike the cubic in v̂ of the report, it stays
   accurate to machine precision near switch-on. The slow-shock quadratic is rewritten for
   t = 1 − v̂ so it stays well conditioned for weak shocks. Fans are integrated with
   OrdinaryDiffEq (Vern9, tolerance 1e-12); where c_s ≈ c_A the slow-fan equations use
   (c_A² − c_s²) = c_A² B_t²/(ρ(c_f² − c_A²)) to avoid cancellation. Derivatives through root finders use the
   implicit function theorem, so ForwardDiff gives exact Jacobians.
4. **Solve.** NonlinearSolve trust region, on a residual scaled by the problem's own
   pressure, speed and field scales. Starts: α_R = α/2, α/2 + π, and 0 and π for coplanar
   data; then a homotopy from the trivial problem.
5. **Checks.** In the physical frame: Rankine–Hugoniot with the full conservative flux
   for every discontinuity, Lax type (1→2 fast, 3→4 slow) and entropy, fan end speeds,
   the eigen relation (A(q) − λI)·dq/dτ = 0 at points inside each fan (no
   eigen-decomposition, so it stays reliable when wave speeds nearly coincide), an
   independent integration along eigenvectors of the numerical Jacobian where the
   family is well separated, continuity at the contact, and wave ordering.

One wave builder (`build_side`) produces the waves used by the residual, the output
and the checks, so they cannot disagree.

## Validation

- Reproduces Tables 1–2 of the report (twist angle 1.5) to the printed digits.
- Brio–Wu (γ = 2), Ryu–Jones 1a and 2a agree with reference solutions.
- 2000 random non-planar problems: all `Success`; 355 of them agree with independently
  computed reference solutions (`test/data`, 8 decimals) to 1.3e-8.
- A harder stress set of 1000 problems (field ratios up to 1:300, density and pressure
  ratios up to 1:100, strong flows, 20 % coplanar, γ ∈ {1.4, 5/3, 2}): 97.8 % `Success`,
  2.2 % `RegularLimit` (vacuum generation), no `NoConvergence`, no wrong answers.
- Small transverse field, 30 random problems per setting, both regimes c_A ≷ a on the
  small side, larger side Bt/√p from 1e-2 down to 1e-4 and smaller side from 1e-4 down
  to 1e-6: for a > c_A at least 29/30 `Success` in every setting; for c_A > a 27/30
  (3 genuine vacuum limits) while the larger side is ≥ 1e-3, and 24–25/30 when it is
  3e-4 or 1e-4, the rest ending as `NoConvergence`. Worst check error 1.2e-11.
- Invariance under Galilean shifts, transverse rotations, B → −B and x → −x.

Run the tests with `julia --project=. test/runtests.jl` (about 60 s after compilation)
or `Pkg.test()`.

## Notes on the report

- Sec. 2.1 states twist angle α = 1/2, and Table 2 lists cos(0.5), sin(0.5); the
  tables and figures correspond to α = 1.5.
- Eq. (54), second line, carries an extra factor (c/c_A)². The correct relation,
  used here, is dvt/ds = ∓ c (Bt/Bn) / (1 − c²/c_A²) · e_t.
- Table 1 prints vz = 0.438321 behind the slow fan; the correct value is 0.438329.

## Roadmap

- Bn → 0: below a small threshold switch to a quasi-Euler solver (fast waves and a
  tangential discontinuity); the regular solver itself stays accurate down to
  Bn/√p ≈ 1e-8.
- Later: intermediate and compound waves.

## License

MIT, see `LICENSE`.
