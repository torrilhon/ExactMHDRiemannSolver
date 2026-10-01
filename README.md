# ExactMHDRiemannSolver.jl

Exact solutions of Riemann problems of one-dimensional ideal magnetohydrodynamics
(γ-law gas), for benchmarking numerical schemes. The algorithm follows

> M. Torrilhon, [*Exact Solver and Uniqueness Conditions for Riemann Problems of Ideal
> Magnetohydrodynamics*](https://www.sam.math.ethz.ch/sam_reports/reports_final/reports2002/2002-06.pdf),
> SAM Research Report 2002-06, ETH Zürich, 2002.

Also see:

> M. Torrilhon, *Uniqueness conditions for Riemann problems of ideal
> magnetohydrodynamics*, Journal of Plasma Physics **69**(3), 253–276 (2003),
> doi:[10.1017/S0022377803002186](https://doi.org/10.1017/S0022377803002186).

Version 0.1.0 computes **regular solutions**: fast and slow Lax shocks, fast and slow
rarefaction fans, rotational (Alfvén) discontinuities and the contact; for a vanishing
normal field it uses a quasi-Euler solver (see below). A solution is reported as
`Success` only after it has passed a set of independent checks; inputs outside the
supported domain are refused with a stated reason. Intermediate and compound waves are
not computed yet.

## Supported initial states

The input is checked before any solve. Field thresholds refer to |B|/√p with p the
pressure of the same state.

| Quantity | Condition | Otherwise |
| --- | --- | --- |
| all values, γ | finite, γ > 1 | `InvalidInput` |
| density, pressure | ρ > 0 and p > 0 on both sides | `InvalidInput` |
| normal field Bn | equal on both sides to within 1e-12·max(\|Bn\|, 1) | `InvalidInput` (`:bn_jump`) |
| jump ratios | 1e-6 ≤ ρ_R/ρ_L ≤ 1e6 and 1e-6 ≤ p_R/p_L ≤ 1e6 | `Unsupported` (`:extreme_ratio`) |
| velocities | no restriction | — |

| Case | Normal field | Transverse field | Solver |
| --- | --- | --- | --- |
| vanishing normal field | \|Bn\|/√p ≤ 1e-10 on both sides | any, including Bt = 0 on one or both sides | quasi-Euler (Bn = 0) |
| regular | \|Bn\|/√p > 1e-10 on at least one side | \|Bt\|/√p ≥ 1e-4 on the larger side and ≥ 1e-6 on the smaller side | regular |
| otherwise | \|Bn\|/√p > 1e-10 on at least one side | either condition violated, e.g. Bt = 0 on one side | `Unsupported` (`:switch_on_off`) |

The thresholds are the defaults of `bn_euler`, `bt_min_larger`, `bt_min_smaller` and
`ratio_max`. A problem that passes can still end as `RegularLimit` when its solution
needs a limit of the regular waves: vacuum, a middle state with |Bt|/√p < 1e-8
(`bt_floor`), or an intermediate or compound wave. With both transverse fields close to
their thresholds a few percent of problems end as `NoConvergence`.

## How to cite

If you use this code, please cite the report and the article above. A
[`CITATION.cff`](CITATION.cff) file is included,
so GitHub shows a "Cite this repository" button. A Zenodo DOI for releases is planned.

## Installation

The package is not registered. Install it from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/torrilhon/ExactMHDRiemannSolver")
```

Julia 1.10 or newer. The first call compiles the nonlinear and ODE solvers, which can
take tens of seconds. After that, typical problems solve in milliseconds; problems that
need the homotopy fallback take longer, bounded by `time_limit`.

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

The same problem as a file (`examples/` also has Brio–Wu); `problem_load` returns the
states, γ and the output time and grid of the file:

```julia
L, R, γ, t, x = problem_load("examples/paper.toml")
sol = solve(RiemannProblem(L, R; γ))
write_csv("paper.csv", sol; t, x)
```

A relative path is looked up in the current directory first, then in the package
directory, so the bundled examples load from anywhere.

From the command line:

```
julia --project=. bin/riemann.jl examples/paper.toml out.csv
```

## Interface

| Function | Purpose |
| --- | --- |
| `RiemannProblem(L, R; γ = 5/3)` | problem; states are primitive 8-vectors (ρ, vx, vy, vz, Bx, By, Bz, p), Bx equal on both sides |
| `solve(prob; opts = SolverOptions(), guess = nothing)` | solve; reports failures through `retcode` instead of throwing |
| `sol.retcode`, `sol.reason` | outcome, see below |
| `sample(sol, ξ)`, `sample(sol, x, t)` | state at ξ = x/t; fan interiors from the dense ODE output |
| `wavetable(sol)` | one row per wave: kind, speeds, state right of it |
| `write_csv(path, sol; t, x)` | sampled solution as CSV |
| `check(sol)` | re-run the independent checks |
| `problem_load(path)` | read a problem file (TOML, format in `examples/`): `(; L, R, γ, t, x)` |

Return codes:

| `retcode` | Meaning |
| --- | --- |
| `Success` | converged, all independent checks passed |
| `InvalidInput` | non-finite values, ρ ≤ 0, p ≤ 0, γ ≤ 1, or Bx differs between the states |
| `Unsupported` | outside the v0.1 domain: `:switch_on_off` (Bt/√p below the thresholds below), `:extreme_ratio` |
| `RegularLimit` | the solution needs a limit of the regular waves: `:vacuum`, `:fast_switch_off`, `:slow_limit` (a compound or intermediate wave would be needed), `:bt_small` |
| `NoConvergence` | no solution found from any start, including the homotopy; for Bn ≈ 0 also a shock too strong for double precision (P*/P ≳ 1e14) |
| `CheckFailed` | converged but failed the independent checks; please report it |

All thresholds are fields of `SolverOptions`, in units of √p of the state concerned:

| Option | Default | Meaning |
| --- | --- | --- |
| `bn_euler` | 1e-10 | Bn/√p (both states) at or below which the quasi-Euler solver is used |
| `bt_min_larger` | 1e-4 | minimum Bt/√p of the input state with the larger value |
| `bt_min_smaller` | 1e-6 | minimum Bt/√p of the input state with the smaller value |
| `bt_floor` | 1e-8 | minimum Bt/√p of every middle state; below it the result is `RegularLimit` |
| `ratio_max` | 1e6 | maximum of ρ_R/ρ_L, p_R/p_L and their inverses |
| `time_limit` | 5 s | wall-clock budget of the homotopy fallback |

So one side may carry an almost vanishing transverse field as long as the other one
does not; with both sides small the solver gives up more often. In the tests below
such cases ended as `NoConvergence`, not as solutions that failed the checks. The Bt thresholds apply only to the regular solver: with
Bn/√p ≤ `bn_euler` any Bt is accepted, including Bt = 0 on one or both sides.

**Vanishing normal field.** For Bn/√p ≤ 1e-10 the problem is solved as Bn = 0
(`sol.method === :quasi_euler`). The slow and Alfvén waves then merge with the contact
into a tangential discontinuity (u and p + Bt²/2 continuous; ρ, p, Bt and vt may jump),
and the two fast waves are gas-dynamic waves with Bt/ρ, the direction of Bt and vt
constant and c_f² = (γp + Bt²)/ρ. This system is genuinely nonlinear for every Bt, so
no regular-wave limit other than vacuum (`RegularLimit :vacuum`) can occur. The single
unknown is the total pressure at the contact, found by bracketing. In the test cases,
the regular solver just above Bn/√p = 1e-10 and the quasi-Euler solver just below
agree to within 1e-8 in the outer wave speeds, u* and the total pressure P*.

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
   transverse field (after dividing out the trivial root); its largest real root in the
   physical range is the fast shock, found by bracketing. Unlike the cubic in v̂ of the
   report, it stays well conditioned near switch-on. The slow-shock quadratic is rewritten for
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

The figures below are from the author's test runs for v0.1. They describe the tested
cases, not guarantees for arbitrary input; please report any problem where a `Success`
result looks wrong.

- Reproduces Tables 1–2 of the report (twist angle 1.5) to the printed digits.
- Brio–Wu (γ = 2), Ryu–Jones 1a and 2a agree with reference solutions.
- 2000 random non-planar problems: all `Success`; 355 of them agree with independently
  computed reference solutions (`test/data`, 8 decimals) to 1.3e-8.
- A harder stress set of 1000 problems (field ratios up to 1:300, density and pressure
  ratios up to 1:100, strong flows, 20 % coplanar, γ ∈ {1.4, 5/3, 2}): 97.8 % `Success`,
  2.2 % `RegularLimit` (vacuum generation), no `NoConvergence` and no `CheckFailed`.
- Small transverse field, 30 random problems per setting, both regimes c_A ≷ a on the
  small side, larger side Bt/√p from 1e-2 down to 1e-4 and smaller side from 1e-4 down
  to 1e-6: for a > c_A at least 29/30 `Success` in every setting; for c_A > a 27/30
  (3 genuine vacuum limits) while the larger side is ≥ 1e-3, and 24–25/30 when it is
  3e-4 or 1e-4, the rest ending as `NoConvergence`. Worst check error 1.2e-11.
- Bn = 0: Sod's problem matches Toro (p* = 0.30313, u* = 0.92745); 3800 random
  perpendicular problems with Bt/√p from 0 to 1e3 and ratios up to 1:10⁴: all
  `Success` except cases that generate vacuum.
- Invariance under Galilean shifts, transverse rotations, B → −B and x → −x.

Run the tests with `julia --project=. test/runtests.jl` or `Pkg.test()`; the test
suite covers a subset of the cases above.

## Notes on the SAM report 2002

- Sec. 2.1 states twist angle α = 1/2, and Table 2 lists cos(0.5), sin(0.5); the
  tables and figures correspond to α = 1.5.
- Eq. (54), second line, carries an extra factor (c/c_A)². The correct relation,
  used here, is dvt/ds = ∓ c (Bt/Bn) / (1 − c²/c_A²) · e_t.
- Table 1 prints vz = 0.438321 behind the slow fan; the correct value is 0.438329.

## Roadmap

- Intermediate and compound waves (Bt → 0 on both sides, switch-on/off, 180° problems).

## License

MIT, see `LICENSE`.
