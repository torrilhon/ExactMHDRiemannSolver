# Legacy C code

The original implementation of the exact MHD Riemann solver, written by M. Torrilhon
(2002–2008, small changes in 2020) alongside the report

> M. Torrilhon, *Exact Solver and Uniqueness Conditions for Riemann Problems of Ideal
> Magnetohydrodynamics*, Research Report 2002-06, SAM, ETH Zürich (2002).

It is kept here for reference and reproducibility only. It is not used by the Julia
package and is superseded by it. The files in `src/` are unchanged from the original
(comments are in German, some files are ISO-8859-1 encoded).

## Files

| File | Content |
| --- | --- |
| `MHD_Main1.c` | main program: hard-coded problem, command-line options, continuation, database I/O |
| `MHD_Residuum.c` | residual of the 5×5 system; `Half()` builds the waves of one side |
| `MHD_fShock.c` | fast shock for given Mach number (cubic Hugoniot, bisection + secant) |
| `MHD_sShock.c` | slow shock for given drop of Bt, including the compound wave |
| `MHD_sShockMax.c` | slow shock of maximal Mach number (limit of the regular slow shock) |
| `MHD_Rarefaction.c` | fast and slow rarefaction fans (ODE integration with dopri5) |
| `MHD_fRareMax.c` | maximal strength of a fast fan (integration to the Bt → 0 singularity) |
| `MHD_OutPut.c` | writes the solution to `Erg.dat` and a wave table to stdout |
| `Newton.c` | Newton's method with finite-difference Jacobian and matrix inversion |
| `DiffGlgn.c`, `dopri5.h` | Dormand–Prince 5(4) integrator of Hairer & Wanner |
| `tools.c`, `tools.h` | command-line parsing, timers |
| `MHD1.h` | shared declarations (`cField` state struct) |
| `input.txt`, `DataBase.dat`, `DataBase.dat.old` | example problems with solved path variables |

The unknowns are the path variables (ψf⁻, ψs⁻, α_R, ψs⁺, ψf⁺) of the report, Sec. 6:
ψ > 0 are shocks, ψ < 0 rarefactions, α_R is the direction of Bt in the middle states.

## Build and run

```
cd src
make            # uses g++ (the sources are compiled as C++)
./rpx
```

- The problem (left and right states and a starting guess for the path variables) is
  hard-coded in `main()` of `MHD_Main1.c`; γ is the global `gam` (5/3) in the same file.
  Changing either means recompiling.
- `-f file` reads a problem (two states and five path variables, format as in
  `input.txt`). Without `-p 1` the program walks from the hard-coded problem to the
  one in the file by continuation in `-n` steps (default 10); with `-p 1` it uses the
  path variables from the file as the starting guess.
- `-db k` starts from the k-th entry of `DataBase.dat`; `-w 1` appends a converged
  solution to `DataBase.dat`.
- Output: the Newton residuals and the wave table on stdout, and `Erg.dat` with
  columns x, ρ, vx, vy, vz, Bx, By, Bz, p. There x = (wave speed)·0.08, that is the
  solution at the fixed time t = 0.08.
- States are given as (ρ, vx, vy, vz, Bx, By, Bz, p); Bx must be positive.

## Known issues (code review, 2026)

The code reproduces Tables 1–3 of the report (with twist angle 1.5; the report's "1/2"
is a typo). A review found these problems:

1. **Spurious intermediate-shock solutions.** When a slow shock turns Bt negative,
   `Half()` in `MHD_Residuum.c` changes `W[3]` but not `W[4]`, and `MHD_OutPut.c` treats
   the rotation differently again. Newton can then converge to a "solution" that is
   discontinuous at the contact. About 6 % of converged random non-planar problems
   were affected, with no warning apart from the message "Intermediate Wave...".
2. **Bx < 0 gives wrong results.** The Alfvén jump assumes Bx > 0, but negative Bx is
   not rejected or transformed; Newton converges to states that violate the jump conditions.
3. **Fast-shock root finder fails at round-off.** `NST()` in `MHD_fShock.c` (and
   `MHD_sShockMax.c`) iterates the secant method to 1e-15; when two function values
   coincide it divides by zero, giving ρ = 0 and p = ∞ at random shock strengths. This
   is the main cause of non-convergence.
4. **NaN counts as converged.** `NormInf()` ignores NaN, so Newton stops and prints a NaN
   result without a warning.
5. **No damping in Newton.** Divergent iterations can reach path variables of order
   10³, after which dopri5 runs for very long.
6. **Degenerate inputs.** Bt = 0 on a side, Bx = 0 and identical states give NaN; the
   switch to the Bt = 0 formulas at Bt/√p < 1e-3 in `fShock()` is discontinuous.
7. **Continuation.** `n/nstep` in `Continuation()` is an integer division, so the
   transverse velocity is not interpolated.
8. **Smaller points.** Hard-coded γ, problem and output time; `%12.8f` output columns can
   merge; `delete` on memory from `new[]`; unchecked `fopen`; command-line options
   matched with `strstr`.

The report itself has a typo in eq. (54), second line (an extra factor (c/c_A)²); the
code uses the correct relation.

## Patch

`patches/review-fixes.patch` fixes items 3, 4 and 7 with minimal changes (guarded secant
step, NaN-aware convergence test, floating-point continuation). With it, the share of
random non-planar problems that converge from the default guess rose from 40 % to about
72 %. Apply it with

```
cd src && patch -p1 < ../patches/review-fixes.patch && make
```

Items 1, 2, 5 and 6 are design issues; they are addressed in the Julia package, not
in this code.
