# u0-anchor: necessary conditions vs. additional empirical requirements

## Conclusion

The image conditions

- `u` is a sufficiently smooth classical solution,
- `f` is strictly positive in `Omega`,
- `g` is sufficiently smooth on `partial Omega`,

look like **necessary baseline conditions** for the current `u0_anchor` code, but they are **not sufficient by themselves**.

The clearest counterexample is `trig_convex`: it is smooth, uses smooth Dirichlet data, and has positive `f`, but it still fails under the current anchored continuation because the convex branch is close to degeneracy and the solution is more oscillatory.

## Formal Table

| Case | Satisfies baseline conditions? | Structural character | Outcome with current method | What this tells us |
| --- | --- | --- | --- | --- |
| `paper_example_5_1` | yes | smooth, strongly convex, radial exponential, `f >= 1` | success, very strong diagnostics | baseline conditions are enough when the convex branch is well separated from degeneracy |
| `paper_example_5_1_noexact` | yes | same geometry as above, but no exact diagnostics | success | same conclusion as above; exact solution availability is not required for convergence |
| `mild_exp_separable` | yes | smooth separable convex solution, mild curvature, `f > 0` | success, very strong diagnostics | smooth positive-data cases with gentle curvature are a very good fit |
| `mixed_polynomial_convex` | yes | smooth convex polynomial with mixed term, `f >= 1` | success, very strong diagnostics | moderate coupling is fine if the branch remains uniformly convex |
| `mild_sextic_convex` | yes | smooth sextic convex polynomial, `f >= 1` | success, very strong diagnostics | higher-order polynomial growth is still fine when convexity stays strong |
| `gaussian_bump_convex` | yes | quadratic background plus smooth local bump, `f > 0` | success, very strong diagnostics | localized smooth perturbations are also well handled |
| `quartic_convex` | yes | smooth quartic convex polynomial, stronger curvature variation | success, but weaker accuracy/PSD diagnostics | baseline conditions are not enough to guarantee equally good quality; curvature variation already stresses the method |
| `scaled_quartic_convex` | yes | same quartic structure with weaker overall scale | success, but weaker accuracy/PSD diagnostics | smaller curvature scale can make the branch quality weaker even when the case still converges |
| `trig_convex` | yes | smooth oscillatory convex solution, `f > 0` but close to zero in places, near-degenerate | failure | baseline conditions alone do not guarantee success; near-degeneracy and oscillation can break the current homotopy path |

## What seems necessary

From the experiments above, the current code appears to need at least:

- smooth classical `u`, so the Legendre-Galerkin discretization and Hessian diagnostics make sense,
- smooth Dirichlet data `g`,
- strictly positive `f` in the interior, so the Monge-Ampere target stays elliptic in the classical convex setting.

These are best viewed as **necessary entry conditions** for the class of examples we tested.

## What still seems to be missing

The successful and failed cases suggest that the method also benefits from several extra properties that are **not written in the image conditions**:

1. **Uniform strict convexity away from degeneracy**

   The most successful cases have Hessians that stay comfortably positive definite along the branch. The failed `trig_convex` case is much closer to degeneracy, even though `f` remains positive.

2. **A branch that is not too oscillatory**

   Mild exponential or polynomial cases behave well. The more oscillatory trigonometric case is much harder for the current anchor and Newton continuation to follow.

3. **A Poisson-style anchor that lands near the correct convex branch**

   The method starts from `Delta u0 = 2 sqrt(f)`. This works well when that anchor is already close to the desired convex family, but it is much less reliable when the target solution has a different geometric character.

4. **Moderate Hessian variation**

   Even among successful cases, quartic-type examples show noticeably weaker quality than the exponential and mild-polynomial cases. So strong variation in curvature can degrade branch quality before full failure occurs.

## Working Rule Of Thumb

At the moment, a safer empirical statement is:

> The current `u0_anchor` method is well suited for smooth manufactured Monge-Ampere problems on the square with smooth boundary data, strictly positive `f`, and a solution that stays uniformly and non-degenerately convex with relatively mild geometric variation.

And the complementary warning is:

> Smoothness plus `f > 0` is not enough on its own; near-degenerate or strongly oscillatory convex branches can still fail.
