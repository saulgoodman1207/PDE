# u0-anchor case summary

| Case | Family | Difficulty | Success | Final t | Residual inf | GridL2 | PSD ratio | Suitability |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| paper_example_5_1 | radial exponential | easy | true | 0 | 1.575e-12 | 2.140e-13 | 1.000 | well suited |
| paper_example_5_1_noexact | radial exponential (no exact diagnostics) | easy | true | 0 | 1.575e-12 | n/a | 1.000 | well suited |
| quartic_convex | quartic polynomial | medium | true | 0 | 9.725e-12 | 1.858e-04 | 0.927 | usable but weaker |
| scaled_quartic_convex | scaled quartic polynomial | medium | true | 0 | 1.945e-12 | 3.716e-05 | 0.927 | usable but weaker |
| mild_exp_separable | separable exponential | easy | true | 0 | 5.232e-12 | 1.981e-11 | 1.000 | well suited |
| mixed_polynomial_convex | mixed polynomial with xy coupling | easy-medium | true | 0 | 1.544e-12 | 2.650e-13 | 1.000 | well suited |
| mild_sextic_convex | separable sextic polynomial | easy-medium | true | 0 | 6.118e-12 | 9.714e-13 | 1.000 | well suited |
| gaussian_bump_convex | quadratic plus Gaussian bump | easy-medium | true | 0 | 1.323e-11 | 1.224e-11 | 1.000 | well suited |
| trig_convex | trigonometric near-degenerate convex | hard | false | 1 | 0.000e+00 | 1.283e-02 | 0.757 | not suitable yet |

## Case notes

- `paper_example_5_1`: smooth strongly convex exact solution Solver mode: default adaptive constrained run. Continuation reaches t = 0 cleanly with strong convexity diagnostics.
- `paper_example_5_1_noexact`: same solution, exact diagnostics disabled Solver mode: default adaptive constrained run. Continuation reaches t = 0 cleanly with strong convexity diagnostics.
- `quartic_convex`: smooth polynomial convex solution Solver mode: default adaptive constrained run. Solver reaches t = 0, but convexity/accuracy diagnostics are noticeably weaker.
- `scaled_quartic_convex`: same quartic shape with weaker curvature scale Solver mode: default adaptive constrained run. Solver reaches t = 0, but convexity/accuracy diagnostics are noticeably weaker.
- `mild_exp_separable`: mild-curvature separable convex solution Solver mode: default adaptive constrained run. Continuation reaches t = 0 cleanly with strong convexity diagnostics.
- `mixed_polynomial_convex`: smooth convex polynomial with mixed term Solver mode: default adaptive constrained run. Continuation reaches t = 0 cleanly with strong convexity diagnostics.
- `mild_sextic_convex`: smooth sextic convex solution Solver mode: default adaptive constrained run. Continuation reaches t = 0 cleanly with strong convexity diagnostics.
- `gaussian_bump_convex`: smooth strongly convex quadratic with local bump Solver mode: default adaptive constrained run. Continuation reaches t = 0 cleanly with strong convexity diagnostics.
- `trig_convex`: oscillatory convex solution near degeneracy Solver mode: fixed-step unconstrained stress test. hard case did not reach t = 0 under the current settings.
