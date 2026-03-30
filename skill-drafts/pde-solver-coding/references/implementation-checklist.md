# Implementation Checklist

Use this file as the final pass before presenting PDE solver code or numerical results.

## Problem setup

- Confirm the PDE signs and coefficient locations.
- Confirm the domain and boundary conditions.
- Confirm whether the grid includes boundaries or interior points only.
- Confirm the norm or diagnostic used to measure success.

## Spatial discretization

- Verify the stencil coefficients or spectral derivative formula.
- Verify scaling by `dx`, `dy`, `dz`, or wave numbers.
- Verify boundary rows, ghost-cell formulas, or basis-specific endpoint handling.
- Verify that multidimensional operators act on the intended axis ordering.

## Time discretization

- Verify the timestep definition and loop bounds.
- Verify whether each term is explicit, implicit, semi-implicit, or split.
- For explicit schemes, check the relevant CFL or diffusive stability restriction.
- For implicit schemes, check matrix dimensions, conditioning, and boundary incorporation.

## Spectral-specific checks

- Verify the grid and mode ordering are consistent.
- Verify transform normalization conventions.
- Verify alias control for nonlinear terms when needed.
- Verify the inverse transform returns a real-valued field when the physics requires it.

## Numerical verification

- Test against an exact or manufactured solution when possible.
- Run at least two or three resolutions and compare observed order.
- Check whether the reported error norm matches the theory claim.
- For conservative PDEs, inspect mass, energy, or another invariant if relevant.

## Presentation

- State all numerical parameters used in the experiment.
- Label figures with the method, resolution, and final time.
- Distinguish spatial error from temporal error when both matter.
- Do not claim convergence without showing the refinement evidence.
