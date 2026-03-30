---
name: pde-solver-coding
description: Write, refactor, and verify code for solving partial differential equations with finite-difference, spectral, pseudospectral, and method-of-lines approaches. Use when the user asks for PDE solver implementations, discretization choices, boundary-condition handling, stability or consistency checks, convergence experiments, manufactured-solution tests, or debugging numerical PDE code for equations such as heat, wave, advection-diffusion, Poisson, Burgers, or reaction-diffusion.
---

# PDE Solver Coding

Write numerical PDE code that is explicit about the model, discretization, solver structure, and verification.

Default to the repository's existing language and numerical stack. If there is no established stack, prefer Python with NumPy, SciPy, and Matplotlib for research prototypes.

## Default workflow

1. Lock the mathematical problem.
Extract the PDE, domain, initial condition, boundary condition, coefficients, source term, and target outputs.
If the request omits a critical item, make the smallest reasonable assumption and state it.

2. Pick the numerical family deliberately.
Use `references/method-selection.md` to decide between finite differences and spectral methods.
State the spatial discretization, the time integrator if applicable, and the expected spatial and temporal orders.

3. Separate the code into stable pieces.
Prefer this structure:
- parameter block
- grid construction
- operator construction
- initial and boundary data
- time stepping or linear solve
- diagnostics and error measurement
- plotting or saved outputs

4. Make indexing and boundary treatment explicit.
Do not hide boundary conditions inside ad hoc slices.
For finite differences, state the stencil and how the boundary enters the discrete system.
For spectral methods, state the basis, collocation points, transform conventions, and any dealiasing rule.

5. Verify before polishing.
Prefer an exact solution, manufactured solution, or conserved quantity check.
Run at least one mesh-refinement or mode-refinement study when the environment allows.
Use `references/implementation-checklist.md` as the final pass.

6. Explain numerical limits.
Call out CFL restrictions, stiffness, conditioning, aliasing, Gibbs effects, or compatibility requirements when they matter.

## Output standard

When writing or revising solver code, aim for the following:

- Keep the discretization readable enough that the scheme can be reconstructed from the code.
- Name arrays by mathematical role, not by vague placeholders.
- Put physical parameters and numerical parameters in one visible block.
- Avoid mixing plotting code into the solver core.
- Expose grid spacing, timestep, and polynomial or Fourier resolution directly.
- Report the norm used for error measurement.
- Prefer short driver functions over monolithic scripts when the code is more than a few dozen lines.

## Finite-difference rules

- State the grid layout: node-centered, cell-centered, or staggered.
- State the stencil order and whether the method is explicit, implicit, or semi-implicit.
- For explicit time stepping, check the timestep against the relevant stability restriction.
- For implicit schemes, assemble the linear system so the role of boundary conditions remains visible.
- For nonlinear terms, avoid silent reuse of stale states; show whether the update is fully explicit, linearized, or iterative.

## Spectral and pseudospectral rules

- Match the basis to the geometry.
Use Fourier for periodic problems.
Use Chebyshev-style grids for nonperiodic smooth problems on simple intervals when spectral accuracy is appropriate.
- Build differentiation operators or FFT-based derivatives in one place.
- For nonlinear Fourier pseudospectral terms, apply a dealiasing rule when the nonlinearity can generate unresolved modes.
- Warn when low regularity or discontinuities make spectral convergence unrealistic.

## Debugging workflow

When the user asks to debug PDE code, inspect in this order:

1. PDE sign conventions and coefficient placement
2. Grid spacing and timestep definitions
3. Boundary-condition enforcement
4. Operator scaling factors
5. Initial-condition projection or interpolation
6. Linear solver conditioning or singularity
7. Error metric and reference solution

Prefer to isolate one failure mode at a time with a minimal reproducible case.

## Reference map

- `references/method-selection.md`: choose between spectral and finite-difference approaches
- `references/implementation-checklist.md`: final coding and verification checklist
