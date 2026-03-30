# Method Selection

Use this file when deciding how to discretize a PDE or when the user asks whether a spectral or finite-difference method is more appropriate.

## Quick choice

- Use finite differences when the geometry is simple, the coefficients or solution have limited regularity, the boundary treatment is important, or the user wants a robust baseline implementation.
- Use Fourier spectral or pseudospectral methods when the domain is periodic and the solution is smooth enough that high-order accuracy is worth the tighter modeling assumptions.
- Use Chebyshev-style spectral methods on simple nonperiodic intervals when the solution is smooth and boundary accuracy matters more than code simplicity.

## Finite differences

Prefer finite differences when:

- the domain is an interval or a tensor-product box and the user wants transparent code
- mixed or nonperiodic boundary conditions are central
- the solution may develop steep gradients, layers, or nonsmooth features
- the goal is a dependable reference solver or a baseline for comparison

Tradeoffs:

- easier boundary-condition handling
- easier debugging of local stencil errors
- lower per-degree-of-freedom accuracy than spectral methods on smooth solutions
- stability restrictions can be severe for explicit schemes

## Fourier spectral or pseudospectral

Prefer Fourier spectral methods when:

- the problem is periodic
- the solution is smooth over the full domain
- the user wants high accuracy with relatively few modes
- FFT-based differentiation keeps the implementation compact

Watch for:

- aliasing in nonlinear terms
- incorrect wave-number ordering
- mismatch between nonperiodic boundaries and Fourier assumptions
- poor behavior when the solution is not smooth

## Chebyshev-style spectral

Prefer Chebyshev-style methods when:

- the domain is a simple bounded interval
- the solution is smooth but not periodic
- boundary resolution and high accuracy are both important

Watch for:

- ill-conditioning at high resolution
- dense differentiation matrices in simple prototype implementations
- endpoint clustering that changes timestep restrictions for explicit schemes

## Time discretization guide

- Use explicit schemes for nonstiff problems and cheap prototypes, but check the stability bound.
- Use Crank-Nicolson, backward Euler, BDF, or IMEX-style updates for diffusive or stiff terms when explicit steps are too restrictive.
- State clearly which terms are treated explicitly and which are treated implicitly.

## Minimum justification to include in an answer

Whenever choosing a method, state:

1. why the geometry and boundary conditions fit the method
2. expected spatial order or spectral convergence assumption
3. timestep restriction or linear-solve consequence
4. one verification strategy
