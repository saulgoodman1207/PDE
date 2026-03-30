function monge_ampere_homotopy_sine_galerkin_demo()
% Prototype solver for a homotopy-regularized Monge-Ampere-type equation
% on Omega = (0,1)^2:
%
%   (1 - t) * ( -eps * Delta^2 u ) + t * det(D^2 u) = f,   0 <= t <= 1.
%
% The code uses:
% 1. tensor-product sine Galerkin basis for homogeneous Dirichlet boundary,
% 2. continuation in t from 0 to 1,
% 3. damped Newton iterations on the spectral coefficients.
%
% This is a practical prototype. It follows the paper's "continuation +
% spectral discretization + Newton" structure, but it uses a sine-Galerkin
% basis instead of the paper's Legendre basis phi_n = L_n - L_{n+2}.
% For a fully faithful implementation of the paper, the next step is to
% replace the sine basis and uniform-grid projection with:
%   - Legendre basis,
%   - Legendre-Gauss quadrature,
%   - L2 projection of nonlinear terms.

    clc;
    close all;

    problem = default_problem();
    result = solve_homotopy_problem(problem);

    fprintf('Final continuation parameter t = %.2f\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    fprintf('Relative L2 grid error         = %.3e\n', result.rel_l2_error);
    fprintf('Max grid error                 = %.3e\n', result.max_error);

    plot_solution(result);
end

function problem = default_problem()
% Manufactured smooth solution for testing:
%   u_exact = sin(pi x) sin(pi y)
% Then the target equation at t = 1 is det(D^2 u) = f.

    problem.N = 8;              % number of sine modes in each direction
    problem.M = 64;             % quadrature/evaluation points per direction
    problem.epsilon = 1e-2;     % biharmonic regularization strength
    problem.t_values = [0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.96, 1.00];
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.fd_step = 1e-7;
    problem.line_search_min = 2^-12;

    problem.exact_u = @(X, Y) sin(pi * X) .* sin(pi * Y);
    problem.f = @(X, Y) exact_rhs(X, Y);
end

function F = exact_rhs(X, Y)
    s = sin(pi * X) .* sin(pi * Y);
    c = cos(pi * X) .* cos(pi * Y);
    F = pi^4 * (s.^2 - c.^2);
end

function result = solve_homotopy_problem(problem)
    N = problem.N;
    M = problem.M;
    eps_reg = problem.epsilon;
    t_values = problem.t_values(:).';

    [grid, basis] = build_sine_galerkin_operators(N, M);
    F_grid = problem.f(grid.X, grid.Y);
    F_hat = project_to_sine_basis(F_grid, basis);

    lambda = basis.lambda;
    U_hat = -F_hat ./ (eps_reg * lambda.^2);
    U_vec = U_hat(:);

    residual_history = zeros(numel(t_values), 1);
    newton_steps = zeros(numel(t_values), 1);

    for k = 1:numel(t_values)
        t = t_values(k);
        if k == 1
            residual_history(k) = norm(residual(U_vec, t, eps_reg, F_hat, basis), inf);
            newton_steps(k) = 0;
            continue;
        end

        [U_vec, iter_count, res_norm] = newton_solve( ...
            U_vec, ...
            @(z) residual(z, t, eps_reg, F_hat, basis), ...
            problem.newton_tol, ...
            problem.newton_maxit, ...
            problem.fd_step, ...
            problem.line_search_min);

        residual_history(k) = res_norm;
        newton_steps(k) = iter_count;

        fprintf('t = %5.2f, Newton steps = %2d, residual = %.3e\n', t, iter_count, res_norm);
    end

    U_hat = reshape(U_vec, N, N);
    U_grid = basis.Sx * U_hat * basis.Sy.';
    U_exact = problem.exact_u(grid.X, grid.Y);
    err = U_grid - U_exact;

    result.grid = grid;
    result.U_grid = U_grid;
    result.U_exact = U_exact;
    result.error = err;
    result.rel_l2_error = norm(err(:)) / norm(U_exact(:));
    result.max_error = norm(err(:), inf);
    result.coefficients = U_hat;
    result.t_values = t_values;
    result.residual_history = residual_history;
    result.newton_steps = newton_steps;
end

function [grid, basis] = build_sine_galerkin_operators(N, M)
    x = (1:M)' / (M + 1);
    y = x;
    [X, Y] = meshgrid(x, y);

    modes = 1:N;
    px = pi * modes;
    py = pi * modes;

    Sx = sin(x * px);
    Sy = sin(y * py);
    Cx = cos(x * px);
    Cy = cos(y * py);

    Sxx = Sx .* (-(px .^ 2));
    Syy = Sy .* (-(py .^ 2));
    Cx1 = Cx .* px;
    Cy1 = Cy .* py;

    [mx, my] = ndgrid(modes, modes);
    lambda = (pi^2 * (mx.^2 + my.^2));

    grid.x = x;
    grid.y = y;
    grid.X = X;
    grid.Y = Y;
    grid.dx = 1 / (M + 1);
    grid.dy = 1 / (M + 1);

    basis.N = N;
    basis.M = M;
    basis.Sx = Sx;
    basis.Sy = Sy;
    basis.Sxx = Sxx;
    basis.Syy = Syy;
    basis.Cx1 = Cx1;
    basis.Cy1 = Cy1;
    basis.lambda = lambda;
    basis.dx = grid.dx;
    basis.dy = grid.dy;
end

function r = residual(U_vec, t, eps_reg, F_hat, basis)
    N = basis.N;
    U_hat = reshape(U_vec, N, N);

    U_xx = basis.Sxx * U_hat * basis.Sy.';
    U_yy = basis.Sx  * U_hat * basis.Syy.';
    U_xy = basis.Cx1 * U_hat * basis.Cy1.';

    det_hessian = U_xx .* U_yy - U_xy .^ 2;
    det_hat = project_to_sine_basis(det_hessian, basis);

    linear_hat = -(1 - t) * eps_reg * (basis.lambda .^ 2) .* U_hat;
    r_hat = linear_hat + t * det_hat - F_hat;
    r = r_hat(:);
end

function G_hat = project_to_sine_basis(G, basis)
% Approximate L2 projection onto the sine basis using a uniform interior grid.
% Since int_0^1 sin(m pi x)^2 dx = 1/2, the coefficient is
%   G_hat(m,n) = 4 * int int G(x,y) sin(m pi x) sin(n pi y) dx dy.

    G_hat = 4 * basis.dx * basis.dy * (basis.Sx.' * G * basis.Sy);
end

function [z, iter_count, res_norm] = newton_solve(z0, F, tol, maxit, fd_step, line_search_min)
    z = z0;
    r = F(z);
    res_norm = norm(r, inf);

    for iter_count = 1:maxit
        if res_norm < tol
            return;
        end

        J = finite_difference_jacobian(F, z, r, fd_step);
        delta = -(J + 1e-12 * eye(size(J))) \ r;

        alpha = 1.0;
        accepted = false;
        while alpha >= line_search_min
            z_trial = z + alpha * delta;
            r_trial = F(z_trial);
            if norm(r_trial, inf) < (1 - 1e-4 * alpha) * res_norm
                z = z_trial;
                r = r_trial;
                res_norm = norm(r, inf);
                accepted = true;
                break;
            end
            alpha = alpha / 2;
        end

        if ~accepted
            z = z + delta;
            r = F(z);
            res_norm = norm(r, inf);
        end
    end
end

function J = finite_difference_jacobian(F, z, Fz, h)
    n = numel(z);
    J = zeros(n, n);

    for j = 1:n
        dz = zeros(n, 1);
        step = h * max(1, abs(z(j)));
        dz(j) = step;
        J(:, j) = (F(z + dz) - Fz) / step;
    end
end

function plot_solution(result)
    figure('Color', 'w', 'Name', 'Homotopy Spectral Solver');

    subplot(2, 2, 1);
    surf(result.grid.X, result.grid.Y, result.U_grid, 'EdgeColor', 'none');
    title('Numerical Solution');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    view(35, 30);
    colorbar;

    subplot(2, 2, 2);
    surf(result.grid.X, result.grid.Y, result.U_exact, 'EdgeColor', 'none');
    title('Exact Solution');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    view(35, 30);
    colorbar;

    subplot(2, 2, 3);
    surf(result.grid.X, result.grid.Y, abs(result.error), 'EdgeColor', 'none');
    title('|Error|');
    xlabel('x');
    ylabel('y');
    zlabel('|u - u_{exact}|');
    view(35, 30);
    colorbar;

    subplot(2, 2, 4);
    semilogy(result.t_values, max(result.residual_history, eps), '-o', 'LineWidth', 1.5);
    title('Continuation Residual');
    xlabel('t');
    ylabel('Residual inf-norm');
    grid on;
end
