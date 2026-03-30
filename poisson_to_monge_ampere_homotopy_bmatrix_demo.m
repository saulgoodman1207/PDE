function poisson_to_monge_ampere_homotopy_bmatrix_demo()
% Fixed-right-hand-side continuation for
%
%   t * Delta u + (1 - t) * det(D^2 u) = f,
%
% on Omega = (-1,1)^2, using the paper's Legendre-Galerkin basis and the
% explicit B-matrix transform chain.
%
% "True solve" setup:
% 1. Choose u_exact(x,y) = sin(pi x) sin(pi y),
% 2. Define the fixed target right-hand side from the nonlinear endpoint:
%       f(x,y) = det(D^2 u_exact),
% 3. Solve the homotopy from t = 1 down to t = 0.
%
% Only the final t = 0 problem has u_exact as an exact solution. For t > 0,
% the continuation problems are auxiliary problems used to reach the target.

    clc;
    close all;

    problem = default_problem();
    result = solve_fixed_rhs_homotopy(problem);

    fprintf('Final continuation parameter t = %.2f\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    fprintf('Final relative L2 error        = %.3e\n', result.rel_l2_error);
    fprintf('Final max error                = %.3e\n', result.max_error);

    plot_solution(result);
end

function problem = default_problem()
    problem.N = 10;
    problem.Q = 32;
    problem.eval_points = 101;
    problem.t_values = [1.00, 0.95, 0.90, 0.82, 0.72, 0.60, 0.48, 0.36, 0.24, 0.12, 0.06, 0.02, 0.00];
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.fd_step = 1e-7;
    problem.line_search_min = 2^-12;

    problem.exact_u = @(X, Y) sin(pi * X) .* sin(pi * Y);
    problem.fixed_f = @(X, Y) det_hessian_exact(X, Y);
end

function D = det_hessian_exact(X, Y)
    s = sin(pi * X) .* sin(pi * Y);
    c = cos(pi * X) .* cos(pi * Y);
    D = pi^4 * (s.^2 - c.^2);
end

function result = solve_fixed_rhs_homotopy(problem)
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points);

    F_vals = problem.fixed_f(ops.Xq, ops.Yq);
    F_vec = ops.T1 * F_vals(:);

    % Use the exact solution projected into the Galerkin space to verify
    % that the discretization and Newton solve behave correctly near the
    % intended branch.
    U_exact_vals = problem.exact_u(ops.Xq, ops.Yq);
    U_mat = project_to_phi(U_exact_vals, ops);
    U_vec = U_mat(:);

    residual_history = zeros(numel(problem.t_values), 1);
    newton_steps = zeros(numel(problem.t_values), 1);

    for k = 1:numel(problem.t_values)
        t = problem.t_values(k);

        [U_vec, iter_count, res_norm] = newton_solve( ...
            U_vec, ...
            @(z) residual(z, t, F_vec, ops), ...
            problem.newton_tol, ...
            problem.newton_maxit, ...
            problem.fd_step, ...
            problem.line_search_min);

        residual_history(k) = res_norm;
        newton_steps(k) = iter_count;
        fprintf('t = %5.2f, Newton steps = %2d, residual = %.3e\n', t, iter_count, res_norm);
    end

    U_hat = reshape(U_vec, problem.N, problem.N);
    U_eval = ops.Phi_eval * U_hat * ops.Phi_eval.';
    U_exact = problem.exact_u(ops.Xe, ops.Ye);
    err = U_eval - U_exact;

    result.coefficients = U_hat;
    result.t_values = problem.t_values(:);
    result.residual_history = residual_history;
    result.newton_steps = newton_steps;
    result.X_eval = ops.Xe;
    result.Y_eval = ops.Ye;
    result.U_eval = U_eval;
    result.U_exact = U_exact;
    result.error = err;
    result.rel_l2_error = norm(err(:)) / max(norm(U_exact(:)), eps);
    result.max_error = norm(err(:), inf);
end

function r = residual(U_vec, t, F_vec, ops)
    u_xx_vals = ops.T2 * U_vec;
    u_yy_vals = ops.T3 * U_vec;
    u_xy_vals = ops.T4 * U_vec;

    lap_vec = ops.L * U_vec;
    det_vec = ops.T1 * (u_xx_vals .* u_yy_vals - u_xy_vals .^ 2);

    r = t * lap_vec + (1 - t) * det_vec - F_vec;
end

function coeffs = project_to_phi(U_vals, ops)
    Pphi = ops.B3.' * ops.B2 * ops.B1;
    coeffs = Pphi * U_vals * Pphi.';
end

function ops = build_b_matrix_operators(N, Q, eval_points)
    full_count = N + 2;
    full_deg = full_count - 1;

    [xq, wq] = legendre_gauss(Q);
    Pq = legendre_values(xq, full_deg);
    gamma = 2 ./ (2 * (0:full_deg)' + 1);

    B1 = diag(1 ./ gamma) * (Pq.' * diag(wq));
    B2 = diag(gamma);

    B3 = zeros(full_count, N);
    for j = 0:N-1
        B3(j + 1, j + 1) = 1;
        B3(j + 3, j + 1) = -1;
    end

    B4 = Pq;

    B5 = zeros(full_count, full_count);
    for j = 0:full_deg
        for k = 0:full_deg
            if mod(k + j, 2) == 0 && k <= j - 2
                B5(k + 1, j + 1) = (k + 0.5) * (j * (j + 1) - k * (k + 1));
            end
        end
    end

    B6 = zeros(full_count, full_count);
    for j = 0:full_deg
        for k = 0:full_deg
            if mod(k + j, 2) == 1 && k <= j - 1
                B6(k + 1, j + 1) = 2 * k + 1;
            end
        end
    end

    coeff_to_legendre = kron(B3, B3);
    nodal_from_legendre = kron(B4, B4);
    Pphi = B3.' * B2 * B1;

    T1 = kron(Pphi, Pphi);
    T2 = nodal_from_legendre * kron(eye(full_count), B5) * coeff_to_legendre;
    T3 = nodal_from_legendre * kron(B5, eye(full_count)) * coeff_to_legendre;
    T4 = nodal_from_legendre * kron(B6, B6) * coeff_to_legendre;

    % Linear Laplacian operator in weak Galerkin form:
    %   (Delta u_N, v_N) = (u_xx + u_yy, v_N)
    L = T1 * (T2 + T3);

    xe = linspace(-1, 1, eval_points).';
    Pe = legendre_values(xe, full_deg);
    Phi_eval = Pe * B3;
    [Xq, Yq] = meshgrid(xq, xq);
    [Xe, Ye] = meshgrid(xe, xe);

    ops.N = N;
    ops.B1 = B1;
    ops.B2 = B2;
    ops.B3 = B3;
    ops.B4 = B4;
    ops.B5 = B5;
    ops.B6 = B6;
    ops.T1 = T1;
    ops.T2 = T2;
    ops.T3 = T3;
    ops.T4 = T4;
    ops.L = L;
    ops.Phi_eval = Phi_eval;
    ops.Xq = Xq;
    ops.Yq = Yq;
    ops.Xe = Xe;
    ops.Ye = Ye;
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

function [x, w] = legendre_gauss(n)
    beta = 0.5 ./ sqrt(1 - (2 * (1:n-1)) .^ (-2));
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    x = diag(D);
    [x, idx] = sort(x);
    V = V(:, idx);
    w = 2 * (V(1, :) .^ 2).';
end

function P = legendre_values(x, max_degree)
    m = numel(x);
    P = zeros(m, max_degree + 1);
    P(:, 1) = 1;
    if max_degree == 0
        return;
    end

    P(:, 2) = x;
    for n = 1:max_degree-1
        P(:, n + 2) = ((2 * n + 1) * x .* P(:, n + 1) - n * P(:, n)) / (n + 1);
    end
end

function plot_solution(result)
    figure('Color', 'w', 'Name', 'Poisson to Monge-Ampere Homotopy');

    subplot(2, 2, 1);
    surf(result.X_eval, result.Y_eval, result.U_eval, 'EdgeColor', 'none');
    title('Final Numerical Solution');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    view(35, 30);
    colorbar;

    subplot(2, 2, 2);
    surf(result.X_eval, result.Y_eval, result.U_exact, 'EdgeColor', 'none');
    title('Exact Solution at t = 0');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    view(35, 30);
    colorbar;

    subplot(2, 2, 3);
    surf(result.X_eval, result.Y_eval, abs(result.error), 'EdgeColor', 'none');
    title('|Final Error|');
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
