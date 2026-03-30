function monge_ampere_homotopy_legendre_galerkin_demo()
% Legendre-Galerkin prototype for the homotopy problem
%
%   (1 - t) * ( -eps * Delta^2 u ) + t * det(D^2 u) = f
%
% on Omega = (-1,1)^2.
%
% This version uses the paper's basis and, unlike the previous direct
% assembly prototype, it now builds the discrete operators explicitly from
% the paper's B-matrix chain:
%
%   phi_n(x) = L_n(x) - L_{n+2}(x),
%   B1: forward Legendre transform,
%   B3: map from phi-coefficients to Legendre coefficients,
%   B4: backward Legendre transform,
%   B5: second-derivative coefficient map,
%   B6: first-derivative coefficient map.
%
% The nonlinear term is assembled through
%
%   T1[(T2 U).*(T3 U)] - T1[(T4 U).^2],
%
% mirroring the paper's matrix notation.

    clc;
    close all;

    problem = default_problem();
    result = solve_homotopy_problem(problem);

    fprintf('Final continuation parameter t = %.2f\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    fprintf('Relative L2 evaluation error   = %.3e\n', result.rel_l2_error);
    fprintf('Max evaluation error           = %.3e\n', result.max_error);

    plot_solution(result);
end

function problem = default_problem()
% Manufactured smooth solution:
%   u_exact(x,y) = -(1 - x^2)(1 - y^2),
% which lies exactly in span{phi_0(x) phi_0(y)}.

    problem.N = 8;          % number of phi-modes per direction
    problem.Q = 24;         % Legendre-Gauss points per direction
    problem.epsilon = 1e-2;
    problem.t_values = [0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.96, 1.00];
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.fd_step = 1e-7;
    problem.line_search_min = 2^-12;
    problem.eval_points = 81;

    problem.exact_u = @(X, Y) -(1 - X.^2) .* (1 - Y.^2);
    problem.f = @(X, Y) 4 * (1 - X.^2 - Y.^2 - 3 * X.^2 .* Y.^2);
end

function result = solve_homotopy_problem(problem)
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points);

    f_vals = problem.f(ops.Xq, ops.Yq);
    F_vec = ops.T1 * f_vals(:);

    U_vec = solve_linear_biharmonic(problem.epsilon, F_vec, ops);

    residual_history = zeros(numel(problem.t_values), 1);
    newton_steps = zeros(numel(problem.t_values), 1);

    for k = 1:numel(problem.t_values)
        t = problem.t_values(k);
        if k == 1
            residual_history(k) = norm(residual(U_vec, t, problem.epsilon, F_vec, ops), inf);
            newton_steps(k) = 0;
            fprintf('t = %5.2f, Newton steps = %2d, residual = %.3e\n', t, 0, residual_history(k));
            continue;
        end

        [U_vec, iter_count, res_norm] = newton_solve( ...
            U_vec, ...
            @(z) residual(z, t, problem.epsilon, F_vec, ops), ...
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

function ops = build_b_matrix_operators(N, Q, eval_points)
    full_count = N + 2;     % L_0,...,L_{N+1}
    full_deg = full_count - 1;

    [xq, wq] = legendre_gauss(Q);
    Pq = legendre_values(xq, full_deg);
    gamma = 2 ./ (2 * (0:full_deg)' + 1);

    % B1: forward discrete Legendre transform
    B1 = diag(1 ./ gamma) * (Pq.' * diag(wq));

    % B2: diagonal gamma matrix from Legendre orthogonality
    B2 = diag(gamma);

    % B3: phi-coefficients -> full Legendre coefficients
    B3 = zeros(full_count, N);
    for j = 0:N-1
        B3(j + 1, j + 1) = 1;
        B3(j + 3, j + 1) = -1;
    end

    % B4: backward Legendre transform
    B4 = Pq;

    % B5: second-derivative coefficient map in Legendre space
    B5 = zeros(full_count, full_count);
    for j = 0:full_deg
        for k = 0:full_deg
            if mod(k + j, 2) == 0 && k <= j - 2
                B5(k + 1, j + 1) = (k + 0.5) * (j * (j + 1) - k * (k + 1));
            end
        end
    end

    % B6: first-derivative coefficient map in Legendre space
    B6 = zeros(full_count, full_count);
    for j = 0:full_deg
        for k = 0:full_deg
            if mod(k + j, 2) == 1 && k <= j - 1
                B6(k + 1, j + 1) = 2 * k + 1;
            end
        end
    end

    % 1D phi evaluation and nodal-to-phi projection.
    Phi_q = B4 * B3;
    Pphi = B3.' * B2 * B1;

    % Linear biharmonic operator T = A⊗B + 2C⊗C + B⊗A
    [A, B, C] = build_abc_matrices(N);
    T = kron(A, B) + 2 * kron(C, C) + kron(B, A);

    % Nonlinear transform chain from phi-coefficients to nodal derivatives.
    coeff_to_legendre = kron(B3, B3);
    nodal_from_legendre = kron(B4, B4);

    T1 = kron(Pphi, Pphi);
    T2 = nodal_from_legendre * kron(eye(full_count), B5) * coeff_to_legendre;
    T3 = nodal_from_legendre * kron(B5, eye(full_count)) * coeff_to_legendre;
    T4 = nodal_from_legendre * kron(B6, B6) * coeff_to_legendre;

    % Evaluation grid for plots.
    xe = linspace(-1, 1, eval_points).';
    ye = xe;
    Pe = legendre_values(xe, full_deg);
    Phi_eval = Pe * B3;
    [Xe, Ye] = meshgrid(xe, ye);
    [Xq, Yq] = meshgrid(xq, xq);

    ops.N = N;
    ops.Q = Q;
    ops.full_count = full_count;
    ops.B1 = B1;
    ops.B2 = B2;
    ops.B3 = B3;
    ops.B4 = B4;
    ops.B5 = B5;
    ops.B6 = B6;
    ops.A = A;
    ops.B = B;
    ops.C = C;
    ops.T = T;
    ops.T1 = T1;
    ops.T2 = T2;
    ops.T3 = T3;
    ops.T4 = T4;
    ops.Phi_q = Phi_q;
    ops.Phi_eval = Phi_eval;
    ops.xq = xq;
    ops.wq = wq;
    ops.Xq = Xq;
    ops.Yq = Yq;
    ops.Xe = Xe;
    ops.Ye = Ye;
end

function [A, B, C] = build_abc_matrices(N)
    A = zeros(N, N);
    B = zeros(N, N);
    C = zeros(N, N);

    for j = 0:N-1
        C(j + 1, j + 1) = 4 * j + 6;

        B(j + 1, j + 1) = 2 / (2 * j + 1) + 2 / (2 * j + 5);
        if j + 2 <= N - 1
            B(j + 3, j + 1) = -2 / (2 * j + 5);
            B(j + 1, j + 3) = B(j + 3, j + 1);
        end

        for k = 0:N-1
            if mod(k + j, 2) == 0
                A(k + 1, j + 1) = 2 * (2 * k + 3) * (2 * j + 3) * (k + 1) * (k + 2);
            end
        end
    end
end

function U_vec = solve_linear_biharmonic(eps_reg, F_vec, ops)
    U_vec = -(1 / eps_reg) * (ops.T \ F_vec);
end

function r = residual(U_vec, t, eps_reg, F_vec, ops)
    u_xx_vals = ops.T2 * U_vec;
    u_yy_vals = ops.T3 * U_vec;
    u_xy_vals = ops.T4 * U_vec;

    nonlin_vec = ops.T1 * (u_xx_vals .* u_yy_vals - u_xy_vals .^ 2);
    r = -(1 - t) * eps_reg * (ops.T * U_vec) + t * nonlin_vec - F_vec;
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
% Golub-Welsch Gauss-Legendre nodes and weights on [-1,1].

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
    figure('Color', 'w', 'Name', 'Legendre-Galerkin B-Matrix Solver');

    subplot(2, 2, 1);
    surf(result.X_eval, result.Y_eval, result.U_eval, 'EdgeColor', 'none');
    title('Numerical Solution');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    view(35, 30);
    colorbar;

    subplot(2, 2, 2);
    surf(result.X_eval, result.Y_eval, result.U_exact, 'EdgeColor', 'none');
    title('Exact Solution');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    view(35, 30);
    colorbar;

    subplot(2, 2, 3);
    surf(result.X_eval, result.Y_eval, abs(result.error), 'EdgeColor', 'none');
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
