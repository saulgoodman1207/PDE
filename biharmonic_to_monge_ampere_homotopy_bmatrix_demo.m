function biharmonic_to_monge_ampere_homotopy_bmatrix_demo()
% Fixed-right-hand-side continuation for
%
%   t * Delta^2 u + (1 - t) * det(D^2 u) = f,
%
% on Omega = (-1,1)^2, using the paper's Legendre-Galerkin basis and the
% explicit B-matrix transform chain.
%
% Setup:
% 1. Choose u_exact(x,y) = sin(pi x) sin(pi y),
% 2. Fix the target right-hand side from the nonlinear endpoint:
%       f(x,y) = det(D^2 u_exact),
% 3. Continue from t = 1 down to t = 0.
%
% Only the final t = 0 problem has u_exact as an exact solution. The
% intermediate t-level problems are continuation problems.

    clc;
    close all;

    problem = default_problem();
    result = solve_fixed_rhs_homotopy(problem);

    fprintf('Final continuation parameter t = %.2f\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    fprintf('Final relative L2 error        = %.3e\n', result.rel_l2_error);
    fprintf('Final max error                = %.3e\n', result.max_error);
    fprintf('Final Jacobian sigma_min       = %.3e\n', result.jacobian_sigma_min(end));
    fprintf('Final Jacobian cond            = %.3e\n', result.jacobian_cond(end));
    fprintf('Final Hessian lambda_min min   = %.3e\n', result.hessian_lambda_min(end));
    fprintf('Final Hessian lambda_max max   = %.3e\n', result.hessian_lambda_max(end));
    fprintf('Final Hessian PSD ratio        = %.3f\n', result.hessian_psd_ratio(end));
    fprintf('Final Hessian indefinite ratio = %.3f\n', result.hessian_indef_ratio(end));

    plot_solution(result);
end

function problem = default_problem()
    problem.N = 10;
    problem.Q = 32;
    problem.eval_points = 101;
    problem.t_start = 1.00;
    problem.t_end = 0.00;
    problem.initial_step = 0.05;
    problem.min_step = 0.0025;
    problem.max_step = 0.10;
    problem.step_growth = 1.20;
    problem.accept_tol = 1e-8;
    problem.max_continuation_steps = 120;
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.line_search_min = 2^-12;

    problem.exact_u = @(X, Y) sin(pi * X) .* sin(pi * Y);
    problem.fixed_f = @(X, Y) det_hessian_exact(X, Y);
    problem.boundary_extension = @(X, Y) sin(pi * X) .* sin(pi * Y);
end

function D = det_hessian_exact(X, Y)
    s = sin(pi * X) .* sin(pi * Y);
    c = cos(pi * X) .* cos(pi * Y);
    D = pi^4 * (s.^2 - c.^2);
end

function result = solve_fixed_rhs_homotopy(problem)
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points);
    lift = build_boundary_lift(problem, ops);

    F_vals = problem.fixed_f(ops.Xq, ops.Yq);
    F_vec = ops.T1 * F_vals(:);

    % Paper-style initialization: start the first homotopy level from zero.
    U_vec = zeros(problem.N^2, 1);

    t_history = zeros(problem.max_continuation_steps, 1);
    residual_history = zeros(problem.max_continuation_steps, 1);
    newton_steps = zeros(problem.max_continuation_steps, 1);
    jacobian_sigma_min = zeros(problem.max_continuation_steps, 1);
    jacobian_cond = zeros(problem.max_continuation_steps, 1);
    hessian_lambda_min = zeros(problem.max_continuation_steps, 1);
    hessian_lambda_max = zeros(problem.max_continuation_steps, 1);
    hessian_psd_ratio = zeros(problem.max_continuation_steps, 1);
    hessian_nsd_ratio = zeros(problem.max_continuation_steps, 1);
    hessian_indef_ratio = zeros(problem.max_continuation_steps, 1);

    [U_vec, iter_count, res_norm] = newton_solve( ...
        U_vec, ...
        @(z) residual(z, problem.t_start, F_vec, ops, lift), ...
        @(z) jacobian(z, problem.t_start, ops, lift), ...
        problem.newton_tol, ...
        problem.newton_maxit, ...
        problem.line_search_min);

    accepted = 1;
    current_t = problem.t_start;
    step = problem.initial_step;
    [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, lift));
    [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, lift);
    t_history(accepted) = current_t;
    residual_history(accepted) = res_norm;
    newton_steps(accepted) = iter_count;
    jacobian_sigma_min(accepted) = sigma_min;
    jacobian_cond(accepted) = cond_val;
    hessian_lambda_min(accepted) = lam_min;
    hessian_lambda_max(accepted) = lam_max;
    hessian_psd_ratio(accepted) = psd_ratio;
    hessian_nsd_ratio(accepted) = nsd_ratio;
    hessian_indef_ratio(accepted) = indef_ratio;
    fprintf(['t = %5.3f, Newton steps = %2d, residual = %.3e, sigma_min = %.3e, cond = %.3e, ' ...
        'lam_min = %.3e, lam_max = %.3e, psd = %.3f, indef = %.3f\n'], ...
        current_t, iter_count, res_norm, sigma_min, cond_val, lam_min, lam_max, psd_ratio, indef_ratio);

    while current_t > problem.t_end && accepted < problem.max_continuation_steps
        target_t = max(problem.t_end, current_t - step);
        U_try0 = U_vec;

        [U_try, iter_count, res_norm] = newton_solve( ...
            U_try0, ...
            @(z) residual(z, target_t, F_vec, ops, lift), ...
            @(z) jacobian(z, target_t, ops, lift), ...
            problem.newton_tol, ...
            problem.newton_maxit, ...
            problem.line_search_min);

        if res_norm <= problem.accept_tol
            accepted = accepted + 1;
            current_t = target_t;
            U_vec = U_try;
            [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, lift));
            [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, lift);
            t_history(accepted) = current_t;
            residual_history(accepted) = res_norm;
            newton_steps(accepted) = iter_count;
            jacobian_sigma_min(accepted) = sigma_min;
            jacobian_cond(accepted) = cond_val;
            hessian_lambda_min(accepted) = lam_min;
            hessian_lambda_max(accepted) = lam_max;
            hessian_psd_ratio(accepted) = psd_ratio;
            hessian_nsd_ratio(accepted) = nsd_ratio;
            hessian_indef_ratio(accepted) = indef_ratio;
            fprintf(['t = %5.3f, Newton steps = %2d, residual = %.3e, step = %.4f, sigma_min = %.3e, cond = %.3e, ' ...
                'lam_min = %.3e, lam_max = %.3e, psd = %.3f, indef = %.3f\n'], ...
                current_t, iter_count, res_norm, step, sigma_min, cond_val, lam_min, lam_max, psd_ratio, indef_ratio);

            if current_t > problem.t_end
                step = min(problem.max_step, step * problem.step_growth);
            end
        else
            new_step = step / 2;
            fprintf('reject t = %5.3f, residual = %.3e, step %.4f -> %.4f\n', target_t, res_norm, step, new_step);
            step = new_step;
            if step < problem.min_step
                warning('ContinuationFailed:MinStep', ...
                    'Adaptive continuation stopped before reaching t_end; required step %.4g below min_step %.4g.', ...
                    step, problem.min_step);
                break;
            end
        end
    end

    U_hat = reshape(U_vec, problem.N, problem.N);
    U_eval = ops.Phi_eval * U_hat * ops.Phi_eval.';
    U_exact = problem.exact_u(ops.Xe, ops.Ye);
    u_total_eval = U_eval + lift.eval_vals;
    err = u_total_eval - U_exact;

    result.coefficients = U_hat;
    result.t_values = t_history(1:accepted);
    result.residual_history = residual_history(1:accepted);
    result.newton_steps = newton_steps(1:accepted);
    result.jacobian_sigma_min = jacobian_sigma_min(1:accepted);
    result.jacobian_cond = jacobian_cond(1:accepted);
    result.hessian_lambda_min = hessian_lambda_min(1:accepted);
    result.hessian_lambda_max = hessian_lambda_max(1:accepted);
    result.hessian_psd_ratio = hessian_psd_ratio(1:accepted);
    result.hessian_nsd_ratio = hessian_nsd_ratio(1:accepted);
    result.hessian_indef_ratio = hessian_indef_ratio(1:accepted);
    result.X_eval = ops.Xe;
    result.Y_eval = ops.Ye;
    result.U_eval = u_total_eval;
    result.U_exact = U_exact;
    result.error = err;
    result.rel_l2_error = norm(err(:)) / max(norm(U_exact(:)), eps);
    result.max_error = norm(err(:), inf);
    result.final_t = current_t;
end

function lift = build_boundary_lift(problem, ops)
    if isfield(problem, "boundary_extension") && isa(problem.boundary_extension, "function_handle")
        Ub_vals = problem.boundary_extension(ops.Xq, ops.Yq);
        eval_vals = problem.boundary_extension(ops.Xe, ops.Ye);
    else
        Ub_vals = zeros(size(ops.Xq));
        eval_vals = zeros(size(ops.Xe));
    end

    Ub_vec = Ub_vals(:);
    Ub_leg_vec = ops.forward_legendre * Ub_vec;
    Id = ops.Id_full;

    Ub_xx_vals = ops.nodal_from_legendre * kron(Id, ops.B5) * Ub_leg_vec;
    Ub_yy_vals = ops.nodal_from_legendre * kron(ops.B5, Id) * Ub_leg_vec;
    Ub_xy_vals = ops.nodal_from_legendre * kron(ops.B6, ops.B6) * Ub_leg_vec;

    Ub_lap_vals = Ub_xx_vals + Ub_yy_vals;
    Ub_lap_leg_vec = ops.forward_legendre * Ub_lap_vals;
    Ub_lap_xx_vals = ops.nodal_from_legendre * kron(Id, ops.B5) * Ub_lap_leg_vec;
    Ub_lap_yy_vals = ops.nodal_from_legendre * kron(ops.B5, Id) * Ub_lap_leg_vec;
    Ub_biharm_vals = Ub_lap_xx_vals + Ub_lap_yy_vals;
    Ub_biharm_phi = ops.T1 * Ub_biharm_vals;

    lift.xx_vals = Ub_xx_vals;
    lift.yy_vals = Ub_yy_vals;
    lift.xy_vals = Ub_xy_vals;
    lift.biharm_phi = Ub_biharm_phi;
    lift.eval_vals = eval_vals;
end

function [sigma_min, cond_val] = jacobian_diagnostics(J)
    s = svd(J);
    sigma_min = s(end);
    if sigma_min == 0
        cond_val = inf;
    else
        cond_val = s(1) / sigma_min;
    end
end

function [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, lift)
    u_xx = ops.T2 * U_vec + lift.xx_vals;
    u_yy = ops.T3 * U_vec + lift.yy_vals;
    u_xy = ops.T4 * U_vec + lift.xy_vals;

    trace_h = u_xx + u_yy;
    disc = sqrt(max((u_xx - u_yy) .^ 2 + 4 * u_xy .^ 2, 0));
    lambda1 = 0.5 * (trace_h - disc);
    lambda2 = 0.5 * (trace_h + disc);

    tol = 1e-10;
    psd_mask = (lambda1 >= -tol) & (lambda2 >= -tol);
    nsd_mask = (lambda1 <= tol) & (lambda2 <= tol);
    indef_mask = ~(psd_mask | nsd_mask);

    lam_min = min(lambda1);
    lam_max = max(lambda2);
    psd_ratio = mean(psd_mask);
    nsd_ratio = mean(nsd_mask);
    indef_ratio = mean(indef_mask);
end

function r = residual(U_vec, t, F_vec, ops, lift)
    u_xx_vals = ops.T2 * U_vec + lift.xx_vals;
    u_yy_vals = ops.T3 * U_vec + lift.yy_vals;
    u_xy_vals = ops.T4 * U_vec + lift.xy_vals;

    biharm_vec = ops.T * U_vec;
    det_vec = ops.T1 * (u_xx_vals .* u_yy_vals - u_xy_vals .^ 2);

    r = t * biharm_vec + (1 - t) * det_vec - F_vec + t * lift.biharm_phi;
end

function J = jacobian(U_vec, t, ops, lift)
    u_xx_vals = ops.T2 * U_vec + lift.xx_vals;
    u_yy_vals = ops.T3 * U_vec + lift.yy_vals;
    u_xy_vals = ops.T4 * U_vec + lift.xy_vals;

    nonlin_jac = ops.T1 * ( ...
        diag(u_yy_vals) * ops.T2 + ...
        diag(u_xx_vals) * ops.T3 - ...
        2 * diag(u_xy_vals) * ops.T4);

    J = t * ops.T + (1 - t) * nonlin_jac;
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
    forward_legendre = kron(B1, B1);
    Id_full = speye(full_count);
    Pphi = B3.' * B2 * B1;

    T1 = kron(Pphi, Pphi);
    T2 = nodal_from_legendre * kron(Id_full, B5) * coeff_to_legendre;
    T3 = nodal_from_legendre * kron(B5, Id_full) * coeff_to_legendre;
    T4 = nodal_from_legendre * kron(B6, B6) * coeff_to_legendre;

    [A, B, C] = build_abc_matrices(N);
    T = kron(A, B) + 2 * kron(C, C) + kron(B, A);

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
    ops.T = T;
    ops.T1 = T1;
    ops.T2 = T2;
    ops.T3 = T3;
    ops.T4 = T4;
    ops.Phi_eval = Phi_eval;
    ops.Xq = Xq;
    ops.Yq = Yq;
    ops.Xe = Xe;
    ops.Ye = Ye;
    ops.coeff_to_legendre = coeff_to_legendre;
    ops.nodal_from_legendre = nodal_from_legendre;
    ops.forward_legendre = forward_legendre;
    ops.full_count = full_count;
    ops.Id_full = Id_full;
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

function [z, iter_count, res_norm] = newton_solve(z0, F, Jfun, tol, maxit, line_search_min)
    z = z0;
    r = F(z);
    res_norm = norm(r, inf);

    for iter_count = 1:maxit
        if res_norm < tol
            return;
        end

        J = Jfun(z);
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
    figure('Color', 'w', 'Name', 'Biharmonic to Monge-Ampere Homotopy');

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
