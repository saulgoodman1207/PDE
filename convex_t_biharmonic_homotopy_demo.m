function convex_t_biharmonic_homotopy_demo()
% Fixed-right-hand-side continuation for
%
%   t * Delta^2 u + (1 - t) * det(D^2 u) = f
%
% on the physical square (0,1)^2, mapped to the reference square [-1,1]^2.
% The nonlinear endpoint uses the convex manufactured solution
%
%   u(x,y) = (x^2 + y^2)/2 + x^4 + y^4,
%
% so
%
%   Delta^2 u = 48,    det(D^2 u) = (1 + 12x^2)(1 + 12y^2).
%
% We freeze the right-hand side at the nonlinear endpoint:
%
%   f(x,y) = det(D^2 u_exact) = (1 + 12x^2)(1 + 12y^2).
%
% Hence only the final t = 0 problem has u_exact as an exact solution;
% intermediate t-level problems are genuine continuation states.
%
% To mimic the paper's vanishing-moment boundary treatment, we use
%
%   u = g   on dOmega,      Delta u = t   on dOmega
%
% in physical coordinates. After mapping to the reference square, the
% auxiliary boundary condition becomes Delta_xi u = t / 4 on dOmega.

    clc;
    close all;

    problem = default_problem();
    result = solve_t_homotopy(problem);
    output_dir = fullfile(pwd, 'output', 'convex_t_homotopy_fixed_rhs_paper_bc');
    plot_result(problem, result, output_dir);

    fprintf('Fixed-RHS convex t-homotopy with paper-like BC finished with N = %d\n', problem.N);
    fprintf('Final continuation parameter t = %.2e\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    fprintf('Final L2 error                 = %.3e\n', result.l2_error(end));
    fprintf('Final H1 seminorm error        = %.3e\n', result.h1_error(end));
    fprintf('Final H2 seminorm error        = %.3e\n', result.h2_error(end));
    fprintf('Final Jacobian sigma_min       = %.3e\n', result.jacobian_sigma_min(end));
    fprintf('Final Jacobian cond            = %.3e\n', result.jacobian_cond(end));
    fprintf('Final Hessian lambda_min min   = %.3e\n', result.hessian_lambda_min(end));
    fprintf('Final Hessian lambda_max max   = %.3e\n', result.hessian_lambda_max(end));
    fprintf('Final Hessian PSD ratio        = %.3f\n', result.hessian_psd_ratio(end));
    fprintf('Final Hessian indefinite ratio = %.3f\n', result.hessian_indef_ratio(end));
    fprintf('Figures saved under            = %s\n', output_dir);

    disp(table(result.t_values, result.l2_error, result.h1_error, result.h2_error, ...
        'VariableNames', {'t', 'L2Error', 'H1Error', 'H2Error'}));
end

function problem = default_problem()
    problem.N = 15;
    problem.Q = 32;
    problem.eval_points = 121;
    problem.t_start = 1.00;
    problem.t_end = 0.00;
    problem.initial_step = 0.10;
    problem.min_step = 0.0025;
    problem.max_step = 0.10;
    problem.step_growth = 1.20;
    problem.accept_tol = 1e-8;
    problem.max_continuation_steps = 120;
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.line_search_min = 2^-12;
    problem.convexity_psd_drop_tol = 0.01;
    problem.convexity_indef_rise_tol = 0.01;
    problem.convexity_lam_drop_tol = 2.0;
    problem.convexity_score_tol = 1e-3;
    problem.convexity_target_psd_ratio = 0.995;
    problem.convexity_target_lam_min = 1e-8;

    problem.u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + x.^4 + y.^4;
    problem.fixed_f = @(x, y) (1 + 12 * x.^2) .* (1 + 12 * y.^2);
end

function result = solve_t_homotopy(problem)
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points);
    data = build_example_data(problem, ops);

    U_vec = zeros(problem.N^2, 1);
    t_history = zeros(problem.max_continuation_steps, 1);
    residual_history = zeros(problem.max_continuation_steps, 1);
    newton_steps = zeros(problem.max_continuation_steps, 1);
    l2_error = zeros(problem.max_continuation_steps, 1);
    h1_error = zeros(problem.max_continuation_steps, 1);
    h2_error = zeros(problem.max_continuation_steps, 1);
    jacobian_sigma_min = zeros(problem.max_continuation_steps, 1);
    jacobian_cond = zeros(problem.max_continuation_steps, 1);
    hessian_lambda_min = zeros(problem.max_continuation_steps, 1);
    hessian_lambda_max = zeros(problem.max_continuation_steps, 1);
    hessian_psd_ratio = zeros(problem.max_continuation_steps, 1);
    hessian_nsd_ratio = zeros(problem.max_continuation_steps, 1);
    hessian_indef_ratio = zeros(problem.max_continuation_steps, 1);

    current_t = problem.t_start;
    step = problem.initial_step;
    rhs_vec = rhs_for_t(current_t, problem, data, ops);
    [U_vec, iter_count, res_norm] = newton_solve( ...
        U_vec, ...
        @(z) residual(z, current_t, rhs_vec, ops, data), ...
        @(z) jacobian(z, current_t, ops, data), ...
        @(z) convexity_metrics(z, ops, data), ...
        problem, ...
        problem.newton_tol, ...
        problem.newton_maxit, ...
        problem.line_search_min);

    accepted = 1;
    t_history(accepted) = current_t;
    residual_history(accepted) = res_norm;
    newton_steps(accepted) = iter_count;
    [l2_error(accepted), h1_error(accepted), h2_error(accepted)] = compute_errors(U_vec, ops, data);
    [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, data));
    [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, data);
    jacobian_sigma_min(accepted) = sigma_min;
    jacobian_cond(accepted) = cond_val;
    hessian_lambda_min(accepted) = lam_min;
    hessian_lambda_max(accepted) = lam_max;
    hessian_psd_ratio(accepted) = psd_ratio;
    hessian_nsd_ratio(accepted) = nsd_ratio;
    hessian_indef_ratio(accepted) = indef_ratio;
    fprintf(['t = %.2e, Newton steps = %2d, residual = %.3e, sigma_min = %.3e, cond = %.3e, ' ...
        'lam_min = %.3e, lam_max = %.3e, psd = %.3f, indef = %.3f, L2 = %.3e\n'], ...
        current_t, iter_count, res_norm, sigma_min, cond_val, lam_min, lam_max, psd_ratio, indef_ratio, l2_error(accepted));

    while current_t > problem.t_end && accepted < problem.max_continuation_steps
        target_t = max(problem.t_end, current_t - step);
        rhs_vec = rhs_for_t(target_t, problem, data, ops);

        [U_try, iter_count, res_norm] = newton_solve( ...
            U_vec, ...
            @(z) residual(z, target_t, rhs_vec, ops, data), ...
            @(z) jacobian(z, target_t, ops, data), ...
            @(z) convexity_metrics(z, ops, data), ...
            problem, ...
            problem.newton_tol, ...
            problem.newton_maxit, ...
            problem.line_search_min);

        if res_norm <= problem.accept_tol
            U_vec = U_try;
            current_t = target_t;
            accepted = accepted + 1;
            t_history(accepted) = current_t;
            residual_history(accepted) = res_norm;
            newton_steps(accepted) = iter_count;
            [l2_error(accepted), h1_error(accepted), h2_error(accepted)] = compute_errors(U_vec, ops, data);
            [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, data));
            [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, data);
            jacobian_sigma_min(accepted) = sigma_min;
            jacobian_cond(accepted) = cond_val;
            hessian_lambda_min(accepted) = lam_min;
            hessian_lambda_max(accepted) = lam_max;
            hessian_psd_ratio(accepted) = psd_ratio;
            hessian_nsd_ratio(accepted) = nsd_ratio;
            hessian_indef_ratio(accepted) = indef_ratio;
            fprintf(['t = %.2e, Newton steps = %2d, residual = %.3e, step = %.4f, sigma_min = %.3e, cond = %.3e, ' ...
                'lam_min = %.3e, lam_max = %.3e, psd = %.3f, indef = %.3f, L2 = %.3e\n'], ...
                current_t, iter_count, res_norm, step, sigma_min, cond_val, lam_min, lam_max, psd_ratio, indef_ratio, l2_error(accepted));

            if current_t > problem.t_end
                step = min(problem.max_step, step * problem.step_growth);
            end
        else
            new_step = step / 2;
            fprintf('reject t = %.2e, residual = %.3e, step %.4f -> %.4f\n', target_t, res_norm, step, new_step);
            step = new_step;
            if step < problem.min_step
                warning('ContinuationFailed:MinStep', ...
                    'Adaptive continuation stopped before reaching t = 0; required step %.4g below min_step %.4g.', ...
                    step, problem.min_step);
                break;
            end
        end
    end

    result.t_values = t_history(1:accepted);
    result.residual_history = residual_history(1:accepted);
    result.newton_steps = newton_steps(1:accepted);
    result.l2_error = l2_error(1:accepted);
    result.h1_error = h1_error(1:accepted);
    result.h2_error = h2_error(1:accepted);
    result.jacobian_sigma_min = jacobian_sigma_min(1:accepted);
    result.jacobian_cond = jacobian_cond(1:accepted);
    result.hessian_lambda_min = hessian_lambda_min(1:accepted);
    result.hessian_lambda_max = hessian_lambda_max(1:accepted);
    result.hessian_psd_ratio = hessian_psd_ratio(1:accepted);
    result.hessian_nsd_ratio = hessian_nsd_ratio(1:accepted);
    result.hessian_indef_ratio = hessian_indef_ratio(1:accepted);
    result.x_eval = map_to_physical(ops.Xe);
    result.y_eval = map_to_physical(ops.Ye);
    result.u_exact_eval = data.exact_eval;
    result.u_num_eval = total_solution_on_eval(U_vec, ops, data);
    result.error_eval = result.u_num_eval - result.u_exact_eval;
end

function rhs_vec = rhs_for_t(t_val, problem, data, ops)
    f_q = problem.fixed_f(data.xq, data.yq) / 16;
    rhs_vec = ops.T1 * f_q(:) - t_val * data.lift_lap_proj + (t_val^2 / 4) * data.boundary_vec;
end

function data = build_example_data(problem, ops)
    xq = map_to_physical(ops.Xq);
    yq = map_to_physical(ops.Yq);
    xe = map_to_physical(ops.Xe);
    ye = map_to_physical(ops.Ye);

    exact_eval = problem.u_exact(xe, ye);
    Ub_q = boundary_extension_from_exact(xq, yq, problem.u_exact);
    Ub_eval = boundary_extension_from_exact(xe, ye, problem.u_exact);
    lift = differentiate_nodal_field(Ub_q(:), ops);
    lift_lap = lift.xx_vals + lift.yy_vals;

    data.xq = xq;
    data.yq = yq;
    data.lift_xx = lift.xx_vals;
    data.lift_yy = lift.yy_vals;
    data.lift_xy = lift.xy_vals;
    data.lift_lap = lift_lap;
    data.lift_lap_proj = ops.T1 * lift_lap;
    data.exact_eval = exact_eval;
    data.Ub_eval = Ub_eval;
    data.boundary_vec = build_boundary_vector(ops);
end

function vals = boundary_extension_from_exact(x, y, ufun)
    g0y = ufun(0 * y, y);
    g1y = ufun(1 + 0 * y, y);
    gx0 = ufun(x, 0 * x);
    gx1 = ufun(x, 1 + 0 * x);
    g00 = ufun(0, 0);
    g10 = ufun(1, 0);
    g01 = ufun(0, 1);
    g11 = ufun(1, 1);

    vals = (1 - x) .* g0y + x .* g1y ...
        + (1 - y) .* (gx0 - (1 - x) .* g00 - x .* g10) ...
        + y .* (gx1 - (1 - x) .* g01 - x .* g11);
end

function out = differentiate_nodal_field(field_vec, ops)
    leg_vec = ops.forward_legendre * field_vec;
    out.xx_vals = ops.nodal_from_legendre * kron(ops.Id_full, ops.B5) * leg_vec;
    out.yy_vals = ops.nodal_from_legendre * kron(ops.B5, ops.Id_full) * leg_vec;
    out.xy_vals = ops.nodal_from_legendre * kron(ops.B6, ops.B6) * leg_vec;
end

function bvec = build_boundary_vector(ops)
    dphi_left = basis_endpoint_derivative(ops.N, -1);
    dphi_right = basis_endpoint_derivative(ops.N, 1);
    int_phi = basis_integrals(ops.N);

    bmat = (dphi_right - dphi_left) * int_phi.' + int_phi * (dphi_right - dphi_left).';
    bvec = bmat(:);
end

function vals = basis_endpoint_derivative(N, endpoint)
    vals = zeros(N, 1);
    for n = 0:N-1
        vals(n + 1) = legendre_derivative_endpoint(n, endpoint) - legendre_derivative_endpoint(n + 2, endpoint);
    end
end

function val = legendre_derivative_endpoint(n, endpoint)
    base = n * (n + 1) / 2;
    if endpoint > 0
        val = base;
    else
        val = (-1)^(n + 1) * base;
    end
end

function vals = basis_integrals(N)
    vals = zeros(N, 1);
    vals(1) = 2;
end

function r = residual(U_vec, t_val, rhs_vec, ops, data)
    u = second_derivatives(U_vec, ops, data);
    biharm_vec = ops.T * U_vec;
    det_vec = ops.T1 * (u.xx .* u.yy - u.xy .^ 2);
    r = t_val * biharm_vec + (1 - t_val) * det_vec - rhs_vec;
end

function J = jacobian(U_vec, t_val, ops, data)
    u = second_derivatives(U_vec, ops, data);
    nonlin_jac = ops.T1 * ( ...
        diag(u.yy) * ops.T2 + ...
        diag(u.xx) * ops.T3 - ...
        2 * diag(u.xy) * ops.T4);
    J = t_val * ops.T + (1 - t_val) * nonlin_jac;
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

function [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, data)
    u = second_derivatives(U_vec, ops, data);
    trace_h = u.xx + u.yy;
    disc = sqrt(max((u.xx - u.yy) .^ 2 + 4 * u.xy .^ 2, 0));
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

function metrics = convexity_metrics(U_vec, ops, data)
    [lam_min, lam_max, psd_ratio, ~, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, data);
    metrics.lam_min = lam_min;
    metrics.lam_max = lam_max;
    metrics.psd_ratio = psd_ratio;
    metrics.indef_ratio = indef_ratio;
    metrics.score = 2 * indef_ratio - psd_ratio + max(0, -lam_min) / (1 + max(1, abs(lam_max)));
end

function ok = convexity_acceptable(trial_metrics, current_metrics, problem)
    strongly_convex = ...
        trial_metrics.psd_ratio >= problem.convexity_target_psd_ratio && ...
        trial_metrics.lam_min >= problem.convexity_target_lam_min;

    nonworsening = ...
        trial_metrics.psd_ratio >= current_metrics.psd_ratio - problem.convexity_psd_drop_tol && ...
        trial_metrics.indef_ratio <= current_metrics.indef_ratio + problem.convexity_indef_rise_tol && ...
        trial_metrics.lam_min >= current_metrics.lam_min - problem.convexity_lam_drop_tol;

    improving = trial_metrics.score <= current_metrics.score - problem.convexity_score_tol;
    ok = strongly_convex || nonworsening || improving;
end

function u = second_derivatives(U_vec, ops, data)
    u.xx = ops.T2 * U_vec + data.lift_xx;
    u.yy = ops.T3 * U_vec + data.lift_yy;
    u.xy = ops.T4 * U_vec + data.lift_xy;
end

function [l2e, h1e, h2e] = compute_errors(U_vec, ops, data)
    total_eval = total_solution_on_eval(U_vec, ops, data);
    err = total_eval - data.exact_eval;

    hx = 1 / (size(err, 1) - 1);
    [ex_x, ex_y] = gradient(err, hx, hx);
    [ex_xx, ex_xy] = gradient(ex_x, hx, hx);
    [~, ex_yy] = gradient(ex_y, hx, hx);

    l2e = sqrt(sum(err(:) .^ 2) * hx * hx);
    h1e = sqrt(sum((ex_x(:) .^ 2 + ex_y(:) .^ 2)) * hx * hx);
    h2e = sqrt(sum((ex_xx(:) .^ 2 + 2 * ex_xy(:) .^ 2 + ex_yy(:) .^ 2)) * hx * hx);
end

function total_eval = total_solution_on_eval(U_vec, ops, data)
    U_hat = reshape(U_vec, ops.N, ops.N);
    U_eval = ops.Phi_eval * U_hat * ops.Phi_eval.';
    total_eval = U_eval + data.Ub_eval;
end

function plot_result(problem, result, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fig1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1400, 420]);

    subplot(1, 3, 1);
    surf(result.x_eval, result.y_eval, result.u_exact_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u(x,y)');
    title('Exact Solution');
    colorbar;

    subplot(1, 3, 2);
    surf(result.x_eval, result.y_eval, result.u_num_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u_{t,N}(x,y)');
    title(sprintf('Numerical Solution, t = %.0e', result.t_values(end)));
    colorbar;

    subplot(1, 3, 3);
    surf(result.x_eval, result.y_eval, result.error_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('error');
    title('Pointwise Error');
    colorbar;

    exportgraphics(fig1, fullfile(output_dir, 'convex_t_homotopy_surfaces.png'), 'Resolution', 180);
    close(fig1);

    fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [120, 120, 700, 520]);
    semilogy(result.t_values + 1e-14, result.l2_error, '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    semilogy(result.t_values + 1e-14, result.h1_error, '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
    semilogy(result.t_values + 1e-14, result.h2_error, '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('t');
    ylabel('error');
    title(sprintf('Fixed-RHS t-Homotopy Error Curves, N = %d', problem.N));
    legend('L2', 'H1', 'H2', 'Location', 'southwest');
    exportgraphics(fig2, fullfile(output_dir, 'convex_t_homotopy_fixed_rhs_error_curves.png'), 'Resolution', 180);
    close(fig2);
end

function x = map_to_physical(xi)
    x = (xi + 1) / 2;
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
    ops.forward_legendre = forward_legendre;
    ops.nodal_from_legendre = nodal_from_legendre;
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

function [z, iter_count, res_norm] = newton_solve(z0, F, Jfun, convexity_fun, problem, tol, maxit, line_search_min)
    z = z0;
    r = F(z);
    res_norm = norm(r, inf);
    current_metrics = convexity_fun(z);

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
            trial_metrics = convexity_fun(z_trial);
            if norm(r_trial, inf) < (1 - 1e-4 * alpha) * res_norm && ...
                    convexity_acceptable(trial_metrics, current_metrics, problem)
                z = z_trial;
                r = r_trial;
                res_norm = norm(r, inf);
                current_metrics = trial_metrics;
                accepted = true;
                break;
            end
            alpha = alpha / 2;
        end

        if ~accepted
            return;
        end
    end
end
