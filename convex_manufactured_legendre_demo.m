function convex_manufactured_legendre_demo()
% Convex manufactured-solution verification on the reference square [-1,1]^2.
%
% Physical domain: (x,y) in (0,1)^2 with x = (xi+1)/2, y = (eta+1)/2.
% Exact convex solution:
%   u(x,y) = (x^2 + y^2)/2 + x^4 + y^4
% Hessian:
%   D^2 u = diag(1 + 12x^2, 1 + 12y^2) > 0
% Source term:
%   f(x,y) = det(D^2 u) = (1 + 12x^2)(1 + 12y^2)
%
% We solve the transformed vanishing-moment problem
%   -eps * Delta_xi^2 U + det(D_xi^2 (U + Ub)) = f_hat/16
% with boundary lifting Ub built from the Dirichlet data g = u|_{dOmega}.

    clc;
    close all;

    problem = default_problem();
    result = solve_example(problem);
    output_dir = fullfile(pwd, 'output', 'convex_manufactured');
    plot_example(problem, result, output_dir);

    fprintf('Convex manufactured test finished with N = %d\n', problem.N);
    fprintf('Smallest epsilon solved        = %.1e\n', result.eps_values(end));
    fprintf('Final L2 error                 = %.3e\n', result.l2_error(end));
    fprintf('Final H1 seminorm error        = %.3e\n', result.h1_error(end));
    fprintf('Final H2 seminorm error        = %.3e\n', result.h2_error(end));
    fprintf('Figures saved under            = %s\n', output_dir);

    disp(table(result.eps_values, result.l2_error, result.h1_error, result.h2_error, ...
        'VariableNames', {'epsilon', 'L2Error', 'H1Error', 'H2Error'}));
end

function problem = default_problem()
    problem.N = 15;
    problem.Q = 32;
    problem.eval_points = 121;
    problem.eps_values = 10 .^ (-(0:6));
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.line_search_min = 2^-12;

    problem.u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + x.^4 + y.^4;
    problem.f = @(x, y) (1 + 12 * x.^2) .* (1 + 12 * y.^2);
end

function result = solve_example(problem)
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points);
    data = build_example_data(problem, ops);

    U_vec = zeros(problem.N^2, 1);
    num_eps = numel(problem.eps_values);
    l2_error = zeros(num_eps, 1);
    h1_error = zeros(num_eps, 1);
    h2_error = zeros(num_eps, 1);

    for idx = 1:num_eps
        eps_val = problem.eps_values(idx);
        rhs_vec = data.f_proj + eps_val * data.lift_lap_proj - (eps_val^2 / 4) * data.boundary_vec;

        [U_vec, iter_count, res_norm] = newton_solve( ...
            U_vec, ...
            @(z) residual(z, eps_val, rhs_vec, ops, data), ...
            @(z) jacobian(z, eps_val, ops, data), ...
            problem.newton_tol, ...
            problem.newton_maxit, ...
            problem.line_search_min);

        [l2_error(idx), h1_error(idx), h2_error(idx)] = compute_errors(U_vec, ops, data);
        fprintf('eps = %.1e, Newton steps = %2d, residual = %.3e, L2 = %.3e, H1 = %.3e, H2 = %.3e\n', ...
            eps_val, iter_count, res_norm, l2_error(idx), h1_error(idx), h2_error(idx));
    end

    result.eps_values = problem.eps_values(:);
    result.l2_error = l2_error;
    result.h1_error = h1_error;
    result.h2_error = h2_error;
    result.x_eval = map_to_physical(ops.Xe);
    result.y_eval = map_to_physical(ops.Ye);
    result.u_exact_eval = data.exact_eval;
    result.u_num_eval = total_solution_on_eval(U_vec, ops, data);
    result.error_eval = result.u_num_eval - result.u_exact_eval;
end

function data = build_example_data(problem, ops)
    xq = map_to_physical(ops.Xq);
    yq = map_to_physical(ops.Yq);
    xe = map_to_physical(ops.Xe);
    ye = map_to_physical(ops.Ye);

    exact_q = problem.u_exact(xq, yq);
    exact_eval = problem.u_exact(xe, ye);
    f_q = problem.f(xq, yq) / 16;

    Ub_q = boundary_extension_from_exact(xq, yq, problem.u_exact);
    Ub_eval = boundary_extension_from_exact(xe, ye, problem.u_exact);

    lift = differentiate_nodal_field(Ub_q(:), ops);

    data.f_proj = ops.T1 * f_q(:);
    data.lift_xx = lift.xx_vals;
    data.lift_yy = lift.yy_vals;
    data.lift_xy = lift.xy_vals;
    data.lift_lap = lift.xx_vals + lift.yy_vals;
    data.lift_lap_proj = ops.T1 * data.lift_lap;
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

function r = residual(U_vec, eps_val, rhs_vec, ops, data)
    u = second_derivatives(U_vec, ops, data);
    det_vec = ops.T1 * (u.xx .* u.yy - u.xy .^ 2);
    r = -eps_val * (ops.T * U_vec) + det_vec - rhs_vec;
end

function J = jacobian(U_vec, eps_val, ops, data)
    u = second_derivatives(U_vec, ops, data);
    nonlin_jac = ops.T1 * ( ...
        diag(u.yy) * ops.T2 + ...
        diag(u.xx) * ops.T3 - ...
        2 * diag(u.xy) * ops.T4);
    J = -eps_val * ops.T + nonlin_jac;
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

function plot_example(problem, result, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    x_eval = result.x_eval;
    y_eval = result.y_eval;
    u_num = result.u_num_eval;
    u_exact = result.u_exact_eval;
    u_err = result.error_eval;

    fig1 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1400, 420]);

    subplot(1, 3, 1);
    surf(x_eval, y_eval, u_exact, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u(x,y)');
    title('Exact Solution');
    colorbar;

    subplot(1, 3, 2);
    surf(x_eval, y_eval, u_num, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u_{\epsilon,N}(x,y)');
    title(sprintf('Numerical Solution, \\epsilon = %.0e', result.eps_values(end)));
    colorbar;

    subplot(1, 3, 3);
    surf(x_eval, y_eval, u_err, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('error');
    title('Pointwise Error');
    colorbar;

    exportgraphics(fig1, fullfile(output_dir, 'convex_manufactured_surfaces.png'), 'Resolution', 180);
    close(fig1);

    fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [120, 120, 700, 520]);
    semilogy(result.eps_values, result.l2_error, '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    semilogy(result.eps_values, result.h1_error, '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
    semilogy(result.eps_values, result.h2_error, '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('\epsilon');
    ylabel('error');
    title(sprintf('Convex Manufactured Error Curves, N = %d', problem.N));
    legend('L2', 'H1', 'H2', 'Location', 'southwest');
    exportgraphics(fig2, fullfile(output_dir, 'convex_manufactured_error_curves.png'), 'Resolution', 180);
    close(fig2);
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
