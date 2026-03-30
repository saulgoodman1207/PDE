function semilinear_to_monge_ampere_homotopy_demo(case_name)
% Legendre-Galerkin continuation for the semilinear-to-Monge-Ampere homotopy
%
%   t * (Delta u + u^3) + (1 - t) * det(D^2 u) = f,
%
% on the physical square (0,1)^2 with Dirichlet boundary data inherited
% from the target Monge-Ampere problem. The parameter t decreases from 1 to
% 0, so the solve starts from the semilinear endpoint
%
%   Delta u + u^3 = f,
%
% and ends at
%
%   det(D^2 u) = f.
%
% We reuse the paper-style Legendre-Galerkin basis, B-matrix transform
% chain, boundary lifting, continuation strategy, and convexity diagnostics
% from the stable u0-anchored prototype, but keep the new homotopy in a
% standalone file so the stable code remains untouched.
%
% Required problem inputs for the solver are:
%   f : right-hand side function on (0,1)^2
%   g : Dirichlet boundary data
%
% Optional manufactured-data diagnostics can additionally provide
%   u_exact : exact solution used only for error curves and comparison plots.

    if nargin < 1
        case_name = 'paper_example_5_1';
    end

    clc;
    close all;

    problem = default_problem(case_name);
    result = solve_t_homotopy(problem);
    output_dir = fullfile(pwd, 'output', ['semilinear_to_monge_ampere_homotopy_' problem.case_name]);
    plot_result(problem, result, output_dir);
    export_history(result, output_dir);

    fprintf('Manufactured case              = %s\n', problem.case_name);
    fprintf('Semilinear homotopy finished with N = %d\n', problem.N);
    fprintf('Initial guess construction     = %s\n', problem.initial_guess_label);
    fprintf('Bootstrap acceptance rule      = %s\n', bootstrap_acceptance_label(problem));
    fprintf('Bootstrap convexity relaxation = %s\n', bootstrap_relaxation_label(problem));
    fprintf('Solve success                  = %s\n', logical_text(result.success));
    fprintf('Reached target t = %.2e        = %s\n', problem.t_end, logical_text(result.reached_t_end));
    fprintf('Initial guess RHS min f on quadrature = %.3e\n', result.initial_rhs_min_f);
    fprintf('Initial guess RHS clipped to zero     = %s\n', logical_text(result.initial_rhs_clipped_to_zero));
    fprintf('Final continuation parameter t = %.2e\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    fprintf('Final auxiliary residual inf   = %.3e\n', result.aux_residual_history(end));
    fprintf('Final MA residual inf          = %.3e\n', result.ma_residual_history(end));
    if problem.has_exact
        fprintf('Final grid L2 indicator        = %.3e\n', result.l2_error_indicator(end));
        fprintf('Final FD-H1 indicator          = %.3e\n', result.h1_fd_error_indicator(end));
        fprintf('Final FD-H2 indicator          = %.3e\n', result.h2_fd_error_indicator(end));
    else
        fprintf('Final error indicators         = unavailable (u_exact not supplied)\n');
    end
    if ~result.success
        fprintf('Failure reason                 = %s\n', result.failure_reason);
    end
    fprintf('Final Jacobian sigma_min       = %.3e\n', result.jacobian_sigma_min(end));
    fprintf('Final Jacobian cond            = %.3e\n', result.jacobian_cond(end));
    fprintf('Final auxiliary Jacobian sigma = %.3e\n', result.aux_jacobian_sigma_min(end));
    fprintf('Final Hessian lambda_min min   = %.3e\n', result.hessian_lambda_min(end));
    fprintf('Final Hessian lambda_max max   = %.3e\n', result.hessian_lambda_max(end));
    fprintf('Final Hessian PSD ratio        = %.3f\n', result.hessian_psd_ratio(end));
    fprintf('Final Hessian indefinite ratio = %.3f\n', result.hessian_indef_ratio(end));
    fprintf('Final det-positive ratio       = %.3f\n', result.det_positive_ratio(end));
    fprintf('Final solution min/max         = [%.3e, %.3e]\n', ...
        result.solution_min(end), result.solution_max(end));
    fprintf('Figures saved under            = %s\n', output_dir);

    if problem.has_exact
        disp(table(result.t_values, result.l2_error_indicator, ...
            result.h1_fd_error_indicator, result.h2_fd_error_indicator, ...
            'VariableNames', {'t', 'GridL2Indicator', 'FDH1Indicator', 'FDH2Indicator'}));
    end
end

function problem = default_problem(case_name)
    problem.N = 32;
    problem.Q = 32;
    problem.eval_points = 121;
    problem.t_values = [1.0, 10 .^ (-(1:2:13)), 0.0];
    problem.t_start = problem.t_values(1);
    problem.t_end = problem.t_values(end);
    problem.initial_step = 0.10;
    problem.min_step = 0.0025;
    problem.max_step = 0.10;
    problem.step_growth = 1.20;
    problem.accept_tol = 1e-8;
    problem.max_continuation_steps = 120;
    problem.newton_tol = 1e-10;
    problem.newton_maxit = 20;
    problem.line_search_min = 2^-12;
    problem.semilinear_bootstrap_values = 0.0:0.1:1.0;
    problem.semilinear_bootstrap_stage_tol = 1e-5;
    problem.bootstrap_acceptance_mode = 'semilinear_health';
    problem.bootstrap_jacobian_rcond_min = 1e-12;
    problem.bootstrap_jacobian_cond_max = 1e8;
    problem.bootstrap_jacobian_rcond_drop_factor = 0.20;
    problem.bootstrap_jacobian_cond_growth_factor = 10.0;
    problem.bootstrap_solution_max_abs_growth = 1.50;
    problem.bootstrap_solution_max_abs_margin = 0.25;
    problem.show_figures = true;
    problem.require_convexity = true;
    problem.bootstrap_relax_convexity_cases = {'quartic_convex', 'mild_sextic_convex'};
    problem.bootstrap_relax_convexity_until_s = 0.30;
    problem.convexity_psd_drop_tol = 0.01;
    problem.convexity_indef_rise_tol = 0.01;
    problem.convexity_lam_drop_tol = 2.0;
    problem.convexity_score_tol = 1e-3;
    problem.convexity_target_psd_ratio = 0.995;
    problem.convexity_target_lam_min = 1e-8;
    problem.main_convexity_schedule_t = [0.50, 0.20, 0.10, 1e-2];
    problem.main_convexity_schedule_psd_floor = [0.05, 0.25, 0.90, 0.995];
    problem.main_convexity_schedule_lam_floor = [-1.0, -0.10, -1e-3, 1e-8];
    problem.main_convexity_schedule_det_floor = [0.10, 0.50, 0.95, 0.995];
    problem.det_positive_tol = 1e-10;
    problem.f_negative_tol = 1e-12;
    problem.diagnostics = struct();
    problem.has_exact = false;
    problem = configure_case(problem, case_name);
    problem.initial_guess_label = ...
        'Poisson warm start, then bootstrap solve of \Delta u + s u^3 = f';
end

function problem = configure_case(problem, case_name)
    % All current manufactured cases use the same semilinear endpoint
    % Delta u + u^3 = f, and differ only in f, g, and optional exact data.
    switch lower(case_name)
        case {'mild_exp_separable', 'separable_exponential_convex'}
            problem.case_name = 'mild_exp_separable';
            alpha = 0.2;
            u_exact = @(x, y) exp(alpha * x) + exp(alpha * y) - 2;
            u_xx = @(x, y) alpha^2 * exp(alpha * x);
            u_yy = @(x, y) alpha^2 * exp(alpha * y);
            problem.f = @(x, y) u_xx(x, y) .* u_yy(x, y);
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case {'mixed_polynomial_convex', 'sextic_mixed_convex'}
            problem.case_name = 'mixed_polynomial_convex';
            a = 0.01;
            b = 0.01;
            u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + a * (x.^6 + y.^6) + b * (x.^2 .* y.^2);
            u_xx = @(x, y) 1 + 30 * a * x.^4 + 2 * b * y.^2;
            u_yy = @(x, y) 1 + 30 * a * y.^4 + 2 * b * x.^2;
            u_xy = @(x, y) 4 * b * x .* y;
            problem.f = @(x, y) u_xx(x, y) .* u_yy(x, y) - u_xy(x, y).^2;
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case {'gaussian_bump_convex', 'quadratic_gaussian_convex'}
            problem.case_name = 'gaussian_bump_convex';
            amplitude = 0.03;
            beta = 6.0;
            center = 0.5;
            r2 = @(x, y) (x - center).^2 + (y - center).^2;
            bump = @(x, y) exp(-beta * r2(x, y));
            u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + amplitude * bump(x, y);
            u_xx = @(x, y) 1 + amplitude * bump(x, y) .* ...
                (4 * beta^2 * (x - center).^2 - 2 * beta);
            u_yy = @(x, y) 1 + amplitude * bump(x, y) .* ...
                (4 * beta^2 * (y - center).^2 - 2 * beta);
            u_xy = @(x, y) amplitude * bump(x, y) .* ...
                (4 * beta^2 * (x - center) .* (y - center));
            problem.f = @(x, y) u_xx(x, y) .* u_yy(x, y) - u_xy(x, y).^2;
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case {'mild_sextic_convex', 'separable_sextic_convex'}
            problem.case_name = 'mild_sextic_convex';
            a = 0.02;
            u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + a * (x.^6 + y.^6);
            u_xx = @(x, y) 1 + 30 * a * x.^4;
            u_yy = @(x, y) 1 + 30 * a * y.^4;
            problem.f = @(x, y) u_xx(x, y) .* u_yy(x, y);
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case {'trig_convex', 'quadratic_trig_convex'}
            problem.case_name = 'trig_convex';
            alpha = 0.05;
            u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + ...
                alpha * (sin(pi * x).^2 + sin(pi * y).^2);
            hx = @(x) 1 + 2 * alpha * pi^2 * cos(2 * pi * x);
            hy = @(y) 1 + 2 * alpha * pi^2 * cos(2 * pi * y);
            problem.f = @(x, y) hx(x) .* hy(y);
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case {'scaled_quartic_convex', 'quartic_convex_scaled'}
            problem.case_name = 'scaled_quartic_convex';
            scale = 0.20;
            u_exact = @(x, y) scale * (0.5 * (x.^2 + y.^2) + x.^4 + y.^4);
            problem.f = @(x, y) (scale^2) * (1 + 12 * x.^2) .* (1 + 12 * y.^2);
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case 'quartic_convex'
            problem.case_name = 'quartic_convex';
            problem.f = @(x, y) (1 + 12 * x.^2) .* (1 + 12 * y.^2);
            problem.g = @(x, y) 0.5 * (x.^2 + y.^2) + x.^4 + y.^4;
            problem = attach_exact_diagnostics(problem, problem.g);
        case {'paper_example_5_1', 'exp_radial'}
            problem.case_name = 'paper_example_5_1';
            problem.f = @(x, y) (1 + x.^2 + y.^2) .* exp(x.^2 + y.^2);
            problem.g = @(x, y) exp((x.^2 + y.^2) / 2);
            problem = attach_exact_diagnostics(problem, problem.g);
        case {'paper_example_5_1_noexact', 'exp_radial_noexact'}
            problem.case_name = 'paper_example_5_1_noexact';
            problem.f = @(x, y) (1 + x.^2 + y.^2) .* exp(x.^2 + y.^2);
            problem.g = @(x, y) exp((x.^2 + y.^2) / 2);
        otherwise
            error('Unknown manufactured case "%s".', case_name);
    end
end

function problem = attach_exact_diagnostics(problem, u_exact_fun)
    problem.diagnostics.u_exact = u_exact_fun;
    problem.has_exact = true;
end

function result = solve_t_homotopy(problem)
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points);
    data = build_example_data(problem, ops);

    U_vec = data.U_init_vec;
    num_record_levels = numel(problem.t_values);
    t_history = zeros(num_record_levels, 1);
    residual_history = zeros(num_record_levels, 1);
    aux_residual_history = zeros(num_record_levels, 1);
    ma_residual_history = zeros(num_record_levels, 1);
    lap_term_history = zeros(num_record_levels, 1);
    cubic_term_history = zeros(num_record_levels, 1);
    det_term_history = zeros(num_record_levels, 1);
    newton_steps = zeros(num_record_levels, 1);
    l2_error_indicator = zeros(num_record_levels, 1);
    h1_fd_error_indicator = zeros(num_record_levels, 1);
    h2_fd_error_indicator = zeros(num_record_levels, 1);
    jacobian_sigma_min = zeros(num_record_levels, 1);
    jacobian_cond = zeros(num_record_levels, 1);
    aux_jacobian_sigma_min = zeros(num_record_levels, 1);
    hessian_lambda_min = zeros(num_record_levels, 1);
    hessian_lambda_max = zeros(num_record_levels, 1);
    hessian_psd_ratio = zeros(num_record_levels, 1);
    hessian_nsd_ratio = zeros(num_record_levels, 1);
    hessian_indef_ratio = zeros(num_record_levels, 1);
    det_positive_ratio = zeros(num_record_levels, 1);
    solution_min = zeros(num_record_levels, 1);
    solution_max = zeros(num_record_levels, 1);

    current_t = problem.t_start;
    step = problem.initial_step;
    [U_vec, iter_count, res_norm, level_status] = solve_semilinear_bootstrap(U_vec, ops, data, problem);
    res_norm = norm(residual(U_vec, current_t, ops, data), inf);
    failure_reason = '';

    accepted = 1;
    t_history(accepted) = current_t;
    residual_history(accepted) = res_norm;
    newton_steps(accepted) = iter_count;
    [aux_residual_history(accepted), ma_residual_history(accepted), ...
        lap_term_history(accepted), cubic_term_history(accepted), det_term_history(accepted)] = ...
        residual_component_diagnostics(U_vec, current_t, ops, data, problem);
    [l2_error_indicator(accepted), h1_fd_error_indicator(accepted), h2_fd_error_indicator(accepted)] = ...
        compute_error_indicators(U_vec, ops, data, problem);
    [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, data));
    jacobian_sigma_min(accepted) = sigma_min;
    jacobian_cond(accepted) = cond_val;
    [aux_sigma_min, ~] = jacobian_diagnostics(auxiliary_jacobian(U_vec, ops, data));
    aux_jacobian_sigma_min(accepted) = aux_sigma_min;
    [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio, det_pos_ratio, u_min, u_max] = ...
        hessian_branch_diagnostics(U_vec, ops, data, problem);
    hessian_lambda_min(accepted) = lam_min;
    hessian_lambda_max(accepted) = lam_max;
    hessian_psd_ratio(accepted) = psd_ratio;
    hessian_nsd_ratio(accepted) = nsd_ratio;
    hessian_indef_ratio(accepted) = indef_ratio;
    det_positive_ratio(accepted) = det_pos_ratio;
    solution_min(accepted) = u_min;
    solution_max(accepted) = u_max;
    print_continuation_status(current_t, iter_count, res_norm, [], sigma_min, cond_val, ...
        aux_residual_history(accepted), ma_residual_history(accepted), ...
        lam_min, lam_max, psd_ratio, indef_ratio, det_pos_ratio, ...
        l2_error_indicator(accepted), problem.has_exact);

    if ~level_status.converged
        failure_reason = sprintf('Initial semilinear endpoint solve at t = %.2e did not converge (%s).', ...
            current_t, level_status.reason);
    end

    next_record_index = 2;
    internal_step_count = 0;
    while level_status.converged && next_record_index <= numel(problem.t_values)
        target_record_t = problem.t_values(next_record_index);
        target_t = max(target_record_t, current_t - step);
        if abs(target_t - target_record_t) <= 10 * eps(max(1, abs(current_t)))
            target_t = target_record_t;
        end

        [U_try, iter_count, res_norm, level_status] = solve_level(U_vec, target_t, ops, data, problem);
        level_accepted = level_status.converged && res_norm <= problem.accept_tol;

        if level_accepted
            U_vec = U_try;
            current_t = target_t;
            internal_step_count = internal_step_count + 1;
            if current_t <= target_record_t + 10 * eps(max(1, abs(target_record_t)))
                accepted = accepted + 1;
                t_history(accepted) = current_t;
                residual_history(accepted) = res_norm;
                newton_steps(accepted) = iter_count;
                [aux_residual_history(accepted), ma_residual_history(accepted), ...
                    lap_term_history(accepted), cubic_term_history(accepted), det_term_history(accepted)] = ...
                    residual_component_diagnostics(U_vec, current_t, ops, data, problem);
                [l2_error_indicator(accepted), h1_fd_error_indicator(accepted), h2_fd_error_indicator(accepted)] = ...
                    compute_error_indicators(U_vec, ops, data, problem);
                [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, data));
                jacobian_sigma_min(accepted) = sigma_min;
                jacobian_cond(accepted) = cond_val;
                [aux_sigma_min, ~] = jacobian_diagnostics(auxiliary_jacobian(U_vec, ops, data));
                aux_jacobian_sigma_min(accepted) = aux_sigma_min;
                [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio, det_pos_ratio, u_min, u_max] = ...
                    hessian_branch_diagnostics(U_vec, ops, data, problem);
                hessian_lambda_min(accepted) = lam_min;
                hessian_lambda_max(accepted) = lam_max;
                hessian_psd_ratio(accepted) = psd_ratio;
                hessian_nsd_ratio(accepted) = nsd_ratio;
                hessian_indef_ratio(accepted) = indef_ratio;
                det_positive_ratio(accepted) = det_pos_ratio;
                solution_min(accepted) = u_min;
                solution_max(accepted) = u_max;
                print_continuation_status(current_t, iter_count, res_norm, step, sigma_min, cond_val, ...
                    aux_residual_history(accepted), ma_residual_history(accepted), ...
                    lam_min, lam_max, psd_ratio, indef_ratio, det_pos_ratio, ...
                    l2_error_indicator(accepted), problem.has_exact);
                next_record_index = next_record_index + 1;
            end
            if current_t > problem.t_end
                step = min(problem.max_step, step * problem.step_growth);
            end
        else
            new_step = step / 2;
            fprintf('reject t = %.2e, residual = %.3e, reason = %s, step %.4f -> %.4f\n', ...
                target_t, res_norm, level_status.reason, step, new_step);
            step = new_step;
            if step < problem.min_step
                failure_reason = sprintf(['Adaptive continuation stopped before reaching record t = %.2e; ' ...
                    'last rejected internal level t = %.2e (%s).'], ...
                    target_record_t, target_t, level_status.reason);
                warning('ContinuationFailed:MinStep', ...
                    ['Adaptive continuation stopped before reaching record t = %.2e; ' ...
                     'required step %.4g below min_step %.4g.'], ...
                    target_record_t, step, problem.min_step);
                break;
            end
        end

        if internal_step_count >= problem.max_continuation_steps
            failure_reason = sprintf('Maximum internal continuation steps (%d) reached before t = %.2e.', ...
                problem.max_continuation_steps, target_record_t);
            break;
        end
    end

    reached_t_end = abs(current_t - problem.t_end) <= 10 * eps(max(1, abs(current_t)));
    if reached_t_end
        success = true;
        failure_reason = '';
    else
        success = false;
        if isempty(failure_reason)
            if ~level_status.converged
                failure_reason = sprintf('Newton solve failed before reaching t = 0 (%s).', level_status.reason);
            elseif internal_step_count >= problem.max_continuation_steps
                failure_reason = sprintf('Maximum internal continuation steps (%d) reached before t = 0.', ...
                    problem.max_continuation_steps);
            else
                failure_reason = 'Continuation terminated before reaching t = 0.';
            end
        end
    end

    result.t_values = t_history(1:accepted);
    result.residual_history = residual_history(1:accepted);
    result.aux_residual_history = aux_residual_history(1:accepted);
    result.ma_residual_history = ma_residual_history(1:accepted);
    result.lap_term_history = lap_term_history(1:accepted);
    result.cubic_term_history = cubic_term_history(1:accepted);
    result.det_term_history = det_term_history(1:accepted);
    result.newton_steps = newton_steps(1:accepted);
    result.l2_error_indicator = l2_error_indicator(1:accepted);
    result.h1_fd_error_indicator = h1_fd_error_indicator(1:accepted);
    result.h2_fd_error_indicator = h2_fd_error_indicator(1:accepted);
    result.l2_error = result.l2_error_indicator;
    result.h1_error = result.h1_fd_error_indicator;
    result.h2_error = result.h2_fd_error_indicator;
    result.jacobian_sigma_min = jacobian_sigma_min(1:accepted);
    result.jacobian_cond = jacobian_cond(1:accepted);
    result.aux_jacobian_sigma_min = aux_jacobian_sigma_min(1:accepted);
    result.hessian_lambda_min = hessian_lambda_min(1:accepted);
    result.hessian_lambda_max = hessian_lambda_max(1:accepted);
    result.hessian_psd_ratio = hessian_psd_ratio(1:accepted);
    result.hessian_nsd_ratio = hessian_nsd_ratio(1:accepted);
    result.hessian_indef_ratio = hessian_indef_ratio(1:accepted);
    result.det_positive_ratio = det_positive_ratio(1:accepted);
    result.solution_min = solution_min(1:accepted);
    result.solution_max = solution_max(1:accepted);
    result.x_eval = map_to_physical(ops.Xe);
    result.y_eval = map_to_physical(ops.Ye);
    result.has_exact = problem.has_exact;
    result.success = success;
    result.reached_t_end = reached_t_end;
    result.failure_reason = failure_reason;
    result.initial_rhs_min_f = data.initial_rhs_min_f;
    result.initial_rhs_clipped_to_zero = data.initial_rhs_clipped_to_zero;
    result.u_exact_eval = data.exact_eval;
    result.u_init_eval = data.u_init_eval;
    result.u_num_eval = total_solution_on_eval(U_vec, ops, data);
    if problem.has_exact
        result.error_eval = result.u_num_eval - result.u_exact_eval;
    else
        result.error_eval = nan(size(result.u_num_eval));
    end
end

function [U_vec, iter_count, res_norm, status] = solve_level(U_init, t_val, ops, data, problem)
    res_norm = norm(residual(U_init, t_val, ops, data), inf);
    if res_norm <= problem.newton_tol
        U_vec = U_init;
        iter_count = 0;
        status = struct('converged', true, 'reason', 'initial_residual_below_tol');
        return;
    end

    [U_vec, iter_count, res_norm, status] = newton_solve( ...
        U_init, ...
        @(z) residual(z, t_val, ops, data), ...
        @(z) jacobian(z, t_val, ops, data), ...
        @(z) convexity_metrics(z, t_val, ops, data, problem), ...
        problem, ...
        problem.newton_tol, ...
        problem.newton_maxit, ...
        problem.line_search_min);
end

function [U_vec, total_iter_count, res_norm, status] = solve_semilinear_bootstrap(U_init, ops, data, problem)
    U_vec = U_init;
    total_iter_count = 0;
    status = struct('converged', true, 'reason', 'bootstrap_completed');

    for k = 1:numel(problem.semilinear_bootstrap_values)
        s_val = problem.semilinear_bootstrap_values(k);
        res_norm = norm(auxiliary_stage_residual(U_vec, s_val, ops, data), inf);
        if res_norm <= problem.newton_tol
            fprintf('bootstrap s = %.2f, Newton steps = %2d, residual = %.3e\n', s_val, 0, res_norm);
            continue;
        end

        stage_problem = bootstrap_stage_problem(problem, s_val);
        [U_vec, iter_count, res_norm, status] = newton_solve( ...
            U_vec, ...
            @(z) auxiliary_stage_residual(z, s_val, ops, data), ...
            @(z) auxiliary_stage_jacobian(z, s_val, ops, data), ...
            @(z) bootstrap_metrics(z, s_val, ops, data, problem), ...
            stage_problem, ...
            problem.newton_tol, ...
            max(problem.newton_maxit, 50), ...
            problem.line_search_min);

        if ~status.converged && res_norm <= problem.semilinear_bootstrap_stage_tol
            status = struct('converged', true, 'reason', 'bootstrap_stage_tol_reached');
        end
        total_iter_count = total_iter_count + iter_count;
        fprintf('bootstrap s = %.2f, Newton steps = %2d, residual = %.3e\n', s_val, iter_count, res_norm);
        if ~status.converged
            status.reason = sprintf('bootstrap_failed_at_s_%.2f_%s', s_val, status.reason);
            return;
        end
    end

    res_norm = norm(residual(U_vec, problem.t_start, ops, data), inf);
    if res_norm > problem.semilinear_bootstrap_stage_tol
        status = struct('converged', false, 'reason', 'bootstrap_finished_but_t1_residual_not_small');
    end
end

function stage_problem = bootstrap_stage_problem(problem, s_val)
    stage_problem = problem;
    stage_problem.acceptance_mode = 'convexity';
    if isfield(problem, 'bootstrap_acceptance_mode') && ...
            strcmpi(problem.bootstrap_acceptance_mode, 'semilinear_health')
        stage_problem.acceptance_mode = 'bootstrap_semilinear_health';
        stage_problem.require_convexity = false;
        return;
    end

    relax_this_case = any(strcmpi(problem.case_name, problem.bootstrap_relax_convexity_cases));
    if relax_this_case && s_val <= problem.bootstrap_relax_convexity_until_s + 10 * eps
        % Early semilinear bootstrap stages are not yet close to the
        % convex Monge-Ampere branch, so enforcing convexity too
        % aggressively can block even the linear Delta u = f solve.
        stage_problem.require_convexity = false;
    end
end

function data = build_example_data(problem, ops)
    xq = map_to_physical(ops.Xq);
    yq = map_to_physical(ops.Yq);
    xe = map_to_physical(ops.Xe);
    ye = map_to_physical(ops.Ye);

    f_phys_q = problem.f(xq, yq);
    min_f_phys = min(f_phys_q(:));
    if min_f_phys < -problem.f_negative_tol
        error('InitialData:NegativeRHS', ...
            ['Poisson warm start requires nonnegative f on the quadrature grid. ' ...
             'Detected min(f) = %.3e below tolerance %.3e.'], ...
            min_f_phys, problem.f_negative_tol);
    end
    clipped_to_zero = any(f_phys_q(:) < 0);
    if clipped_to_zero
        warning('InitialData:ClippedRHS', ...
            ['Small negative f values detected on the quadrature grid (min(f) = %.3e). ' ...
             'Clipping them to zero before sqrt in the warm-start solve.'], min_f_phys);
    end
    sqrt_f_phys = sqrt(max(f_phys_q(:), 0));
    f_q = f_phys_q / 16;
    Ub_q = boundary_extension_from_fun(xq, yq, problem.g);
    Ub_eval = boundary_extension_from_fun(xe, ye, problem.g);

    lift = differentiate_nodal_field(Ub_q(:), ops);
    poisson_source_vec = 0.5 * sqrt_f_phys - (lift.xx_vals + lift.yy_vals);
    U_init_vec = solve_initial_poisson(ops, poisson_source_vec);
    U_init_hat = reshape(U_init_vec, ops.N, ops.N);
    u_init_eval = ops.Phi_eval * U_init_hat * ops.Phi_eval.' + Ub_eval;

    data.f_proj = ops.T1 * f_q(:);
    data.Ub_q_vec = Ub_q(:);
    data.U_init_vec = U_init_vec;
    data.initial_rhs_min_f = min_f_phys;
    data.initial_rhs_clipped_to_zero = clipped_to_zero;
    data.lift_xx = lift.xx_vals;
    data.lift_yy = lift.yy_vals;
    data.lift_xy = lift.xy_vals;
    if problem.has_exact
        data.exact_eval = problem.diagnostics.u_exact(xe, ye);
    else
        data.exact_eval = nan(size(xe));
    end
    data.u_init_eval = u_init_eval;
    data.Ub_eval = Ub_eval;
end

function U_init_vec = solve_initial_poisson(ops, source_vec)
    U_init_vec = (ops.L + 1e-12 * speye(size(ops.L))) \ (ops.T1 * source_vec);
end

function vals = boundary_extension_from_fun(x, y, ufun)
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

function [aux_vec, ma_vec, lap_vec, cubic_vec, det_vec, u_vals] = residual_components(U_vec, ops, data)
    u = second_derivatives(U_vec, ops, data);
    u_vals = ops.T0 * U_vec + data.Ub_q_vec;
    lap_vec = ops.T1 * (u.xx + u.yy);
    cubic_vec = ops.T1 * (u_vals .^ 3);
    det_vec = ops.T1 * (u.xx .* u.yy - u.xy .^ 2);
    aux_vec = 0.25 * lap_vec + (1 / 16) * cubic_vec - data.f_proj;
    ma_vec = det_vec - data.f_proj;
end

function r = auxiliary_stage_residual(U_vec, s_val, ops, data)
    [~, ~, lap_vec, cubic_vec] = residual_components(U_vec, ops, data);
    r = 0.25 * lap_vec + (s_val / 16) * cubic_vec - data.f_proj;
end

function r = residual(U_vec, t_val, ops, data)
    [aux_vec, ma_vec] = residual_components(U_vec, ops, data);
    r = t_val * aux_vec + (1 - t_val) * ma_vec;
end

function J_aux = auxiliary_jacobian(U_vec, ops, data)
    u_vals = ops.T0 * U_vec + data.Ub_q_vec;
    cubic_jac = ops.T1 * diag(3 * u_vals .^ 2) * ops.T0;
    J_aux = 0.25 * ops.L + (1 / 16) * cubic_jac;
end

function J_stage = auxiliary_stage_jacobian(U_vec, s_val, ops, data)
    u_vals = ops.T0 * U_vec + data.Ub_q_vec;
    cubic_jac = ops.T1 * diag(3 * u_vals .^ 2) * ops.T0;
    J_stage = 0.25 * ops.L + (s_val / 16) * cubic_jac;
end

function J = jacobian(U_vec, t_val, ops, data)
    u = second_derivatives(U_vec, ops, data);
    det_jac = ops.T1 * ( ...
        diag(u.yy) * ops.T2 + ...
        diag(u.xx) * ops.T3 - ...
        2 * diag(u.xy) * ops.T4);
    J = t_val * auxiliary_jacobian(U_vec, ops, data) + (1 - t_val) * det_jac;
end

function [aux_res, ma_res, lap_norm, cubic_norm, det_norm] = residual_component_diagnostics(U_vec, t_val, ops, data, problem)
    [aux_vec, ma_vec, lap_vec, cubic_vec, det_vec] = residual_components(U_vec, ops, data);
    aux_res = norm(aux_vec, inf);
    ma_res = norm(ma_vec, inf);
    lap_norm = norm((t_val / 4) * lap_vec, inf);
    cubic_norm = norm((t_val / 16) * cubic_vec, inf);
    det_norm = norm((1 - t_val) * det_vec, inf);
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

function [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio, det_pos_ratio, u_min, u_max] = ...
        hessian_branch_diagnostics(U_vec, ops, data, problem)
    u = second_derivatives(U_vec, ops, data);
    u_vals = ops.T0 * U_vec + data.Ub_q_vec;
    trace_h = u.xx + u.yy;
    disc = sqrt(max((u.xx - u.yy) .^ 2 + 4 * u.xy .^ 2, 0));
    lambda1 = 0.5 * (trace_h - disc);
    lambda2 = 0.5 * (trace_h + disc);
    det_h = u.xx .* u.yy - u.xy .^ 2;

    tol = 1e-10;
    psd_mask = (lambda1 >= -tol) & (lambda2 >= -tol);
    nsd_mask = (lambda1 <= tol) & (lambda2 <= tol);
    indef_mask = ~(psd_mask | nsd_mask);

    lam_min = min(lambda1);
    lam_max = max(lambda2);
    psd_ratio = mean(psd_mask);
    nsd_ratio = mean(nsd_mask);
    indef_ratio = mean(indef_mask);
    det_pos_ratio = mean(det_h >= -problem.det_positive_tol);
    u_min = min(u_vals);
    u_max = max(u_vals);
end

function metrics = convexity_metrics(U_vec, t_val, ops, data, problem)
    [lam_min, lam_max, psd_ratio, ~, indef_ratio, det_pos_ratio] = ...
        hessian_branch_diagnostics(U_vec, ops, data, problem);
    stage_targets = main_convexity_targets(problem, t_val);
    metrics.lam_min = lam_min;
    metrics.lam_max = lam_max;
    metrics.psd_ratio = psd_ratio;
    metrics.indef_ratio = indef_ratio;
    metrics.det_pos_ratio = det_pos_ratio;
    metrics.score = 2 * indef_ratio - psd_ratio + max(0, -lam_min) / (1 + max(1, abs(lam_max)));
    metrics.t_val = t_val;
    metrics.stage_psd_floor = stage_targets.psd_floor;
    metrics.stage_lam_floor = stage_targets.lam_floor;
    metrics.stage_det_floor = stage_targets.det_floor;
end

function metrics = bootstrap_metrics(U_vec, s_val, ops, data, problem)
    J = auxiliary_stage_jacobian(U_vec, s_val, ops, data);
    J_dense = full(J + 1e-12 * speye(size(J)));
    jacobian_rcond = rcond(J_dense);
    if jacobian_rcond == 0
        jacobian_cond = inf;
    else
        jacobian_cond = 1 / jacobian_rcond;
    end

    u_vals = ops.T0 * U_vec + data.Ub_q_vec;
    metrics.jacobian_rcond = jacobian_rcond;
    metrics.jacobian_cond = jacobian_cond;
    metrics.max_abs_u = max(abs(u_vals));
    metrics.lam_min = nan;
    metrics.lam_max = nan;
    metrics.psd_ratio = nan;
    metrics.indef_ratio = nan;
    metrics.det_pos_ratio = nan;
    metrics.score = nan;
end

function ok = convexity_acceptable(trial_metrics, current_metrics, problem)
    stage_guard = ...
        trial_metrics.psd_ratio >= trial_metrics.stage_psd_floor && ...
        trial_metrics.lam_min >= trial_metrics.stage_lam_floor && ...
        trial_metrics.det_pos_ratio >= trial_metrics.stage_det_floor;

    strongly_convex = ...
        trial_metrics.psd_ratio >= problem.convexity_target_psd_ratio && ...
        trial_metrics.lam_min >= problem.convexity_target_lam_min;

    nonworsening = ...
        trial_metrics.psd_ratio >= current_metrics.psd_ratio - problem.convexity_psd_drop_tol && ...
        trial_metrics.indef_ratio <= current_metrics.indef_ratio + problem.convexity_indef_rise_tol && ...
        trial_metrics.lam_min >= current_metrics.lam_min - problem.convexity_lam_drop_tol;

    improving = trial_metrics.score <= current_metrics.score - problem.convexity_score_tol;
    ok = stage_guard && (strongly_convex || nonworsening || improving);
end

function ok = bootstrap_health_acceptable(trial_metrics, current_metrics, problem)
    healthy_absolute = ...
        trial_metrics.jacobian_rcond >= problem.bootstrap_jacobian_rcond_min && ...
        trial_metrics.jacobian_cond <= problem.bootstrap_jacobian_cond_max;

    healthy_relative = ...
        trial_metrics.jacobian_rcond >= ...
            problem.bootstrap_jacobian_rcond_drop_factor * current_metrics.jacobian_rcond && ...
        trial_metrics.jacobian_cond <= ...
            problem.bootstrap_jacobian_cond_growth_factor * max(current_metrics.jacobian_cond, 1);

    solution_bounded = ...
        trial_metrics.max_abs_u <= ...
        problem.bootstrap_solution_max_abs_growth * current_metrics.max_abs_u + ...
        problem.bootstrap_solution_max_abs_margin;

    ok = solution_bounded && (healthy_absolute || healthy_relative);
end

function targets = main_convexity_targets(problem, t_val)
    targets.psd_floor = 0;
    targets.lam_floor = -inf;
    targets.det_floor = 0;

    schedule_t = problem.main_convexity_schedule_t;
    for k = 1:numel(schedule_t)
        if t_val <= schedule_t(k)
            targets.psd_floor = problem.main_convexity_schedule_psd_floor(k);
            targets.lam_floor = problem.main_convexity_schedule_lam_floor(k);
            targets.det_floor = problem.main_convexity_schedule_det_floor(k);
        end
    end
end

function u = second_derivatives(U_vec, ops, data)
    u.xx = ops.T2 * U_vec + data.lift_xx;
    u.yy = ops.T3 * U_vec + data.lift_yy;
    u.xy = ops.T4 * U_vec + data.lift_xy;
end

function [l2e, h1e, h2e] = compute_error_indicators(U_vec, ops, data, problem)
    if ~problem.has_exact
        l2e = nan;
        h1e = nan;
        h2e = nan;
        return;
    end

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

function print_continuation_status(t_val, iter_count, res_norm, step, sigma_min, cond_val, ...
        aux_res, ma_res, lam_min, lam_max, psd_ratio, indef_ratio, det_pos_ratio, l2_error, has_exact)
    if isempty(step)
        prefix = 't = %.2e, Newton steps = %2d, residual = %.3e, sigma_min = %.3e, cond = %.3e, ';
        values = {t_val, iter_count, res_norm, sigma_min, cond_val};
    else
        prefix = ['t = %.2e, Newton steps = %2d, residual = %.3e, step = %.4f, ' ...
            'sigma_min = %.3e, cond = %.3e, '];
        values = {t_val, iter_count, res_norm, step, sigma_min, cond_val};
    end

    if has_exact
        suffix = ['aux = %.3e, ma = %.3e, lam_min = %.3e, lam_max = %.3e, ' ...
            'psd = %.3f, indef = %.3f, det+ = %.3f, GridL2 = %.3e\n'];
        values = [values, {aux_res, ma_res, lam_min, lam_max, psd_ratio, indef_ratio, det_pos_ratio, l2_error}];
    else
        suffix = ['aux = %.3e, ma = %.3e, lam_min = %.3e, lam_max = %.3e, ' ...
            'psd = %.3f, indef = %.3f, det+ = %.3f\n'];
        values = [values, {aux_res, ma_res, lam_min, lam_max, psd_ratio, indef_ratio, det_pos_ratio}];
    end

    fprintf([prefix suffix], values{:});
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

    fig1 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [100, 100, 1400, 820]);

    subplot(2, 2, 1);
    surf(result.x_eval, result.y_eval, result.u_init_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u_{init}(x,y)');
    title('Poisson Warm Start');
    colorbar;

    subplot(2, 2, 2);
    surf(result.x_eval, result.y_eval, result.u_num_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u_{t,N}(x,y)');
    title(sprintf('Numerical Solution, t = %.0e', result.t_values(end)));
    colorbar;

    if problem.has_exact
        subplot(2, 2, 3);
        surf(result.x_eval, result.y_eval, result.u_exact_eval, 'EdgeColor', 'none');
        view(45, 30);
        xlabel('x');
        ylabel('y');
        zlabel('u(x,y)');
        title('Exact Target');
        colorbar;

        subplot(2, 2, 4);
        surf(result.x_eval, result.y_eval, result.error_eval, 'EdgeColor', 'none');
        view(45, 30);
        xlabel('x');
        ylabel('y');
        zlabel('error');
        title('Pointwise Error');
        colorbar;
    else
        subplot(2, 2, 3);
        surf(result.x_eval, result.y_eval, result.u_num_eval - result.u_init_eval, 'EdgeColor', 'none');
        view(45, 30);
        xlabel('x');
        ylabel('y');
        zlabel('u-u_{init}');
        title('Homotopy Correction');
        colorbar;
    end

    hide_axes_toolbars(fig1);
    exportgraphics(fig1, fullfile(output_dir, 'semilinear_homotopy_surfaces.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig1);
    end

    fig2 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [120, 120, 760, 520]);
    positive_mask = result.t_values > 0;
    neg_log10_t_values = -log10(result.t_values(positive_mask));
    if problem.has_exact
        semilogy(neg_log10_t_values, result.l2_error(positive_mask), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
        hold on;
        semilogy(neg_log10_t_values, result.h1_error(positive_mask), '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
        semilogy(neg_log10_t_values, result.h2_error(positive_mask), '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
        grid on;
        xlabel('-log_{10}(t)');
        ylabel('indicator');
        title(sprintf('Semilinear Homotopy Error Indicators vs -log_{10}(t), N = %d', problem.N));
        legend('Grid L2', 'FD-H1', 'FD-H2', 'Location', 'southwest');
    else
        semilogy(neg_log10_t_values, max(result.residual_history(positive_mask), eps), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
        grid on;
        xlabel('-log_{10}(t)');
        ylabel('residual');
        title(sprintf('Semilinear Homotopy Residual Curves vs -log_{10}(t), N = %d', problem.N));
    end
    hide_axes_toolbars(fig2);
    exportgraphics(fig2, fullfile(output_dir, 'semilinear_homotopy_error_curves.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig2);
    end

    fig3 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [140, 140, 920, 760]);

    subplot(2, 1, 1);
    semilogy(result.t_values + 1e-14, max(result.residual_history, eps), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    semilogy(result.t_values + 1e-14, max(result.aux_residual_history, eps), '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
    semilogy(result.t_values + 1e-14, max(result.ma_residual_history, eps), '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('t');
    ylabel('residual');
    title('Residual Components vs t');
    legend('Total', 'Auxiliary', 'Monge-Ampere', 'Location', 'southwest');

    subplot(2, 1, 2);
    yyaxis left;
    semilogy(result.t_values + 1e-14, max(result.residual_history, eps), '-^', 'LineWidth', 1.6, 'MarkerSize', 7);
    ylabel('residual');
    yyaxis right;
    plot(result.t_values, result.newton_steps, '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    ylabel('Newton steps');
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('t');
    title('Solver Diagnostics vs t');

    hide_axes_toolbars(fig3);
    exportgraphics(fig3, fullfile(output_dir, 'semilinear_homotopy_solver_diagnostics.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig3);
    end

    fig4 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [160, 160, 920, 760]);

    subplot(2, 1, 1);
    semilogy(result.t_values + 1e-14, max(result.lap_term_history, eps), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    semilogy(result.t_values + 1e-14, max(result.cubic_term_history, eps), '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
    semilogy(result.t_values + 1e-14, max(result.det_term_history, eps), '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('t');
    ylabel('weighted term norm');
    title('Term Contributions vs t');
    legend('t/4 Laplacian', 't/16 cubic', '(1-t) determinant', 'Location', 'southwest');

    subplot(2, 1, 2);
    yyaxis left;
    semilogy(result.t_values + 1e-14, max(result.jacobian_sigma_min, eps), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    semilogy(result.t_values + 1e-14, max(result.aux_jacobian_sigma_min, eps), '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
    ylabel('\sigma_{min}');
    yyaxis right;
    plot(result.t_values, result.hessian_psd_ratio, '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
    hold on;
    plot(result.t_values, result.det_positive_ratio, '-^', 'LineWidth', 1.6, 'MarkerSize', 7);
    ylabel('ratio');
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('t');
    title('Jacobian and Branch Diagnostics vs t');
    legend('J sigma_{min}', 'J_{aux} sigma_{min}', 'PSD ratio', 'det-positive ratio', ...
        'Location', 'southwest');

    hide_axes_toolbars(fig4);
    exportgraphics(fig4, fullfile(output_dir, 'semilinear_homotopy_branch_diagnostics.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig4);
    end
    drawnow;
end

function export_history(result, output_dir)
    if result.has_exact
        history_table = table( ...
            result.t_values, ...
            result.residual_history, ...
            result.aux_residual_history, ...
            result.ma_residual_history, ...
            result.lap_term_history, ...
            result.cubic_term_history, ...
            result.det_term_history, ...
            result.newton_steps, ...
            result.l2_error_indicator, ...
            result.h1_fd_error_indicator, ...
            result.h2_fd_error_indicator, ...
            result.jacobian_sigma_min, ...
            result.jacobian_cond, ...
            result.aux_jacobian_sigma_min, ...
            result.hessian_lambda_min, ...
            result.hessian_lambda_max, ...
            result.hessian_psd_ratio, ...
            result.hessian_indef_ratio, ...
            result.det_positive_ratio, ...
            result.solution_min, ...
            result.solution_max, ...
            'VariableNames', { ...
                't', ...
                'ResidualInf', ...
                'AuxResidualInf', ...
                'MAResidualInf', ...
                'LapTermInf', ...
                'CubicTermInf', ...
                'DetTermInf', ...
                'NewtonSteps', ...
                'GridL2Indicator', ...
                'FDH1Indicator', ...
                'FDH2Indicator', ...
                'JacobianSigmaMin', ...
                'JacobianCond', ...
                'AuxJacobianSigmaMin', ...
                'HessianLambdaMin', ...
                'HessianLambdaMax', ...
                'HessianPSDRatio', ...
                'HessianIndefRatio', ...
                'DetPositiveRatio', ...
                'SolutionMin', ...
                'SolutionMax'});
    else
        history_table = table( ...
            result.t_values, ...
            result.residual_history, ...
            result.aux_residual_history, ...
            result.ma_residual_history, ...
            result.lap_term_history, ...
            result.cubic_term_history, ...
            result.det_term_history, ...
            result.newton_steps, ...
            result.jacobian_sigma_min, ...
            result.jacobian_cond, ...
            result.aux_jacobian_sigma_min, ...
            result.hessian_lambda_min, ...
            result.hessian_lambda_max, ...
            result.hessian_psd_ratio, ...
            result.hessian_indef_ratio, ...
            result.det_positive_ratio, ...
            result.solution_min, ...
            result.solution_max, ...
            'VariableNames', { ...
                't', ...
                'ResidualInf', ...
                'AuxResidualInf', ...
                'MAResidualInf', ...
                'LapTermInf', ...
                'CubicTermInf', ...
                'DetTermInf', ...
                'NewtonSteps', ...
                'JacobianSigmaMin', ...
                'JacobianCond', ...
                'AuxJacobianSigmaMin', ...
                'HessianLambdaMin', ...
                'HessianLambdaMax', ...
                'HessianPSDRatio', ...
                'HessianIndefRatio', ...
                'DetPositiveRatio', ...
                'SolutionMin', ...
                'SolutionMax'});
    end

    writetable(history_table, fullfile(output_dir, 'semilinear_homotopy_history.csv'));
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

    Phi_q = B4 * B3;
    coeff_to_legendre = kron(B3, B3);
    nodal_from_legendre = kron(B4, B4);
    forward_legendre = kron(B1, B1);
    Id_full = speye(full_count);
    Pphi = B3.' * B2 * B1;

    T0 = kron(Phi_q, Phi_q);
    T1 = kron(Pphi, Pphi);
    T2 = nodal_from_legendre * kron(Id_full, B5) * coeff_to_legendre;
    T3 = nodal_from_legendre * kron(B5, Id_full) * coeff_to_legendre;
    T4 = nodal_from_legendre * kron(B6, B6) * coeff_to_legendre;
    M = T1 * T0;
    L = T1 * (T2 + T3);

    xe = linspace(-1, 1, eval_points).';
    Pe = legendre_values(xe, full_deg);
    Phi_eval = Pe * B3;
    [Xq, Yq] = meshgrid(xq, xq);
    [Xe, Ye] = meshgrid(xe, xe);

    ops.N = N;
    ops.B5 = B5;
    ops.B6 = B6;
    ops.T0 = T0;
    ops.T1 = T1;
    ops.T2 = T2;
    ops.T3 = T3;
    ops.T4 = T4;
    ops.M = M;
    ops.L = L;
    ops.Phi_eval = Phi_eval;
    ops.Xq = Xq;
    ops.Yq = Yq;
    ops.Xe = Xe;
    ops.Ye = Ye;
    ops.forward_legendre = forward_legendre;
    ops.nodal_from_legendre = nodal_from_legendre;
    ops.Id_full = Id_full;
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

function [z, iter_count, res_norm, status] = newton_solve(z0, F, Jfun, metrics_fun, problem, tol, maxit, line_search_min)
    z = z0;
    r = F(z);
    res_norm = norm(r, inf);
    if res_norm < tol
        iter_count = 0;
        status = struct('converged', true, 'reason', 'initial_residual_below_tol');
        return;
    end

    current_metrics = metrics_fun(z);

    for iter_count = 1:maxit
        J = Jfun(z);
        delta = -(J + 1e-12 * eye(size(J))) \ r;

        alpha = 1.0;
        accepted = false;
        while alpha >= line_search_min
            z_trial = z + alpha * delta;
            r_trial = F(z_trial);
            trial_metrics = metrics_fun(z_trial);
            acceptance_ok = true;
            if isfield(problem, 'acceptance_mode') && ...
                    strcmpi(problem.acceptance_mode, 'bootstrap_semilinear_health')
                acceptance_ok = bootstrap_health_acceptable(trial_metrics, current_metrics, problem);
            elseif isfield(problem, 'require_convexity') && problem.require_convexity
                acceptance_ok = convexity_acceptable(trial_metrics, current_metrics, problem);
            end
            if norm(r_trial, inf) < (1 - 1e-4 * alpha) * res_norm && acceptance_ok
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
            status = struct('converged', false, 'reason', 'line_search_failed');
            return;
        end

        if res_norm < tol
            status = struct('converged', true, 'reason', 'newton_tol_reached');
            return;
        end
    end

    status = struct('converged', false, 'reason', 'max_iterations_reached');
end

function vis = figure_visibility(problem)
    if problem.show_figures
        vis = 'on';
    else
        vis = 'off';
    end
end

function hide_axes_toolbars(fig)
    axes_list = findall(fig, 'Type', 'axes');
    for k = 1:numel(axes_list)
        if isprop(axes_list(k), 'Toolbar') && ~isempty(axes_list(k).Toolbar)
            axes_list(k).Toolbar.Visible = 'off';
        end
    end
end

function txt = logical_text(tf)
    if tf
        txt = 'true';
    else
        txt = 'false';
    end
end

function txt = bootstrap_relaxation_label(problem)
    if isfield(problem, 'bootstrap_acceptance_mode') && ...
            strcmpi(problem.bootstrap_acceptance_mode, 'semilinear_health')
        txt = 'handled by semilinear-health acceptance rule';
        return;
    end

    if isempty(problem.bootstrap_relax_convexity_cases)
        txt = 'disabled';
        return;
    end

    if any(strcmpi(problem.case_name, problem.bootstrap_relax_convexity_cases))
        txt = sprintf('enabled for %s up to s <= %.2f', ...
            problem.case_name, problem.bootstrap_relax_convexity_until_s);
    else
        txt = sprintf('disabled for %s', problem.case_name);
    end
end

function txt = bootstrap_acceptance_label(problem)
    if isfield(problem, 'bootstrap_acceptance_mode') && ...
            strcmpi(problem.bootstrap_acceptance_mode, 'semilinear_health')
        txt = 'semilinear residual + Jacobian health';
    else
        txt = 'legacy convexity-aware line search';
    end
end
