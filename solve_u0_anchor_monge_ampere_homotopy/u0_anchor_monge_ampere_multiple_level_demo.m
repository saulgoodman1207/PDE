function output = u0_anchor_monge_ampere_multiple_level_demo(case_name, options)
% Legendre-Galerkin multiple-level framework for the anchored homotopy
%
%   t * (u - u0) + (1 - t) * (det(D^2 u) - f) = 0,
%
% on the physical square (0,1)^2 with Dirichlet boundary data inherited
% from the target Monge-Ampere problem. The parameter t decreases from 1 to
% 0, so the solve starts from the anchor state u = u0 and ends at
%
%   det(D^2 u) = f.
%
% We keep the mature u0-anchor discretization and Newton machinery, but
% organize the solve in a multiple-level framework:
%   (N_1, t_1) -> (N_2, t_2) -> ... -> (N_L, t_L),
% with t_{k+1} < t_k and the solution from level k projected into the
% larger N_{k+1} space as the initial guess for level k+1.
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
    if nargin < 2 || isempty(options)
        options = struct();
    end
    clc;
    close all;

    base_problem = default_problem(case_name);
    if has_n_sweep_request(options)
        output = run_multiple_level_framework(base_problem, options);
        return;
    end

    problem = apply_runtime_options(base_problem, options);
    output_dir = single_run_output_dir(problem.case_name, options);
    result = run_single_case(problem, output_dir, struct());
    print_run_summary(problem, result, output_dir);
    output = result;
end

function result = run_single_case(problem, output_dir, initial_state)
    if nargin < 3
        initial_state = struct();
    end
    result = solve_t_homotopy(problem, initial_state);
    plot_result(problem, result, output_dir);
    export_history(result, output_dir);
    export_terminal_state(result, problem, output_dir);
end

function print_run_summary(problem, result, output_dir)
    fprintf('Manufactured case             = %s\n', problem.case_name);
    fprintf('Anchored u0 homotopy finished with N = %d\n', problem.N);
    fprintf('Anchor construction           = %s\n', problem.anchor_label);
    fprintf('Solve success                 = %s\n', logical_text(result.success));
    fprintf('Reached target t = %.2e       = %s\n', problem.t_end, logical_text(result.reached_t_end));
    fprintf('Anchor RHS min f on quadrature = %.3e\n', result.anchor_rhs_min_f);
    fprintf('Anchor RHS clipped to zero    = %s\n', logical_text(result.anchor_rhs_clipped_to_zero));
    fprintf('Final continuation parameter t = %.2e\n', result.t_values(end));
    fprintf('Final residual inf-norm        = %.3e\n', result.residual_history(end));
    if problem.has_exact
        fprintf('Final grid L2 indicator        = %.3e\n', result.l2_error_indicator(end));
        fprintf('Final FD-H1 indicator          = %.3e\n', result.h1_fd_error_indicator(end));
        fprintf('Final FD-H2 indicator          = %.3e\n', result.h2_fd_error_indicator(end));
        fprintf('Final continuous L2 norm       = %.3e\n', result.cont_l2_error(end));
        fprintf('Final continuous H1 norm       = %.3e\n', result.cont_h1_error(end));
    else
        fprintf('Final error indicators         = unavailable (u_exact not supplied)\n');
    end
    if ~result.success
        fprintf('Failure reason                = %s\n', result.failure_reason);
    end
    fprintf('Final Jacobian sigma_min       = %.3e\n', result.jacobian_sigma_min(end));
    fprintf('Final Jacobian cond            = %.3e\n', result.jacobian_cond(end));
    fprintf('Final Hessian lambda_min min   = %.3e\n', result.hessian_lambda_min(end));
    fprintf('Final Hessian lambda_max max   = %.3e\n', result.hessian_lambda_max(end));
    fprintf('Final Hessian PSD ratio        = %.3f\n', result.hessian_psd_ratio(end));
    fprintf('Final Hessian indefinite ratio = %.3f\n', result.hessian_indef_ratio(end));
    fprintf('Figures saved under            = %s\n', output_dir);

    if problem.has_exact
        disp(table(result.t_values, result.l2_error_indicator, ...
            result.h1_fd_error_indicator, result.h2_fd_error_indicator, ...
            result.cont_l2_error, result.cont_h1_error, ...
            'VariableNames', {'t', 'GridL2Indicator', 'FDH1Indicator', 'FDH2Indicator', ...
            'ContL2Norm', 'ContH1Norm'}));
    end
end

function tf = has_n_sweep_request(options)
    tf = isstruct(options) && isfield(options, 'N_values') && numel(options.N_values) > 1;
end

function problem = apply_runtime_options(problem, options)
    if ~isstruct(options)
        return;
    end

    if isfield(options, 'N_values') && numel(options.N_values) == 1
        problem.N = max(2, round(options.N_values));
    end
    if isfield(options, 'N')
        problem.N = max(2, round(options.N));
    end
    if isfield(options, 'Q')
        problem.Q = max(4, round(options.Q));
    end
    if isfield(options, 'eval_points')
        problem.eval_points = max(21, round(options.eval_points));
    end
    if isfield(options, 'error_quad_points')
        problem.error_quad_points = max(12, round(options.error_quad_points));
    end
    if isfield(options, 'show_figures')
        problem.show_figures = logical(options.show_figures);
    end
    if isfield(options, 'acceptance_policy') && ~isempty(options.acceptance_policy)
        policy = lower(char(string(options.acceptance_policy)));
        switch policy
            case 'strict_convexity'
                problem.acceptance_policy = 'strict_convexity';
                problem.convexity_mode = 'strict';
            case 'diagnostic_only'
                problem.acceptance_policy = 'diagnostic_only';
                problem.convexity_mode = 'none';
                problem.line_search_allow_fallback = false;
            otherwise
                error('RuntimeOptions:UnknownAcceptancePolicy', ...
                    'Unknown acceptance_policy: %s', policy);
        end
    end
    if isfield(options, 'convexity_mode') && ~isempty(options.convexity_mode)
        problem.convexity_mode = char(string(options.convexity_mode));
    end
    if isfield(options, 'line_search_allow_fallback')
        problem.line_search_allow_fallback = logical(options.line_search_allow_fallback);
    end

    [external_anchor_fun, external_anchor_label, has_external_anchor] = ...
        resolve_external_anchor_options(options);
    if has_external_anchor
        problem.anchor_fun = external_anchor_fun;
        problem.anchor_label = external_anchor_label;
    end
end

function [anchor_fun, anchor_label, has_external_anchor] = resolve_external_anchor_options(options)
    anchor_fun = [];
    anchor_label = '';
    has_external_anchor = false;

    if isfield(options, 'u0_fun') && ~isempty(options.u0_fun)
        validate_external_anchor_fun(options.u0_fun, 'u0_fun');
        anchor_fun = options.u0_fun;
        has_external_anchor = true;
    end

    if isfield(options, 'anchor_fun') && ~isempty(options.anchor_fun)
        validate_external_anchor_fun(options.anchor_fun, 'anchor_fun');
        if ~has_external_anchor
            anchor_fun = options.anchor_fun;
            has_external_anchor = true;
        end
    end

    if has_external_anchor
        if isfield(options, 'anchor_label') && ~isempty(options.anchor_label)
            anchor_label = char(string(options.anchor_label));
        else
            anchor_label = 'external u0 anchor';
        end
    end
end

function validate_external_anchor_fun(candidate_fun, option_name)
    if ~isa(candidate_fun, 'function_handle')
        error('ExternalAnchor:InvalidFunctionHandle', ...
            '%s must be a function handle when supplied.', option_name);
    end
end

function output_dir = single_run_output_dir(case_name, options)
    folder_name = ['u0_anchor_monge_ampere_homotopy_' case_name];
    if isstruct(options) && isfield(options, 'run_tag') && ~isempty(options.run_tag)
        folder_name = [folder_name '_' char(string(options.run_tag))];
    end
    output_dir = fullfile(pwd, 'output', folder_name);
end

function framework_result = run_multiple_level_framework(base_problem, options)
    [N_values, target_p_values, target_t_values] = parse_multiple_level_targets(options);
    output_dir = multiple_level_output_dir(base_problem.case_name, options);
    per_level_dir = fullfile(output_dir, 'per_level');
    if ~exist(per_level_dir, 'dir')
        mkdir(per_level_dir);
    end

    fprintf('Running multiple-level framework for case %s\n', base_problem.case_name);
    fprintf('Levels N = [%s]\n', sprintf('%d ', N_values));
    if any(~isnan(target_p_values))
        fprintf('Using paired p = [%s]\n', sprintf('%d ', target_p_values));
    end

    num_levels = numel(N_values);
    rows = repmat(empty_multiple_level_row(), num_levels, 1);
    level_results = cell(num_levels, 1);

    prev_solution = [];
    prev_N = NaN;
    prev_t = base_problem.t_start;

    for k = 1:num_levels
        problem = apply_runtime_options(base_problem, options);
        problem.N = N_values(k);
        [problem.Q, problem.eval_points] = sweep_resolution_settings(base_problem, problem.N);
        problem.show_figures = false;
        problem = set_multiple_level_t_interval(problem, prev_t, target_t_values(k));

        initial_state = struct();
        if ~isempty(prev_solution)
            initial_state.U_init = prolong_homogeneous_coeffs(prev_solution, prev_N, problem.N);
        end

        run_output_dir = fullfile(per_level_dir, sprintf('level_%02d_N_%03d', k, problem.N));
        fprintf('\n===== MULTI-LEVEL RUN: level = %d, N = %d, start t = %.0e, target t = %.0e =====\n', ...
            k, problem.N, problem.t_start, problem.t_end);
        cpu_start = tic;
        result = run_single_case(problem, run_output_dir, initial_state);
        cpu_seconds = toc(cpu_start);

        rows(k) = collect_multiple_level_row(k, problem, result, target_p_values(k), cpu_seconds);
        level_results{k} = result;
        prev_solution = result.U_final;
        prev_N = problem.N;
        prev_t = problem.t_end;

        if ~result.success
            warning('MultipleLevel:StoppedEarly', ...
                'Stopping multiple-level framework at level %d because the level did not converge.', k);
            rows = rows(1:k);
            level_results = level_results(1:k);
            break;
        end
    end

    summary_table = struct2table(rows);
    writetable(summary_table, fullfile(output_dir, 'u0_anchor_multiple_level_summary.csv'));
    write_multiple_level_markdown(summary_table, base_problem, fullfile(output_dir, 'u0_anchor_multiple_level_summary.md'));
    plot_multiple_level_result(base_problem, summary_table, output_dir);

    disp(summary_table);
    fprintf('Multiple-level summary written under %s\n', output_dir);

    framework_result = struct();
    framework_result.case_name = base_problem.case_name;
    framework_result.output_dir = output_dir;
    framework_result.level_summary = summary_table;
    framework_result.level_results = level_results;
end

function [N_values, target_p_values, target_t_values] = parse_multiple_level_targets(options)
    if ~(isstruct(options) && isfield(options, 'N_values') && numel(options.N_values) > 1)
        error('Multiple-level mode requires N_values with at least two entries.');
    end

    N_values = round(options.N_values(:).');
    if any(~isfinite(N_values)) || any(N_values < 2)
        error('N_values must be finite integers >= 2.');
    end

    if any(diff(N_values) <= 0)
        error('Multiple-level N_values must be strictly increasing.');
    end

    if isfield(options, 'p_values') && ~isempty(options.p_values)
        target_p_values = round(options.p_values(:).');
        if numel(target_p_values) ~= numel(N_values)
            error('p_values must have the same length as N_values.');
        end
        target_t_values = 10 .^ (-target_p_values);
    elseif isfield(options, 'target_t_values') && ~isempty(options.target_t_values)
        target_t_values = double(options.target_t_values(:).');
        if numel(target_t_values) ~= numel(N_values)
            error('target_t_values must have the same length as N_values.');
        end
        target_p_values = NaN(size(target_t_values));
    else
        error('Multiple-level mode requires either p_values or target_t_values.');
    end

    if any(~isfinite(target_t_values)) || any(target_t_values <= 0) || any(target_t_values > 1)
        error('Every target t must satisfy 0 < t <= 1.');
    end
    if any(diff(target_t_values) >= 0)
        error('Multiple-level target t values must be strictly decreasing.');
    end
end

function output_dir = multiple_level_output_dir(case_name, options)
    folder_name = ['u0_anchor_monge_ampere_multiple_level_' case_name];
    if isstruct(options) && isfield(options, 'run_tag') && ~isempty(options.run_tag)
        folder_name = [folder_name '_' char(string(options.run_tag))];
    end
    output_dir = fullfile(pwd, 'output', folder_name);
end

function problem = set_multiple_level_t_interval(problem, start_t, target_t)
    if ~(target_t > 0 && target_t <= start_t && start_t <= 1)
        error('Expected 0 < target_t <= start_t <= 1 in the multiple-level framework.');
    end

    start_tol = t_value_tolerance(start_t);
    target_tol = t_value_tolerance(target_t);
    t_values = unique([problem.t_values(:); start_t; target_t]);
    t_values = sort(t_values, 'descend');
    keep_mask = (t_values <= start_t + start_tol) ...
        & (t_values >= target_t - target_tol);
    t_values = t_values(keep_mask);
    if isempty(t_values) || abs(t_values(1) - start_t) > start_tol
        t_values = [start_t; t_values(:)];
    end
    if abs(t_values(end) - target_t) > target_tol
        t_values = [t_values(:); target_t];
    end

    t_values = unique(t_values, 'stable');
    problem.t_values = t_values(:).';
    problem.t_start = problem.t_values(1);
    problem.t_end = problem.t_values(end);
end

function U_new = prolong_homogeneous_coeffs(U_prev, N_prev, N_new)
    U_prev_hat = reshape(U_prev, N_prev, N_prev);
    U_new_hat = zeros(N_new, N_new);
    copy_n = min(N_prev, N_new);
    U_new_hat(1:copy_n, 1:copy_n) = U_prev_hat(1:copy_n, 1:copy_n);
    U_new = U_new_hat(:);
end

function row = empty_multiple_level_row()
    row = struct( ...
        'Level', NaN, ...
        'N', NaN, ...
        'TargetP', NaN, ...
        'StartT', NaN, ...
        'TargetT', NaN, ...
        'Q', NaN, ...
        'EvalPoints', NaN, ...
        'Success', false, ...
        'ReachedTEnd', false, ...
        'FinalT', NaN, ...
        'FinalResidualInf', NaN, ...
        'FinalGridL2', NaN, ...
        'FinalFDH1', NaN, ...
        'FinalFDH2', NaN, ...
        'FinalContL2', NaN, ...
        'FinalContH1', NaN, ...
        'CPUSeconds', NaN, ...
        'FinalJacobianCond', NaN, ...
        'FailureReason', "" ...
    );
end

function row = collect_multiple_level_row(level_idx, problem, result, target_p, cpu_seconds)
    row = empty_multiple_level_row();
    row.Level = level_idx;
    row.N = problem.N;
    row.TargetP = target_p;
    row.StartT = problem.t_start;
    row.TargetT = problem.t_end;
    row.Q = problem.Q;
    row.EvalPoints = problem.eval_points;
    row.Success = result.success;
    row.ReachedTEnd = result.reached_t_end;
    row.FinalT = result.t_values(end);
    row.FinalResidualInf = result.residual_history(end);
    row.FinalGridL2 = final_indicator_or_nan(result, 'l2_error_indicator');
    row.FinalFDH1 = final_indicator_or_nan(result, 'h1_fd_error_indicator');
    row.FinalFDH2 = final_indicator_or_nan(result, 'h2_fd_error_indicator');
    row.FinalContL2 = final_indicator_or_nan(result, 'cont_l2_error');
    row.FinalContH1 = final_indicator_or_nan(result, 'cont_h1_error');
    row.CPUSeconds = cpu_seconds;
    row.FinalJacobianCond = result.jacobian_cond(end);
    row.FailureReason = string(result.failure_reason);
end

function write_multiple_level_markdown(summary_table, problem, output_path)
    fid = fopen(output_path, 'w');
    if fid < 0
        error('Could not open %s for writing.', output_path);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '# u0-anchor multiple-level summary\n\n');
    fprintf(fid, 'Case: `%s`\n\n', problem.case_name);
    fprintf(fid, '| Level | N | p | start t | target t | Success | Final t | Residual inf | ContL2 | ContH1 | GridL2 | FDH1 | cpu(s) |\n');
    fprintf(fid, '| ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n');
    for k = 1:height(summary_table)
        fprintf(fid, '| %d | %d | %s | %.0e | %.0e | %s | %.5g | %.3e | %s | %s | %s | %s | %.4f |\n', ...
            summary_table.Level(k), ...
            summary_table.N(k), ...
            optional_integer_text(summary_table.TargetP(k)), ...
            summary_table.StartT(k), ...
            summary_table.TargetT(k), ...
            logical_text(summary_table.Success(k)), ...
            summary_table.FinalT(k), ...
            summary_table.FinalResidualInf(k), ...
            optional_number_text(summary_table.FinalContL2(k)), ...
            optional_number_text(summary_table.FinalContH1(k)), ...
            optional_number_text(summary_table.FinalGridL2(k)), ...
            optional_number_text(summary_table.FinalFDH1(k)), ...
            summary_table.CPUSeconds(k));
    end
end

function plot_multiple_level_result(problem, summary_table, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fig = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [120, 120, 840, 560]);
    if problem.has_exact
        semilogy(summary_table.N, max(summary_table.FinalContL2, eps), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        hold on;
        semilogy(summary_table.N, max(summary_table.FinalContH1, eps), '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
        failed_mask = ~summary_table.Success;
        if any(failed_mask)
            semilogy(summary_table.N(failed_mask), max(summary_table.FinalContL2(failed_mask), eps), ...
                'rx', 'LineWidth', 1.8, 'MarkerSize', 10);
            legend('Cont L2', 'Cont H1', 'Failed levels', 'Location', 'southwest');
        else
            legend('Cont L2', 'Cont H1', 'Location', 'southwest');
        end
        ylabel('Final continuous norm');
        title(sprintf('u_0-anchor multiple-level summary (%s)', strrep(problem.case_name, '_', '\_')));
    else
        semilogy(summary_table.N, max(summary_table.FinalResidualInf, eps), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        ylabel('Final residual');
        title(sprintf('u_0-anchor multiple-level residual summary (%s)', strrep(problem.case_name, '_', '\_')));
    end
    grid on;
    xlabel('N');
    hide_axes_toolbars(fig);
    exportgraphics(fig, fullfile(output_dir, 'u0_anchor_multiple_level_error_vs_N.png'), 'Resolution', 180);
    drawnow;
end

function sweep_result = run_n_convergence_sweep(base_problem, options)
    N_values = unique(round(options.N_values(:).'));
    N_values = N_values(isfinite(N_values) & N_values >= 2);
    if isempty(N_values)
        error('N_values must contain at least one integer >= 2.');
    end
    [p_values, fixed_target_t, use_paired_targets, use_fixed_target_t] = ...
        parse_n_sweep_targets(options, N_values);

    output_dir = n_sweep_output_dir(base_problem.case_name, options);
    per_run_dir = fullfile(output_dir, 'per_N');
    if ~exist(per_run_dir, 'dir')
        mkdir(per_run_dir);
    end

    num_runs = numel(N_values);
    rows = repmat(empty_n_sweep_row(), num_runs, 1);

    fprintf('Running N-sweep for case %s with N = [%s]\n', ...
        base_problem.case_name, sprintf('%d ', N_values));
    if use_paired_targets
        fprintf('Using paired p = [%s], i.e. target t = 10^{-p} for each N.\n', ...
            sprintf('%d ', p_values));
    elseif use_fixed_target_t
        fprintf('Using fixed target t = %.0e for every N.\n', fixed_target_t);
    end
    for k = 1:num_runs
        N = N_values(k);
        problem = apply_runtime_options(base_problem, options);
        problem.N = N;
        [problem.Q, problem.eval_points] = sweep_resolution_settings(base_problem, N);
        problem.show_figures = false;
        target_p = NaN;
        target_t = problem.t_end;
        if use_paired_targets
            target_p = p_values(k);
            target_t = 10^(-target_p);
            problem = set_n_sweep_target_t(problem, target_t);
        elseif use_fixed_target_t
            target_t = fixed_target_t;
            problem = set_n_sweep_target_t(problem, target_t);
        end

        run_output_dir = fullfile(per_run_dir, sprintf('N_%03d', N));
        fprintf('\n===== N-SWEEP RUN: case = %s, N = %d, Q = %d, target t = %.0e =====\n', ...
            problem.case_name, problem.N, problem.Q, target_t);
        cpu_start = tic;
        result = run_single_case(problem, run_output_dir);
        cpu_seconds = toc(cpu_start);
        rows(k) = collect_n_sweep_row(problem, result, target_p, target_t, cpu_seconds);
    end

    summary_table = struct2table(rows);
    writetable(summary_table, fullfile(output_dir, 'u0_anchor_N_convergence.csv'));
    write_n_sweep_markdown(summary_table, base_problem, fullfile(output_dir, 'u0_anchor_N_convergence.md'));
    plot_n_sweep_result(base_problem, summary_table, output_dir);

    disp(summary_table);
    fprintf('N-sweep summary written under %s\n', output_dir);

    sweep_result = struct();
    sweep_result.case_name = base_problem.case_name;
    sweep_result.output_dir = output_dir;
    sweep_result.summary_table = summary_table;
end

function output_dir = n_sweep_output_dir(case_name, options)
    folder_name = ['u0_anchor_monge_ampere_homotopy_' case_name '_N_sweep'];
    if isstruct(options) && isfield(options, 'run_tag') && ~isempty(options.run_tag)
        folder_name = [folder_name '_' char(string(options.run_tag))];
    end
    output_dir = fullfile(pwd, 'output', folder_name);
end

function [Q, eval_points] = sweep_resolution_settings(base_problem, N)
    Q = max(24, 2 * N + 2);
    eval_points = max(base_problem.eval_points, 6 * N + 1);
end

function [p_values, fixed_target_t, use_paired_targets, use_fixed_target_t] = ...
        parse_n_sweep_targets(options, N_values)
    if isstruct(options) && isfield(options, 'p_values') && ~isempty(options.p_values)
        p_values = round(options.p_values(:).');
        if numel(p_values) ~= numel(N_values)
            error('p_values must have the same length as N_values.');
        end
        fixed_target_t = NaN;
        use_paired_targets = true;
        use_fixed_target_t = false;
        return;
    end

    if isstruct(options) && isfield(options, 'target_t') && ~isempty(options.target_t)
        fixed_target_t = double(options.target_t);
        if ~(isscalar(fixed_target_t) && isfinite(fixed_target_t) && fixed_target_t > 0)
            error('target_t must be a positive finite scalar.');
        end
        p_values = nan(size(N_values));
        use_paired_targets = false;
        use_fixed_target_t = true;
        return;
    end

    p_values = nan(size(N_values));
    fixed_target_t = NaN;
    use_paired_targets = false;
    use_fixed_target_t = false;
end

function problem = set_n_sweep_target_t(problem, target_t)
    if ~strcmp(problem.continuation_mode, 'scheduled')
        error('Paired p-values are only supported for scheduled continuation runs.');
    end
    if ~(target_t > 0 && target_t <= problem.t_start)
        error('Target t must satisfy 0 < t <= t_start.');
    end

    t_values = unique([problem.t_values(:); target_t]);
    t_values = sort(t_values, 'descend');
    problem.t_values = t_values.';
    problem.t_start = problem.t_values(1);
    problem.t_end = target_t;
end

function row = empty_n_sweep_row()
    row = struct( ...
        'N', NaN, ...
        'TargetP', NaN, ...
        'TargetT', NaN, ...
        'Q', NaN, ...
        'EvalPoints', NaN, ...
        'Success', false, ...
        'ReachedTEnd', false, ...
        'FinalT', NaN, ...
        'FinalResidualInf', NaN, ...
        'FinalGridL2', NaN, ...
        'FinalFDH1', NaN, ...
        'FinalFDH2', NaN, ...
        'FinalContL2', NaN, ...
        'FinalContH1', NaN, ...
        'CPUSeconds', NaN, ...
        'FinalJacobianCond', NaN, ...
        'FinalPSDRatio', NaN, ...
        'FinalIndefRatio', NaN, ...
        'FailureReason', "" ...
    );
end

function row = collect_n_sweep_row(problem, result, target_p, target_t, cpu_seconds)
    row = empty_n_sweep_row();
    row.N = problem.N;
    row.TargetP = target_p;
    row.TargetT = target_t;
    row.Q = problem.Q;
    row.EvalPoints = problem.eval_points;
    row.Success = result.success;
    row.ReachedTEnd = result.reached_t_end;
    row.FinalT = result.t_values(end);
    row.FinalResidualInf = result.residual_history(end);
    row.FinalGridL2 = final_indicator_or_nan(result, 'l2_error_indicator');
    row.FinalFDH1 = final_indicator_or_nan(result, 'h1_fd_error_indicator');
    row.FinalFDH2 = final_indicator_or_nan(result, 'h2_fd_error_indicator');
    row.FinalContL2 = final_indicator_or_nan(result, 'cont_l2_error');
    row.FinalContH1 = final_indicator_or_nan(result, 'cont_h1_error');
    row.CPUSeconds = cpu_seconds;
    row.FinalJacobianCond = result.jacobian_cond(end);
    row.FinalPSDRatio = result.hessian_psd_ratio(end);
    row.FinalIndefRatio = result.hessian_indef_ratio(end);
    row.FailureReason = string(result.failure_reason);
end

function value = final_indicator_or_nan(result, field_name)
    if result.has_exact && isfield(result, field_name) && ~isempty(result.(field_name))
        value = result.(field_name)(end);
    else
        value = NaN;
    end
end

function plot_n_sweep_result(problem, summary_table, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fig = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [120, 120, 840, 560]);
    N_values = summary_table.N;
    if problem.has_exact
        semilogy(N_values, max(summary_table.FinalContL2, eps), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        hold on;
        semilogy(N_values, max(summary_table.FinalContH1, eps), '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
        semilogy(nan, nan, '-d', 'LineWidth', 1.8, 'MarkerSize', 7);
        failed_mask = ~summary_table.Success;
        if any(failed_mask)
            semilogy(N_values(failed_mask), max(summary_table.FinalContL2(failed_mask), eps), ...
                'rx', 'LineWidth', 1.8, 'MarkerSize', 10);
            legend('网格 L2', '差分 H1', '差分 H2', '失败点', 'Location', 'southwest');
        else
            legend('网格 L2', '差分 H1', '差分 H2', 'Location', 'southwest');
        end
        if any(failed_mask)
            legend('Cont L2', 'Cont H1', 'Failed runs', 'Location', 'southwest');
        else
            legend('Cont L2', 'Cont H1', 'Location', 'southwest');
        end
        annotate_n_sweep_p_values(summary_table, summary_table.FinalContL2);
        ylabel('Final continuous norm');
        ylabel('最终误差指标');
        ylabel('Final continuous norm');
        title(n_sweep_plot_title(problem, summary_table));
    else
        semilogy(N_values, max(summary_table.FinalResidualInf, eps), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        ylabel('最终残差');
        title(n_sweep_plot_title(problem, summary_table));
    end
    grid on;
    xlabel('N');
    hide_axes_toolbars(fig);
    png_path = fullfile(output_dir, 'u0_anchor_error_vs_N.png');
    exportgraphics(fig, png_path, 'Resolution', 180);
    drawnow;
end

function write_n_sweep_markdown(summary_table, problem, output_path)
    fid = fopen(output_path, 'w');
    if fid < 0
        error('Could not open %s for writing.', output_path);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '# u0-anchor N convergence summary\n\n');
    fprintf(fid, 'Case: `%s`\n\n', problem.case_name);
    fprintf(fid, '| N | p | target t | Q | Success | Final t | Residual inf | ContL2 | ContH1 | GridL2 | FDH1 | cpu(s) | PSD ratio |\n');
    fprintf(fid, '| ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n');
    for k = 1:height(summary_table)
        fprintf(fid, '| %d | %s | %.0e | %d | %s | %.5g | %.3e | %s | %s | %s | %s | %.4f | %.3f |\n', ...
            summary_table.N(k), ...
            optional_integer_text(summary_table.TargetP(k)), ...
            summary_table.TargetT(k), ...
            summary_table.Q(k), ...
            logical_text(summary_table.Success(k)), ...
            summary_table.FinalT(k), ...
            summary_table.FinalResidualInf(k), ...
            optional_number_text(summary_table.FinalContL2(k)), ...
            optional_number_text(summary_table.FinalContH1(k)), ...
            optional_number_text(summary_table.FinalGridL2(k)), ...
            optional_number_text(summary_table.FinalFDH1(k)), ...
            summary_table.CPUSeconds(k), ...
            summary_table.FinalPSDRatio(k));
    end
end

function text = optional_number_text(value)
    if isnan(value)
        text = 'n/a';
    else
        text = sprintf('%.3e', value);
    end
end

function text = optional_integer_text(value)
    if isnan(value)
        text = 'n/a';
    else
        text = sprintf('%d', round(value));
    end
end

function title_text = n_sweep_plot_title(problem, summary_table)
    base_name = strrep(problem.case_name, '_', '\_');
    if any(~isnan(summary_table.TargetP))
        title_text = sprintf('u_0 锚定同伦误差关于 N 的变化（配对 p，%s）', base_name);
    elseif problem.has_exact
        title_text = sprintf('u_0 锚定同伦最终误差指标关于 N 的变化（%s）', base_name);
    else
        title_text = sprintf('u_0 锚定同伦最终残差关于 N 的变化（%s）', base_name);
    end
end

function annotate_n_sweep_p_values(summary_table, y_values)
    if ~any(~isnan(summary_table.TargetP))
        return;
    end

    text_y = max(y_values, eps) * 1.25;
    for k = 1:height(summary_table)
        if isnan(summary_table.TargetP(k))
            continue;
        end
        text(summary_table.N(k), text_y(k), sprintf('p=%d', round(summary_table.TargetP(k))), ...
            'FontSize', 9, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

function problem = default_problem(case_name)
    problem.N = 15;
    problem.Q = 32;
    problem.eval_points = 121;
    problem.error_quad_points = 96;
    problem.t_values = [1.0, 10 .^ (-(1:2:13)), 0.0];
    problem.t_start = problem.t_values(1);
    problem.t_end = problem.t_values(end);
    problem.continuation_mode = 'scheduled';
    problem.fixed_step = 0.10;
    problem.initial_step = 0.10;
    problem.min_step = 0.0025;
    problem.max_step = 0.10;
    problem.step_growth = 1.20;
    problem.accept_tol = 1e-8;
    problem.max_continuation_steps = 120;
    problem.newton_tol = 1e-14;
    problem.newton_maxit = 20;
    problem.line_search_min = 2^-12;
    problem.show_figures = true;
    problem.convexity_psd_drop_tol = 0.01;
    problem.convexity_indef_rise_tol = 0.01;
    problem.convexity_lam_drop_tol = 2.0;
    problem.convexity_score_tol = 1e-3;
    problem.convexity_target_psd_ratio = 0.995;
    problem.convexity_target_lam_min = 1e-8;
    problem.acceptance_policy = 'strict_convexity';
    problem.convexity_mode = 'strict';
    problem.relaxed_convexity_t_breaks = [0.995, 0.99, 0.97, 0.94];
    problem.relaxed_convexity_psd_floors = [0.72, 0.76, 0.82, 0.90];
    problem.relaxed_convexity_indef_ceils = [0.28, 0.24, 0.18, 0.10];
    problem.relaxed_convexity_lam_floors = [-0.20, -0.12, -0.05, -1e-3];
    problem.line_search_allow_fallback = false;
    problem.line_search_fallback_residual_factor = 0.98;
    problem.f_negative_tol = 1e-12;
    problem.anchor_fun = [];
    problem.diagnostics = struct();
    problem.has_exact = false;
    problem = configure_case(problem, case_name);
    if isempty(problem.anchor_fun)
        problem.anchor_label = '\Delta u_0 = 2\sqrt{f}';
    end
end

function problem = configure_case(problem, case_name)
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
        case {'mild_sextic_convex', 'separable_sextic_convex'}
            problem.case_name = 'mild_sextic_convex';
            a = 0.02;
            u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + a * (x.^6 + y.^6);
            u_xx = @(x, y) 1 + 30 * a * x.^4;
            u_yy = @(x, y) 1 + 30 * a * y.^4;
            problem.f = @(x, y) u_xx(x, y) .* u_yy(x, y);
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case {'scaled_quartic_convex', 'quartic_convex_scaled'}
            problem.case_name = 'scaled_quartic_convex';
            scale = 0.20;
            u_exact = @(x, y) scale * (0.5 * (x.^2 + y.^2) + x.^4 + y.^4);
            problem.f = @(x, y) (scale^2) * (1 + 12 * x.^2) .* (1 + 12 * y.^2);
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
        case {'trig_convex', 'quadratic_trig_convex'}
            problem.case_name = 'trig_convex';
            problem.N = 24;
            problem.Q = 48;
            problem.eval_points = 161;
            problem.continuation_mode = 'fixed';
            problem.fixed_step = 0.01;
            problem.initial_step = 0.01;
            problem.min_step = 1e-4;
            problem.max_step = 0.02;
            problem.step_growth = 1.05;
            problem.max_continuation_steps = 1000;
            problem.newton_maxit = 80;
            problem.line_search_min = 2^-20;
            problem.convexity_psd_drop_tol = 0.05;
            problem.convexity_indef_rise_tol = 0.05;
            problem.convexity_lam_drop_tol = 5.0;
            problem.convexity_score_tol = 0.0;
            problem.convexity_mode = 'none';
            problem.line_search_allow_fallback = true;
            problem.line_search_fallback_residual_factor = 1.02;
            alpha = 0.05;
            u_exact = @(x, y) 0.5 * (x.^2 + y.^2) + ...
                alpha * (sin(pi * x).^2 + sin(pi * y).^2);
            hx = @(x) 1 + 2 * alpha * pi^2 * cos(2 * pi * x);
            hy = @(y) 1 + 2 * alpha * pi^2 * cos(2 * pi * y);
            problem.f = @(x, y) hx(x) .* hy(y);
            problem.g = u_exact;
            problem = attach_exact_diagnostics(problem, u_exact);
        case 'quartic_convex'
            problem.case_name = 'quartic_convex';
            problem.f = @(x, y) (1 + 12 * x.^2) .* (1 + 12 * y.^2);
            problem.g = @(x, y) 0.5 * (x.^2 + y.^2) + x.^4 + y.^4;
            problem.anchor_fun = problem.g;
            problem.anchor_label = '分离型 quartic anchor';
            problem = attach_exact_diagnostics(problem, problem.g);
        case {'signed_branch_polynomial', 'signed_polynomial_branch'}
            problem.case_name = 'signed_branch_polynomial';
            amplitude = 0.05;
            branch_u = @(x, y) amplitude * x .* (1 - x) .* y .* (1 - y);
            u_xx = @(x, y) -2 * amplitude * y .* (1 - y);
            u_yy = @(x, y) -2 * amplitude * x .* (1 - x);
            u_xy = @(x, y) amplitude * (1 - 2 * x) .* (1 - 2 * y);
            problem.f = @(x, y) u_xx(x, y) .* u_yy(x, y) - u_xy(x, y).^2;
            problem.g = @(x, y) zeros(size(x));
            problem = attach_exact_diagnostics(problem, branch_u);
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

function result = solve_t_homotopy(problem, initial_state)
    if nargin < 2
        initial_state = struct();
    end
    ops = build_b_matrix_operators(problem.N, problem.Q, problem.eval_points, problem.error_quad_points);
    data = build_example_data(problem, ops);

    if isstruct(initial_state) && isfield(initial_state, 'U_init') && ~isempty(initial_state.U_init)
        U_vec = initial_state.U_init(:);
        if numel(U_vec) ~= problem.N^2
            error('Initial coefficient vector length must equal N^2 at each multiple-level stage.');
        end
    else
        U_vec = data.U0_vec;
    end
    if strcmp(problem.continuation_mode, 'scheduled')
        num_slots = numel(problem.t_values);
    else
        num_slots = problem.max_continuation_steps;
    end

    t_history = zeros(num_slots, 1);
    residual_history = zeros(num_slots, 1);
    newton_steps = zeros(num_slots, 1);
    l2_error_indicator = zeros(num_slots, 1);
    h1_fd_error_indicator = zeros(num_slots, 1);
    h2_fd_error_indicator = zeros(num_slots, 1);
    cont_l2_error = zeros(num_slots, 1);
    cont_h1_error = zeros(num_slots, 1);
    jacobian_sigma_min = zeros(num_slots, 1);
    jacobian_cond = zeros(num_slots, 1);
    hessian_lambda_min = zeros(num_slots, 1);
    hessian_lambda_max = zeros(num_slots, 1);
    hessian_psd_ratio = zeros(num_slots, 1);
    hessian_nsd_ratio = zeros(num_slots, 1);
    hessian_indef_ratio = zeros(num_slots, 1);

    current_t = problem.t_start;
    step = problem.initial_step;
    [U_vec, iter_count, res_norm, level_status] = solve_level(U_vec, current_t, ops, data, problem);
    failure_reason = '';

    accepted = 1;
    t_history(accepted) = current_t;
    residual_history(accepted) = res_norm;
    newton_steps(accepted) = iter_count;
    [l2_error_indicator(accepted), h1_fd_error_indicator(accepted), h2_fd_error_indicator(accepted)] = ...
        compute_error_indicators(U_vec, ops, data, problem);
    [cont_l2_error(accepted), cont_h1_error(accepted)] = ...
        compute_continuous_error_norms(U_vec, ops, data, problem);
    [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, data));
    [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, data);
    jacobian_sigma_min(accepted) = sigma_min;
    jacobian_cond(accepted) = cond_val;
    hessian_lambda_min(accepted) = lam_min;
    hessian_lambda_max(accepted) = lam_max;
    hessian_psd_ratio(accepted) = psd_ratio;
    hessian_nsd_ratio(accepted) = nsd_ratio;
    hessian_indef_ratio(accepted) = indef_ratio;
    print_continuation_status(current_t, iter_count, res_norm, [], sigma_min, cond_val, ...
        lam_min, lam_max, psd_ratio, indef_ratio, l2_error_indicator(accepted), problem.has_exact);

    if ~level_status.converged
        failure_reason = sprintf('Initial level at t = %.2e did not converge (%s).', current_t, level_status.reason);
    end

    while current_t > problem.t_end && accepted < num_slots && isempty(failure_reason)
        if strcmp(problem.continuation_mode, 'scheduled')
            target_t = problem.t_values(accepted + 1);
            current_step = current_t - target_t;
        elseif strcmp(problem.continuation_mode, 'fixed')
            current_step = problem.fixed_step;
            target_t = max(problem.t_end, current_t - current_step);
        else
            current_step = step;
            target_t = max(problem.t_end, current_t - current_step);
        end
        if abs(target_t - problem.t_end) <= t_value_tolerance([target_t, problem.t_end, current_t])
            target_t = problem.t_end;
        end

        [U_try, iter_count, res_norm, level_status] = solve_level(U_vec, target_t, ops, data, problem);
        level_accepted = level_status.converged && res_norm <= problem.accept_tol;

        if level_accepted
            U_vec = U_try;
            current_t = target_t;
            accepted = accepted + 1;
            t_history(accepted) = current_t;
            residual_history(accepted) = res_norm;
            newton_steps(accepted) = iter_count;
            [l2_error_indicator(accepted), h1_fd_error_indicator(accepted), h2_fd_error_indicator(accepted)] = ...
                compute_error_indicators(U_vec, ops, data, problem);
            [cont_l2_error(accepted), cont_h1_error(accepted)] = ...
                compute_continuous_error_norms(U_vec, ops, data, problem);
            [sigma_min, cond_val] = jacobian_diagnostics(jacobian(U_vec, current_t, ops, data));
            [lam_min, lam_max, psd_ratio, nsd_ratio, indef_ratio] = hessian_branch_diagnostics(U_vec, ops, data);
            jacobian_sigma_min(accepted) = sigma_min;
            jacobian_cond(accepted) = cond_val;
            hessian_lambda_min(accepted) = lam_min;
            hessian_lambda_max(accepted) = lam_max;
            hessian_psd_ratio(accepted) = psd_ratio;
            hessian_nsd_ratio(accepted) = nsd_ratio;
            hessian_indef_ratio(accepted) = indef_ratio;
            print_continuation_status(current_t, iter_count, res_norm, current_step, sigma_min, cond_val, ...
                lam_min, lam_max, psd_ratio, indef_ratio, l2_error_indicator(accepted), problem.has_exact);

            if strcmp(problem.continuation_mode, 'adaptive') && current_t > problem.t_end
                step = min(problem.max_step, step * problem.step_growth);
            end
        else
            if strcmp(problem.continuation_mode, 'scheduled')
                failure_reason = sprintf(['Scheduled continuation failed before reaching t = 0; ' ...
                    'level t = %.2e did not converge (%s).'], target_t, level_status.reason);
                warning('ContinuationFailed:ScheduledStep', ...
                    'Scheduled continuation failed at t = %.4g with residual %.3e (%s).', ...
                    target_t, res_norm, level_status.reason);
                break;
            elseif strcmp(problem.continuation_mode, 'fixed')
                failure_reason = sprintf(['Fixed-step continuation failed before reaching t = 0; ' ...
                    'level t = %.2e did not converge (%s).'], target_t, level_status.reason);
                warning('ContinuationFailed:FixedStep', ...
                    'Fixed-step continuation failed at t = %.4g with residual %.3e (%s).', ...
                    target_t, res_norm, level_status.reason);
                break;
            else
                new_step = step / 2;
                fprintf('reject t = %.2e, residual = %.3e, reason = %s, step %.4f -> %.4f\n', ...
                    target_t, res_norm, level_status.reason, step, new_step);
                step = new_step;
                if step < problem.min_step
                failure_reason = sprintf(['Adaptive continuation stopped before reaching t = 0; ' ...
                    'last rejected level t = %.2e (%s).'], target_t, level_status.reason);
                warning('ContinuationFailed:MinStep', ...
                    'Adaptive continuation stopped before reaching t = 0; required step %.4g below min_step %.4g.', ...
                    step, problem.min_step);
                break;
                end
            end
        end
    end

    reached_t_end = abs(current_t - problem.t_end) <= ...
        t_value_tolerance([current_t, problem.t_end]);
    if reached_t_end
        success = true;
        failure_reason = '';
    else
        success = false;
        if isempty(failure_reason)
            if ~level_status.converged
                failure_reason = sprintf('Newton solve failed before reaching t = 0 (%s).', level_status.reason);
            elseif accepted >= problem.max_continuation_steps
                failure_reason = sprintf('Maximum continuation steps (%d) reached before t = 0.', ...
                    problem.max_continuation_steps);
            else
                failure_reason = 'Continuation terminated before reaching t = 0.';
            end
        end
    end

    result.t_values = t_history(1:accepted);
    result.residual_history = residual_history(1:accepted);
    result.newton_steps = newton_steps(1:accepted);
    result.l2_error_indicator = l2_error_indicator(1:accepted);
    result.h1_fd_error_indicator = h1_fd_error_indicator(1:accepted);
    result.h2_fd_error_indicator = h2_fd_error_indicator(1:accepted);
    result.cont_l2_error = cont_l2_error(1:accepted);
    result.cont_h1_error = cont_h1_error(1:accepted);
    result.l2_error = result.l2_error_indicator;
    result.h1_error = result.h1_fd_error_indicator;
    result.h2_error = result.h2_fd_error_indicator;
    result.jacobian_sigma_min = jacobian_sigma_min(1:accepted);
    result.jacobian_cond = jacobian_cond(1:accepted);
    result.hessian_lambda_min = hessian_lambda_min(1:accepted);
    result.hessian_lambda_max = hessian_lambda_max(1:accepted);
    result.hessian_psd_ratio = hessian_psd_ratio(1:accepted);
    result.hessian_nsd_ratio = hessian_nsd_ratio(1:accepted);
    result.hessian_indef_ratio = hessian_indef_ratio(1:accepted);
    result.x_eval = map_to_physical(ops.Xe);
    result.y_eval = map_to_physical(ops.Ye);
    result.has_exact = problem.has_exact;
    result.success = success;
    result.reached_t_end = reached_t_end;
    result.failure_reason = failure_reason;
    result.anchor_rhs_min_f = data.anchor_rhs_min_f;
    result.anchor_rhs_clipped_to_zero = data.anchor_rhs_clipped_to_zero;
    result.u_exact_eval = data.exact_eval;
    result.u0_eval = data.u0_eval;
    result.u_num_eval = total_solution_on_eval(U_vec, ops, data);
    result.U_final = U_vec;
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
        @(z) convexity_metrics(z, ops, data), ...
        t_val, ...
        problem, ...
        problem.newton_tol, ...
        problem.newton_maxit, ...
        problem.line_search_min);
end

function data = build_example_data(problem, ops)
    xq = map_to_physical(ops.Xq);
    yq = map_to_physical(ops.Yq);
    xe = map_to_physical(ops.Xe);
    ye = map_to_physical(ops.Ye);
    xc = map_to_physical(ops.Xc);
    yc = map_to_physical(ops.Yc);

    f_phys_q = problem.f(xq, yq);
    min_f_phys = min(f_phys_q(:));
    clipped_to_zero = false;
    f_q = f_phys_q / 16;
    Ub_q = boundary_extension_from_fun(xq, yq, problem.g);
    Ub_eval = boundary_extension_from_fun(xe, ye, problem.g);
    Ub_cont = boundary_extension_from_fun(xc, yc, problem.g);
    [Ub_cont_x, Ub_cont_y] = function_gradient_from_fun(@(x, y) boundary_extension_from_fun(x, y, problem.g), xc, yc);

    lift = differentiate_nodal_field(Ub_q(:), ops);
    if isa(problem.anchor_fun, 'function_handle')
        anchor_q = problem.anchor_fun(xq, yq);
        anchor_eval = problem.anchor_fun(xe, ye);
        anchor_hom_q = anchor_q - Ub_q;
        U0_hat = project_to_phi(anchor_hom_q, ops);
        U0_vec = U0_hat(:);
        u0_q = ops.T0 * U0_vec + Ub_q(:);
        u0_eval = ops.Phi_eval * U0_hat * ops.Phi_eval.' + Ub_eval;
    else
        if min_f_phys < -problem.f_negative_tol
            error('AnchorData:NegativeRHS', ...
                ['Anchor construction requires nonnegative f on the quadrature grid. ' ...
                 'Detected min(f) = %.3e below tolerance %.3e.'], ...
                min_f_phys, problem.f_negative_tol);
        end
        clipped_to_zero = any(f_phys_q(:) < 0);
        if clipped_to_zero
            warning('AnchorData:ClippedRHS', ...
                ['Small negative f values detected on the quadrature grid (min(f) = %.3e). ' ...
                 'Clipping them to zero before sqrt in the anchor construction.'], min_f_phys);
        end
        sqrt_f_phys = sqrt(max(f_phys_q(:), 0));
        poisson_source_vec = 0.5 * sqrt_f_phys ...
            - (lift.xx_vals + lift.yy_vals);
        U0_vec = solve_anchor_poisson(ops, poisson_source_vec);
        u0_q = ops.T0 * U0_vec + Ub_q(:);
        U0_hat = reshape(U0_vec, ops.N, ops.N);
        u0_eval = ops.Phi_eval * U0_hat * ops.Phi_eval.' + Ub_eval;
    end

    data.f_proj = ops.T1 * f_q(:);
    data.u0_q_vec = u0_q(:);
    data.Ub_q_vec = Ub_q(:);
    data.U0_vec = U0_vec;
    data.anchor_rhs_min_f = min_f_phys;
    data.anchor_rhs_clipped_to_zero = clipped_to_zero;
    data.lift_xx = lift.xx_vals;
    data.lift_yy = lift.yy_vals;
    data.lift_xy = lift.xy_vals;
    if problem.has_exact
        data.exact_eval = problem.diagnostics.u_exact(xe, ye);
        data.exact_cont = problem.diagnostics.u_exact(xc, yc);
        [data.exact_cont_x, data.exact_cont_y] = ...
            function_gradient_from_fun(problem.diagnostics.u_exact, xc, yc);
    else
        data.exact_eval = nan(size(xe));
        data.exact_cont = nan(size(xc));
        data.exact_cont_x = nan(size(xc));
        data.exact_cont_y = nan(size(xc));
    end
    data.u0_eval = u0_eval;
    data.Ub_eval = Ub_eval;
    data.Ub_cont = Ub_cont;
    data.Ub_cont_x = Ub_cont_x;
    data.Ub_cont_y = Ub_cont_y;
    data.cont_weights = ops.Wc;
    data.problem_g = problem.g;
end

function U0_vec = solve_anchor_poisson(ops, source_vec)
    U0_vec = (ops.L + 1e-12 * speye(size(ops.L))) \ (ops.T1 * source_vec);
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

function coeffs = project_to_phi(U_vals, ops)
    coeffs = ops.Pphi * U_vals * ops.Pphi.';
end

function out = differentiate_nodal_field(field_vec, ops)
    leg_vec = ops.forward_legendre * field_vec;
    out.xx_vals = ops.nodal_from_legendre * kron(ops.Id_full, ops.B5) * leg_vec;
    out.yy_vals = ops.nodal_from_legendre * kron(ops.B5, ops.Id_full) * leg_vec;
    out.xy_vals = ops.nodal_from_legendre * kron(ops.B6, ops.B6) * leg_vec;
end

function [fx, fy] = function_gradient_from_fun(ufun, x, y)
    h = 1e-30;
    fx = imag(ufun(x + 1i * h, y)) / h;
    fy = imag(ufun(x, y + 1i * h)) / h;
end

function r = residual(U_vec, t_val, ops, data)
    u = second_derivatives(U_vec, ops, data);
    u_vals = ops.T0 * U_vec + data.Ub_q_vec;
    mass_vec = ops.T1 * (u_vals - data.u0_q_vec);
    det_vec = ops.T1 * (u.xx .* u.yy - u.xy .^ 2);

    r = (t_val / 16) * mass_vec + (1 - t_val) * det_vec - (1 - t_val) * data.f_proj;
end

function J = jacobian(U_vec, t_val, ops, data)
    u = second_derivatives(U_vec, ops, data);
    nonlin_jac = ops.T1 * ( ...
        diag(u.yy) * ops.T2 + ...
        diag(u.xx) * ops.T3 - ...
        2 * diag(u.xy) * ops.T4);

    J = (t_val / 16) * ops.M + (1 - t_val) * nonlin_jac;
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

function ok = convexity_acceptable(trial_metrics, current_metrics, problem, t_val)
    if strcmp(problem.convexity_mode, 'none')
        % In research mode we still compute Hessian diagnostics, but they do not gate step acceptance.
        ok = true;
        return;
    end

    strongly_convex = ...
        trial_metrics.psd_ratio >= problem.convexity_target_psd_ratio && ...
        trial_metrics.lam_min >= problem.convexity_target_lam_min;

    nonworsening = ...
        trial_metrics.psd_ratio >= current_metrics.psd_ratio - problem.convexity_psd_drop_tol && ...
        trial_metrics.indef_ratio <= current_metrics.indef_ratio + problem.convexity_indef_rise_tol && ...
        trial_metrics.lam_min >= current_metrics.lam_min - problem.convexity_lam_drop_tol;

    improving = trial_metrics.score <= current_metrics.score - problem.convexity_score_tol;
    if strcmp(problem.convexity_mode, 'relaxed_trig')
        envelope_ok = convexity_envelope_acceptable(trial_metrics, problem, t_val);
        ok = strongly_convex || nonworsening || improving || envelope_ok;
    else
        ok = strongly_convex || nonworsening || improving;
    end
end

function ok = convexity_envelope_acceptable(trial_metrics, problem, t_val)
    [psd_floor, indef_ceil, lam_floor] = relaxed_convexity_envelope(problem, t_val);
    ok = trial_metrics.psd_ratio >= psd_floor && ...
        trial_metrics.indef_ratio <= indef_ceil && ...
        trial_metrics.lam_min >= lam_floor;
end

function [psd_floor, indef_ceil, lam_floor] = relaxed_convexity_envelope(problem, t_val)
    idx = find(t_val >= problem.relaxed_convexity_t_breaks, 1, 'first');
    if isempty(idx)
        psd_floor = problem.convexity_target_psd_ratio;
        indef_ceil = max(0, 1 - problem.convexity_target_psd_ratio);
        lam_floor = problem.convexity_target_lam_min;
        return;
    end

    psd_floor = problem.relaxed_convexity_psd_floors(idx);
    indef_ceil = problem.relaxed_convexity_indef_ceils(idx);
    lam_floor = problem.relaxed_convexity_lam_floors(idx);
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

function [l2e, h1e] = compute_continuous_error_norms(U_vec, ops, data, problem)
    if ~problem.has_exact
        l2e = nan;
        h1e = nan;
        return;
    end

    [u_num, ux_num, uy_num] = total_solution_on_continuous_grid(U_vec, ops, data);
    err = u_num - data.exact_cont;
    err_x = ux_num - data.exact_cont_x;
    err_y = uy_num - data.exact_cont_y;

    weights = data.cont_weights;
    l2e = sqrt(sum(weights(:) .* abs(err(:)) .^ 2));
    h1e = sqrt(sum(weights(:) .* (abs(err_x(:)) .^ 2 + abs(err_y(:)) .^ 2)));
end

function print_continuation_status(t_val, iter_count, res_norm, step, sigma_min, cond_val, ...
        lam_min, lam_max, psd_ratio, indef_ratio, l2_error, has_exact)
    if isempty(step)
        prefix = 't = %.2e, Newton steps = %2d, residual = %.3e, sigma_min = %.3e, cond = %.3e, ';
        values = {t_val, iter_count, res_norm, sigma_min, cond_val};
    else
        prefix = ['t = %.2e, Newton steps = %2d, residual = %.3e, step = %.4f, ' ...
            'sigma_min = %.3e, cond = %.3e, '];
        values = {t_val, iter_count, res_norm, step, sigma_min, cond_val};
    end

    if has_exact
        suffix = 'lam_min = %.3e, lam_max = %.3e, psd = %.3f, indef = %.3f, GridL2 = %.3e\n';
        values = [values, {lam_min, lam_max, psd_ratio, indef_ratio, l2_error}];
    else
        suffix = 'lam_min = %.3e, lam_max = %.3e, psd = %.3f, indef = %.3f\n';
        values = [values, {lam_min, lam_max, psd_ratio, indef_ratio}];
    end

    fprintf([prefix suffix], values{:});
end

function total_eval = total_solution_on_eval(U_vec, ops, data)
    U_hat = reshape(U_vec, ops.N, ops.N);
    U_eval = ops.Phi_eval * U_hat * ops.Phi_eval.';
    total_eval = U_eval + data.Ub_eval;
end

function [u_cont, ux_cont, uy_cont] = total_solution_on_continuous_grid(U_vec, ops, data)
    x_cont = map_to_physical(ops.Xc);
    y_cont = map_to_physical(ops.Yc);
    h = 1e-30;

    u_cont = real(evaluate_total_solution_tensor(U_vec, x_cont, y_cont, ops, data));
    ux_cont = imag(evaluate_total_solution_tensor(U_vec, x_cont + 1i * h, y_cont, ops, data)) / h;
    uy_cont = imag(evaluate_total_solution_tensor(U_vec, x_cont, y_cont + 1i * h, ops, data)) / h;
end

function total_vals = evaluate_total_solution_tensor(U_vec, x_phys, y_phys, ops, data)
    x_line = x_phys(1, :).';
    y_line = y_phys(:, 1);
    xi_line = 2 * x_line - 1;
    eta_line = 2 * y_line - 1;

    Phi_x = legendre_values(xi_line, ops.full_deg) * ops.B3;
    Phi_y = legendre_values(eta_line, ops.full_deg) * ops.B3;
    U_hat = reshape(U_vec, ops.N, ops.N);
    U_hom = Phi_y * U_hat * Phi_x.';

    total_vals = U_hom + boundary_extension_from_fun(x_phys, y_phys, data.problem_g);
end

function plot_result(problem, result, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fig1 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [100, 100, 1400, 820]);

    subplot(2, 2, 1);
    surf(result.x_eval, result.y_eval, result.u0_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u_0(x,y)');
    title('Poisson 锚定初值');
    colorbar;

    subplot(2, 2, 2);
    surf(result.x_eval, result.y_eval, result.u_num_eval, 'EdgeColor', 'none');
    view(45, 30);
    xlabel('x');
    ylabel('y');
    zlabel('u_{t,N}(x,y)');
    title(sprintf('数值解，t = %.0e', result.t_values(end)));
    colorbar;

    if problem.has_exact
        subplot(2, 2, 3);
        surf(result.x_eval, result.y_eval, result.u_exact_eval, 'EdgeColor', 'none');
        view(45, 30);
        xlabel('x');
        ylabel('y');
        zlabel('u(x,y)');
        title('精确解');
        colorbar;

        subplot(2, 2, 4);
        surf(result.x_eval, result.y_eval, result.error_eval, 'EdgeColor', 'none');
        view(45, 30);
        xlabel('x');
        ylabel('y');
        zlabel('error');
        title('点态误差');
        colorbar;
    else
        subplot(2, 2, 3);
        surf(result.x_eval, result.y_eval, result.u_num_eval - result.u0_eval, 'EdgeColor', 'none');
        view(45, 30);
        xlabel('x');
        ylabel('y');
        zlabel('u-u_0');
        title('同伦修正');
        colorbar;
    end

    hide_axes_toolbars(fig1);
    exportgraphics(fig1, fullfile(output_dir, 'u0_anchor_homotopy_surfaces.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig1);
    end

    fig2 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [120, 120, 760, 520]);
    [curve_x, tick_positions, tick_labels, display_mask] = homotopy_curve_x_values(result.t_values);
    if problem.has_exact
        semilogy(curve_x, result.l2_error(display_mask), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
        hold on;
        semilogy(curve_x, result.h1_error(display_mask), '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
        semilogy(curve_x, result.h2_error(display_mask), '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
        grid on;
        xlabel('-log_{10}(t)');
        ylabel('误差指标');
        title(sprintf('u_0 锚定同伦误差指标关于 -log_{10}(t) 的变化，N = %d', problem.N));
        legend('网格 L2', '差分 H1', '差分 H2', 'Location', 'southwest');
    else
        semilogy(curve_x, max(result.residual_history(display_mask), eps), '-o', ...
            'LineWidth', 1.6, 'MarkerSize', 7);
        grid on;
        xlabel('-log_{10}(t)');
        ylabel('残差');
        title(sprintf('u_0 锚定同伦残差关于 -log_{10}(t) 的变化，N = %d', problem.N));
    end
    xticks(tick_positions);
    xticklabels(tick_labels);
    hide_axes_toolbars(fig2);
    exportgraphics(fig2, fullfile(output_dir, 'u0_anchor_homotopy_error_curves.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig2);
    end

    fig_log = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [140, 140, 760, 520]);
    if problem.has_exact
        semilogy(result.t_values + 1e-14, result.l2_error, '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        hold on;
        semilogy(result.t_values + 1e-14, result.h1_error, '-s', 'LineWidth', 1.8, 'MarkerSize', 7);
        semilogy(result.t_values + 1e-14, result.h2_error, '-d', 'LineWidth', 1.8, 'MarkerSize', 7);
        set(gca, 'XDir', 'reverse');
        grid on;
        xlabel('t');
        ylabel('对数尺度误差指标');
        title('误差指标对 t 的对数图');
        legend('网格 L2', '差分 H1', '差分 H2', 'Location', 'southwest');
    else
        semilogy(result.t_values + 1e-14, max(result.residual_history, eps), '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
        set(gca, 'XDir', 'reverse');
        grid on;
        xlabel('t');
        ylabel('residual');
        title('残差对 t 的对数图');
    end
    hide_axes_toolbars(fig_log);
    exportgraphics(fig_log, fullfile(output_dir, 'u0_anchor_log_error_vs_t.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig_log);
    end

    fig3 = figure('Visible', figure_visibility(problem), 'Color', 'w', 'Position', [140, 140, 820, 720]);

    subplot(2, 1, 1);
    if problem.has_exact
        semilogy(result.t_values + 1e-14, result.l2_error, '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
        hold on;
        semilogy(result.t_values + 1e-14, result.h1_error, '-s', 'LineWidth', 1.6, 'MarkerSize', 7);
        semilogy(result.t_values + 1e-14, result.h2_error, '-d', 'LineWidth', 1.6, 'MarkerSize', 7);
        set(gca, 'XDir', 'reverse');
        grid on;
        xlabel('t');
        ylabel('误差指标');
        title('误差指标关于 t 的变化');
        legend('网格 L2', '差分 H1', '差分 H2', 'Location', 'southwest');
    else
        semilogy(result.t_values + 1e-14, max(result.residual_history, eps), '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
        set(gca, 'XDir', 'reverse');
        grid on;
        xlabel('t');
        ylabel('residual');
        title('残差关于 t 的变化');
    end

    subplot(2, 1, 2);
    yyaxis left;
    semilogy(result.t_values + 1e-14, result.residual_history, '-^', 'LineWidth', 1.6, 'MarkerSize', 7);
    ylabel('residual');
    yyaxis right;
    plot(result.t_values, result.newton_steps, '-o', 'LineWidth', 1.6, 'MarkerSize', 7);
    ylabel('Newton steps');
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('t');
    title('求解器诊断关于 t 的变化');

    hide_axes_toolbars(fig3);
    exportgraphics(fig3, fullfile(output_dir, 'u0_anchor_error_vs_t.png'), 'Resolution', 180);
    if ~problem.show_figures
        close(fig3);
    end
    drawnow;
end

function export_history(result, output_dir)
    if result.has_exact
        history_table = table( ...
            result.t_values, ...
            result.residual_history, ...
            result.newton_steps, ...
            result.l2_error_indicator, ...
            result.h1_fd_error_indicator, ...
            result.h2_fd_error_indicator, ...
            result.cont_l2_error, ...
            result.cont_h1_error, ...
            result.jacobian_sigma_min, ...
            result.jacobian_cond, ...
            result.hessian_lambda_min, ...
            result.hessian_lambda_max, ...
            result.hessian_psd_ratio, ...
            result.hessian_indef_ratio, ...
            'VariableNames', { ...
                't', ...
                'ResidualInf', ...
                'NewtonSteps', ...
                'GridL2Indicator', ...
                'FDH1Indicator', ...
                'FDH2Indicator', ...
                'ContL2Indicator', ...
                'ContH1Indicator', ...
                'JacobianSigmaMin', ...
                'JacobianCond', ...
                'HessianLambdaMin', ...
                'HessianLambdaMax', ...
                'HessianPSDRatio', ...
                'HessianIndefRatio'});
    else
        history_table = table( ...
            result.t_values, ...
            result.residual_history, ...
            result.newton_steps, ...
            result.jacobian_sigma_min, ...
            result.jacobian_cond, ...
            result.hessian_lambda_min, ...
            result.hessian_lambda_max, ...
            result.hessian_psd_ratio, ...
            result.hessian_indef_ratio, ...
            'VariableNames', { ...
                't', ...
                'ResidualInf', ...
                'NewtonSteps', ...
                'JacobianSigmaMin', ...
                'JacobianCond', ...
                'HessianLambdaMin', ...
                'HessianLambdaMax', ...
                'HessianPSDRatio', ...
                'HessianIndefRatio'});
    end

    writetable(history_table, fullfile(output_dir, 'u0_anchor_history.csv'));
end

function export_terminal_state(result, problem, output_dir)
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    anchor_label = char(string(problem.anchor_label));
    terminal_state = struct( ...
        'U_final', result.U_final, ...
        'u_num_eval', result.u_num_eval, ...
        'u0_eval', result.u0_eval, ...
        'x_eval', result.x_eval, ...
        'y_eval', result.y_eval, ...
        't_final', final_scalar_or_nan(result.t_values), ...
        'residual_final', final_scalar_or_nan(result.residual_history), ...
        'jacobian_cond_final', final_scalar_or_nan(result.jacobian_cond), ...
        'psd_ratio_final', final_scalar_or_nan(result.hessian_psd_ratio), ...
        'indef_ratio_final', final_scalar_or_nan(result.hessian_indef_ratio), ...
        'anchor_label', anchor_label, ...
        'acceptance_policy', char(string(problem.acceptance_policy)), ...
        'convexity_mode', char(string(problem.convexity_mode)), ...
        'line_search_allow_fallback', logical(problem.line_search_allow_fallback));
    save(fullfile(output_dir, 'terminal_state.mat'), '-struct', 'terminal_state');

    terminal_metrics = table( ...
        problem.N, ...
        problem.Q, ...
        problem.eval_points, ...
        result.success, ...
        result.reached_t_end, ...
        terminal_state.t_final, ...
        terminal_state.residual_final, ...
        terminal_state.jacobian_cond_final, ...
        terminal_state.psd_ratio_final, ...
        terminal_state.indef_ratio_final, ...
        string(problem.acceptance_policy), ...
        string(problem.convexity_mode), ...
        logical(problem.line_search_allow_fallback), ...
        'VariableNames', { ...
            'N', ...
            'Q', ...
            'EvalPoints', ...
            'Success', ...
            'ReachedTEnd', ...
            't_final', ...
            'residual_final', ...
            'jacobian_cond_final', ...
            'psd_ratio_final', ...
            'indef_ratio_final', ...
            'acceptance_policy', ...
            'convexity_mode', ...
            'line_search_allow_fallback'});
    writetable(terminal_metrics, fullfile(output_dir, 'terminal_metrics.csv'));

    fid = fopen(fullfile(output_dir, 'u0_definition.txt'), 'w');
    if fid < 0
        error('TerminalExport:WriteFailed', 'Unable to open u0_definition.txt for writing.');
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'u0 definition for case %s\n', problem.case_name);
    fprintf(fid, 'Anchor label: %s\n', anchor_label);
    fprintf(fid, 'Acceptance policy: %s\n', char(string(problem.acceptance_policy)));
    fprintf(fid, 'Convexity mode: %s\n', char(string(problem.convexity_mode)));
    fprintf(fid, 'Line-search fallback enabled: %s\n', logical_text(problem.line_search_allow_fallback));
    fprintf(fid, 'Terminal state artifacts exported under %s\n', output_dir);
    clear cleanup;
end

function value = final_scalar_or_nan(vec)
    if isempty(vec)
        value = NaN;
    else
        value = vec(end);
    end
end

function x = map_to_physical(xi)
    x = (xi + 1) / 2;
end

function ops = build_b_matrix_operators(N, Q, eval_points, error_quad_points)
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
    [xc, wc] = legendre_gauss(error_quad_points);
    Pc = legendre_values(xc, full_deg);
    dPc = legendre_derivative_values(xc, full_deg);
    Phi_cont = Pc * B3;
    DPhi_cont = dPc * B3;
    [Xq, Yq] = meshgrid(xq, xq);
    [Xe, Ye] = meshgrid(xe, xe);
    [Xc, Yc] = meshgrid(xc, xc);
    Wc = 0.25 * (wc * wc.');

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
    ops.Pphi = Pphi;
    ops.B3 = B3;
    ops.full_deg = full_deg;
    ops.Phi_eval = Phi_eval;
    ops.Phi_cont = Phi_cont;
    ops.DPhi_cont = DPhi_cont;
    ops.Xq = Xq;
    ops.Yq = Yq;
    ops.Xe = Xe;
    ops.Ye = Ye;
    ops.Xc = Xc;
    ops.Yc = Yc;
    ops.Wc = Wc;
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

function dP = legendre_derivative_values(x, max_degree)
    P = legendre_values(x, max_degree);
    dP = zeros(size(P));
    if max_degree == 0
        return;
    end

    dP(:, 2) = 1;
    denom = x .^ 2 - 1;
    for n = 2:max_degree
        dP(:, n + 1) = n * (x .* P(:, n + 1) - P(:, n)) ./ denom;
    end
end

function [z, iter_count, res_norm, status] = newton_solve(z0, F, Jfun, convexity_fun, t_val, problem, tol, maxit, line_search_min)
    z = z0;
    r = F(z);
    res_norm = norm(r, inf);
    if res_norm < tol
        iter_count = 0;
        status = struct('converged', true, 'reason', 'initial_residual_below_tol');
        return;
    end

    current_metrics = convexity_fun(z);

    for iter_count = 1:maxit
        J = Jfun(z);
        delta = -(J + 1e-12 * eye(size(J))) \ r;

        alpha = 1.0;
        accepted = false;
        fallback_found = false;
        fallback_z = z;
        fallback_r = r;
        fallback_res_norm = res_norm;
        fallback_metrics = current_metrics;
        while alpha >= line_search_min
            z_trial = z + alpha * delta;
            r_trial = F(z_trial);
            trial_res_norm = norm(r_trial, inf);
            trial_metrics = convexity_fun(z_trial);
            if trial_res_norm < (1 - 1e-4 * alpha) * res_norm && ...
                    convexity_acceptable(trial_metrics, current_metrics, problem, t_val)
                z = z_trial;
                r = r_trial;
                res_norm = trial_res_norm;
                current_metrics = trial_metrics;
                accepted = true;
                break;
            end
            if problem.line_search_allow_fallback && ...
                    trial_res_norm < problem.line_search_fallback_residual_factor * res_norm && ...
                    convexity_envelope_acceptable(trial_metrics, problem, t_val)
                if ~fallback_found || trial_res_norm < fallback_res_norm
                    fallback_found = true;
                    fallback_z = z_trial;
                    fallback_r = r_trial;
                    fallback_res_norm = trial_res_norm;
                    fallback_metrics = trial_metrics;
                end
            end
            alpha = alpha / 2;
        end

        if ~accepted
            if fallback_found
                z = fallback_z;
                r = fallback_r;
                res_norm = fallback_res_norm;
                current_metrics = fallback_metrics;
            else
                status = struct('converged', false, 'reason', 'line_search_failed');
                return;
            end
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

function [curve_x, tick_positions, tick_labels, display_mask] = homotopy_curve_x_values(t_values)
    positive_mask = t_values > 0;
    zero_mask = t_values == 0;
    display_mask = positive_mask | zero_mask;

    positive_x = -log10(t_values(positive_mask));
    curve_x = positive_x;
    tick_positions = positive_x(:).';
    tick_labels = arrayfun(@(val) sprintf('%.0f', val), tick_positions, 'UniformOutput', false);

    if any(zero_mask)
        if isempty(positive_x)
            zero_x = 1;
        else
            zero_x = max(positive_x) + 2;
        end
        curve_x = [curve_x; zero_x];
        tick_positions = [tick_positions, zero_x];
        tick_labels = [tick_labels, {'t=0'}];
    end
end

function tol = t_value_tolerance(values)
    finite_values = abs(values(isfinite(values)));
    if isempty(finite_values)
        scale = realmin('double');
    else
        scale = max(finite_values);
        if scale == 0
            scale = realmin('double');
        end
    end
    tol = 10 * eps(scale);
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
