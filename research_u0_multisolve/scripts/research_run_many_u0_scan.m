function cfg = research_run_many_u0_scan(overrides)
%RESEARCH_RUN_MANY_U0_SCAN Run a batched multi-u0 research scan.
%
%   cfg = research_run_many_u0_scan()
%   cfg = research_run_many_u0_scan(overrides)
%
% The scan driver builds a small library of external u0 candidates, runs
% one multiple-level solve per candidate unless dry_run is enabled, and
% exports a scan-level manifest under the research workspace.

    if nargin < 1 || isempty(overrides)
        overrides = struct();
    end
    if ~isstruct(overrides)
        error('ResearchScan:InvalidOverrides', ...
            'research_run_many_u0_scan expects a struct of runtime overrides.');
    end

    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    solver_dir = fullfile(repo_root, 'solve_u0_anchor_monge_ampere_homotopy');
    workspace_dir = fullfile(repo_root, 'research_u0_multisolve');

    assert(isfolder(solver_dir), ...
        'Could not find solver directory: %s', solver_dir);
    addpath(solver_dir);

    cfg = default_scan_config(workspace_dir);
    cfg = apply_scan_overrides(cfg, overrides);
    cfg.anchor_specs = research_u0_candidate_library(cfg.case_name, cfg.anchor_tags);

    ensure_directory(cfg.results_root);
    ensure_directory(fullfile(cfg.results_root, 'runs'));

    manifest_rows = repmat(empty_scan_row(), numel(cfg.anchor_specs), 1);
    run_results = cell(numel(cfg.anchor_specs), 1);

    fprintf('Running research multi-u0 scan\n');
    fprintf('Case name      : %s\n', cfg.case_name);
    fprintf('Experiment tag : %s\n', cfg.experiment_tag);
    fprintf('Dry run        : %s\n', logical_text(cfg.dry_run));
    fprintf('Acceptance policy : %s\n', cfg.acceptance_policy);
    fprintf('Results root   : %s\n', cfg.results_root);
    fprintf('Anchor count   : %d\n', numel(cfg.anchor_specs));

    for k = 1:numel(cfg.anchor_specs)
        anchor_spec = cfg.anchor_specs(k);
        run_workspace = fullfile(cfg.results_root, 'runs', char(anchor_spec.tag));
        ensure_directory(run_workspace);

        row = build_scan_row(cfg, anchor_spec, run_workspace);
        if ~cfg.dry_run
            [solver_output, row] = run_anchor_case(cfg, anchor_spec, row, run_workspace);
            run_results{k} = solver_output;
        end
        manifest_rows(k) = row;
    end

    manifest_table = struct2table(manifest_rows);
    cfg.scan_manifest = manifest_table;
    cfg.manifest_csv_path = fullfile(cfg.results_root, 'scan_manifest.csv');
    cfg.manifest_mat_path = fullfile(cfg.results_root, 'scan_manifest.mat');
    writetable(manifest_table, cfg.manifest_csv_path);
    scan_manifest = manifest_table; %#ok<NASGU>
    save(cfg.manifest_mat_path, 'scan_manifest', 'cfg');

    if ~cfg.dry_run
        cfg.run_results = run_results;
    end

    if nargout == 0
        clear cfg
    end
end

function cfg = default_scan_config(workspace_dir)
    cfg = struct();
    cfg.case_name = 'paper_example_5_1';
    cfg.experiment_tag = 'research_baseline_multi_u0_scan';
    cfg.results_root = fullfile(workspace_dir, 'results', cfg.experiment_tag);
    cfg.dry_run = false;
    cfg.show_figures = false;
    cfg.anchor_tags = strings(0, 1);
    cfg.acceptance_policy = "diagnostic_only";
    cfg.multiple_level_options = struct( ...
        'N_values', [6, 10, 14], ...
        'p_values', [3, 5, 7], ...
        'show_figures', false, ...
        'acceptance_policy', char(cfg.acceptance_policy));
end

function cfg = apply_scan_overrides(cfg, overrides)
    override_fields = {'case_name', 'experiment_tag', 'dry_run', 'show_figures', ...
        'results_root', 'anchor_tags', 'multiple_level_options', 'acceptance_policy'};
    for k = 1:numel(override_fields)
        field_name = override_fields{k};
        if isfield(overrides, field_name) && ~isempty(overrides.(field_name))
            cfg.(field_name) = overrides.(field_name);
        end
    end

    if isfield(cfg, 'case_name')
        cfg.case_name = char(string(cfg.case_name));
    end
    cfg.experiment_tag = char(string(cfg.experiment_tag));
    cfg.results_root = char(string(cfg.results_root));
    cfg.dry_run = logical(cfg.dry_run);
    cfg.show_figures = logical(cfg.show_figures);
    cfg.anchor_tags = string(cfg.anchor_tags(:));
    cfg.acceptance_policy = string(cfg.acceptance_policy);
    if ~isfield(cfg.multiple_level_options, 'acceptance_policy') || ...
            isempty(cfg.multiple_level_options.acceptance_policy)
        cfg.multiple_level_options.acceptance_policy = char(cfg.acceptance_policy);
    else
        cfg.acceptance_policy = string(cfg.multiple_level_options.acceptance_policy);
    end
    if isfield(cfg.multiple_level_options, 'show_figures')
        cfg.multiple_level_options.show_figures = logical(cfg.multiple_level_options.show_figures);
    end
end

function row = build_scan_row(cfg, anchor_spec, run_workspace)
    row = empty_scan_row();
    row.CaseName = string(cfg.case_name);
    row.ExperimentTag = string(cfg.experiment_tag);
    row.AnchorTag = string(anchor_spec.tag);
    row.AnchorFamily = string(anchor_spec.family);
    row.AnchorLabel = string(anchor_spec.anchor_label);
    row.RunTag = string(anchor_spec.tag);
    row.RunDirectory = string(normalize_path(run_workspace));
    row.LevelCount = numel(cfg.multiple_level_options.N_values);
    row.FinalLevelN = cfg.multiple_level_options.N_values(end);
    row.AcceptancePolicy = string(cfg.multiple_level_options.acceptance_policy);
    row.ConvexityMode = "";
    row.LineSearchFallbackEnabled = false;
    row.Success = false;
    row.ReachedTEnd = false;
    row.FinalT = NaN;
    row.FinalResidualInf = NaN;
    row.TerminalStatePath = "";
    row.TerminalMetricsPath = "";
    row.U0DefinitionPath = "";
end

function [solver_output, row] = run_anchor_case(cfg, anchor_spec, row, run_workspace)
    solver_options = cfg.multiple_level_options;
    solver_options.show_figures = cfg.show_figures;
    solver_options.run_tag = char(anchor_spec.tag);
    if isa(anchor_spec.u0_fun, 'function_handle')
        solver_options.u0_fun = anchor_spec.u0_fun;
        solver_options.anchor_label = char(anchor_spec.anchor_label);
    end

    original_dir = pwd;
    cleanup = onCleanup(@() cd(original_dir));
    cd(run_workspace);
    solver_output = u0_anchor_monge_ampere_multiple_level_demo(cfg.case_name, solver_options);
    clear cleanup;
    cd(original_dir);

    row.Success = all(cellfun(@(r) r.success, solver_output.level_results));
    row.ReachedTEnd = all(cellfun(@(r) r.reached_t_end, solver_output.level_results));
    row.LevelCount = height(solver_output.level_summary);
    row.FinalLevelN = solver_output.level_summary.N(end);
    row.FinalT = solver_output.level_summary.FinalT(end);
    row.FinalResidualInf = solver_output.level_summary.FinalResidualInf(end);

    [terminal_state_path, terminal_metrics_path, u0_definition_path] = ...
        infer_terminal_export_paths(solver_output);
    row.TerminalStatePath = string(normalize_path(terminal_state_path));
    row.TerminalMetricsPath = string(normalize_path(terminal_metrics_path));
    row.U0DefinitionPath = string(normalize_path(u0_definition_path));

    terminal_state = load(terminal_state_path);
    if isfield(terminal_state, 'acceptance_policy')
        row.AcceptancePolicy = string(terminal_state.acceptance_policy);
    end
    if isfield(terminal_state, 'convexity_mode')
        row.ConvexityMode = string(terminal_state.convexity_mode);
    end
    if isfield(terminal_state, 'line_search_allow_fallback')
        row.LineSearchFallbackEnabled = logical(terminal_state.line_search_allow_fallback);
    end
end

function [terminal_state_path, terminal_metrics_path, u0_definition_path] = ...
        infer_terminal_export_paths(solver_output)
    last_row = solver_output.level_summary(end, :);
    final_level_dir = fullfile(solver_output.output_dir, 'per_level', ...
        sprintf('level_%02d_N_%03d', last_row.Level, last_row.N));
    terminal_state_path = fullfile(final_level_dir, 'terminal_state.mat');
    terminal_metrics_path = fullfile(final_level_dir, 'terminal_metrics.csv');
    u0_definition_path = fullfile(final_level_dir, 'u0_definition.txt');
end

function ensure_directory(path_str)
    if ~exist(path_str, 'dir')
        mkdir(path_str);
    end
end

function path_str = normalize_path(path_str)
    path_str = strrep(path_str, '\', '/');
end

function row = empty_scan_row()
    row = struct( ...
        'CaseName', "", ...
        'ExperimentTag', "", ...
        'AnchorTag', "", ...
        'AnchorFamily', "", ...
        'AnchorLabel', "", ...
        'RunTag', "", ...
        'RunDirectory', "", ...
        'LevelCount', NaN, ...
        'FinalLevelN', NaN, ...
        'AcceptancePolicy', "", ...
        'ConvexityMode', "", ...
        'LineSearchFallbackEnabled', false, ...
        'Success', false, ...
        'ReachedTEnd', false, ...
        'FinalT', NaN, ...
        'FinalResidualInf', NaN, ...
        'TerminalStatePath', "", ...
        'TerminalMetricsPath', "", ...
        'U0DefinitionPath', "");
end

function txt = logical_text(tf)
    if tf
        txt = 'true';
    else
        txt = 'false';
    end
end
