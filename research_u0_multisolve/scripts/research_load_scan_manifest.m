function manifest = research_load_scan_manifest(scan_root)
%RESEARCH_LOAD_SCAN_MANIFEST Load the scan manifest from a research result root.

    manifest_csv = fullfile(scan_root, 'scan_manifest.csv');
    assert(isfolder(scan_root), ...
        'Results directory does not exist: %s', scan_root);
    assert(isfile(manifest_csv), ...
        'Missing scan manifest: %s', manifest_csv);

    raw = readcell(manifest_csv, 'Delimiter', ',');
    variable_names = matlab.lang.makeValidName(string(raw(1, :)));
    data = raw(2:end, :);
    manifest = cell2table(data, 'VariableNames', cellstr(variable_names));

    string_columns = {'CaseName', 'ExperimentTag', 'AnchorTag', 'AnchorFamily', ...
        'AnchorLabel', 'RunTag', 'RunDirectory', 'AcceptancePolicy', ...
        'ConvexityMode', 'TerminalStatePath', ...
        'TerminalMetricsPath', 'U0DefinitionPath'};
    numeric_columns = {'LevelCount', 'FinalLevelN', 'FinalT', 'FinalResidualInf'};
    logical_columns = {'Success', 'ReachedTEnd', 'LineSearchFallbackEnabled'};

    for k = 1:numel(string_columns)
        manifest.(string_columns{k}) = string(manifest.(string_columns{k}));
    end
    for k = 1:numel(numeric_columns)
        manifest.(numeric_columns{k}) = unwrap_numeric_column(manifest.(numeric_columns{k}));
    end
    for k = 1:numel(logical_columns)
        manifest.(logical_columns{k}) = logical(unwrap_numeric_column(manifest.(logical_columns{k})));
    end
end

function values = unwrap_numeric_column(values)
    if iscell(values)
        values = cellfun(@double, values);
    else
        values = double(values);
    end
end
