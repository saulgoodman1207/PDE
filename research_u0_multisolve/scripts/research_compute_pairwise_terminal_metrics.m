function pairwise_metrics = research_compute_pairwise_terminal_metrics(inventory_table, inventory_states)
%RESEARCH_COMPUTE_PAIRWISE_TERMINAL_METRICS Compute unordered pairwise metrics.

    num_states = numel(inventory_states);
    if num_states < 2
        pairwise_metrics = struct2table(repmat(empty_pairwise_row(), 0, 1));
        return;
    end

    pair_count = nchoosek(num_states, 2);
    rows = repmat(empty_pairwise_row(), pair_count, 1);
    row_idx = 0;

    for i = 1:num_states-1
        state_a = inventory_states{i};
        for j = i+1:num_states
            state_b = inventory_states{j};
            row_idx = row_idx + 1;
            rows(row_idx) = pairwise_row(inventory_table(i, :), inventory_table(j, :), state_a, state_b);
        end
    end

    pairwise_metrics = struct2table(rows);
end

function row = pairwise_row(inv_a, inv_b, state_a, state_b)
    comparable_grid = isequal(size(state_a.u_num_eval), size(state_b.u_num_eval)) ...
        && isequal(size(state_a.x_eval), size(state_b.x_eval)) ...
        && isequal(size(state_a.y_eval), size(state_b.y_eval));
    comparable_coeffs = numel(state_a.U_final) == numel(state_b.U_final);

    row = empty_pairwise_row();
    row.AnchorTagA = inv_a.AnchorTag(1);
    row.AnchorTagB = inv_b.AnchorTag(1);
    row.AnchorFamilyA = inv_a.AnchorFamily(1);
    row.AnchorFamilyB = inv_b.AnchorFamily(1);
    row.AcceptancePolicyA = inv_a.AcceptancePolicy(1);
    row.AcceptancePolicyB = inv_b.AcceptancePolicy(1);
    row.ConvexityModeA = inv_a.ConvexityMode(1);
    row.ConvexityModeB = inv_b.ConvexityMode(1);
    row.LineSearchFallbackA = inv_a.LineSearchFallbackEnabled(1);
    row.LineSearchFallbackB = inv_b.LineSearchFallbackEnabled(1);
    row.AcceptancePolicyMatch = inv_a.AcceptancePolicy(1) == inv_b.AcceptancePolicy(1);
    row.ConvexityModeMatch = inv_a.ConvexityMode(1) == inv_b.ConvexityMode(1);
    row.LineSearchFallbackMatch = inv_a.LineSearchFallbackEnabled(1) == inv_b.LineSearchFallbackEnabled(1);
    row.ComparableGrid = comparable_grid;
    row.ComparableCoefficients = comparable_coeffs;
    row.SuccessA = inv_a.Success(1);
    row.SuccessB = inv_b.Success(1);
    row.ReachedTEndA = inv_a.ReachedTEnd(1);
    row.ReachedTEndB = inv_b.ReachedTEnd(1);

    if comparable_grid
        diff_u = state_a.u_num_eval(:) - state_b.u_num_eval(:);
        row.AbsUInf = max(abs(diff_u));
        row.RelUFro = norm(diff_u) / max(norm(state_a.u_num_eval(:)), eps);
    end

    if comparable_coeffs
        diff_coeffs = state_a.U_final(:) - state_b.U_final(:);
        row.CoeffInf = max(abs(diff_coeffs));
    end

    row.AbsResidualDelta = abs(inv_a.FinalResidualInf(1) - inv_b.FinalResidualInf(1));
    row.AbsJacobianCondDelta = abs(inv_a.FinalJacobianCond(1) - inv_b.FinalJacobianCond(1));
    row.AbsPSDRatioDelta = abs(inv_a.FinalPSDRatio(1) - inv_b.FinalPSDRatio(1));
    row.AbsIndefRatioDelta = abs(inv_a.FinalIndefRatio(1) - inv_b.FinalIndefRatio(1));
end

function row = empty_pairwise_row()
    row = struct( ...
        'AnchorTagA', "", ...
        'AnchorTagB', "", ...
        'AnchorFamilyA', "", ...
        'AnchorFamilyB', "", ...
        'AcceptancePolicyA', "", ...
        'AcceptancePolicyB', "", ...
        'ConvexityModeA', "", ...
        'ConvexityModeB', "", ...
        'LineSearchFallbackA', false, ...
        'LineSearchFallbackB', false, ...
        'AcceptancePolicyMatch', false, ...
        'ConvexityModeMatch', false, ...
        'LineSearchFallbackMatch', false, ...
        'ComparableGrid', false, ...
        'ComparableCoefficients', false, ...
        'SuccessA', false, ...
        'SuccessB', false, ...
        'ReachedTEndA', false, ...
        'ReachedTEndB', false, ...
        'AbsUInf', NaN, ...
        'RelUFro', NaN, ...
        'CoeffInf', NaN, ...
        'AbsResidualDelta', NaN, ...
        'AbsJacobianCondDelta', NaN, ...
        'AbsPSDRatioDelta', NaN, ...
        'AbsIndefRatioDelta', NaN);
end
