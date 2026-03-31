function u0_anchor_multiple_level_tiny_target_t_test()
% Regression test for tiny target t values in the multiple-level driver.

    opts = struct( ...
        'N_values', [2, 4], ...
        'p_values', [15, 19], ...
        'run_tag', 'multiple_level_tiny_target_t_test', ...
        'show_figures', false);

    out = u0_anchor_monge_ampere_multiple_level_demo('paper_example_5_1', opts);
    T = out.level_summary;

    assert(isequal(T.N.', [2, 4]), 'Expected stage degrees N = [2, 4].');
    assert(isequal(T.TargetP.', [15, 19]), 'Expected paired p = [15, 19].');
    assert(all(T.TargetT > 0), 'Tiny target t values must remain strictly positive.');
    assert(all(abs(T.TargetT - [1e-15; 1e-19]) < [1e-27; 1e-31]), ...
        'Expected tiny target t values to be preserved exactly in the summary.');
    assert(all(T.ReachedTEnd), 'Each tiny-target level should reach its requested t.');
end
