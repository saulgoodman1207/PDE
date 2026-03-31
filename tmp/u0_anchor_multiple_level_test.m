function u0_anchor_multiple_level_test()
% Regression test for the new multiple-level u0-anchor driver.

    opts = struct( ...
        'N_values', [2, 4], ...
        'p_values', [3, 5], ...
        'run_tag', 'multiple_level_test', ...
        'show_figures', false);

    out = u0_anchor_monge_ampere_multiple_level_demo('paper_example_5_1', opts);
    T = out.level_summary;

    assert(isequal(T.Level.', [1, 2]), 'Expected two multiple-level stages.');
    assert(isequal(T.N.', [2, 4]), 'Expected stage degrees N = [2, 4].');
    assert(isequal(T.TargetP.', [3, 5]), 'Expected paired p = [3, 5].');
    assert(all(abs(T.StartT - [1; 1e-3]) < 1e-20), ...
        'Expected the second level to start from the previous level target t.');
    assert(all(abs(T.TargetT - [1e-3; 1e-5]) < 1e-20), ...
        'Expected stage targets t = [1e-3, 1e-5].');
    assert(all(T.Success), 'The minimal multiple-level paper_example_5_1 test should succeed.');
end
