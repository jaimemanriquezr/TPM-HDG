function run_suite()
% RUN_SUITE Run every tests/test_*.m and report failures at the end.
% A test fails if it throws or prints "N check(s) failed" with N > 0.
% Assumes the toolbox is on the path (installed, or initpath run). Each
% test runs in its own function workspace so the scripts' clc/clear cannot
% wipe the loop state.
root = fileparts(fileparts(mfilename('fullpath')));
tests = dir(fullfile(root, 'tests', 'test_*.m'));
n_fail = 0; failed = {};
for i = 1:numel(tests)
    name = tests(i).name;
    fprintf('=== %s ===\n', name);
    try
        out = evalc('run_one(fullfile(root, ''tests'', name))');
        fprintf('%s', out);
        counts = regexp(out, '(\d+) check\(s\) failed', 'tokens');
        if any(cellfun(@(c) str2double(c{1}), counts) > 0)
            error('run_suite:checksFailed', 'the script reported failed checks');
        end
    catch err
        fprintf('FAILED: %s\n', err.message);
        n_fail = n_fail + 1; failed{end+1} = name; %#ok<AGROW>
    end
    cd(root);
end
fprintf('SUITE: %d/%d failed\n', n_fail, numel(tests));
if n_fail > 0, disp(failed); end
end

function run_one(script_path)
run(script_path);
end
