%-------------------------------------------------------------------------
%Description: DarcyProblem.get_inverse_tensor against inv.
%
%             The inverse of the diffusion tensor is checked for a general,
%             a symmetric and a diagonal kappa, both as a constant matrix
%             and as a cell of function handles (column-major cell, as the
%             constant form).
%-------------------------------------------------------------------------
n_failed = 0;

x = [0.1 0.4; 0.7 0.2];
y = [0.3 0.9; 0.5 0.6];
tensors = {[2 .5; .3 1], [2 .5; .5 1], [3 0; 0 .5]};
names = ["general", "symmetric", "diagonal"];

for i = 1:numel(tensors)
    K = tensors{i};
    ref = inv(K);

    k_I = DarcyProblem.get_inverse_tensor(K, x, y);
    err = max(cellfun(@(c, r) max(abs(c - r), [], "all"), k_I, num2cell(ref)), [], "all");
    n_failed = n_failed + check(names(i) + ": constant", err < 1e-14);

    K_fun = arrayfun(@(c) @(x, y) c + 0*x, K, UniformOutput=false);
    k_I = DarcyProblem.get_inverse_tensor(K_fun, x, y);
    err = max(cellfun(@(c, r) max(abs(c - r), [], "all"), k_I, num2cell(ref)), [], "all");
    n_failed = n_failed + check(names(i) + ": function handles", err < 1e-14);
end

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
