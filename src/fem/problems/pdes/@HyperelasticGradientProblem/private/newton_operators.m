function operators = newton_operators(pde, mesh, k, H)
% The assembly of get_matrices at the iterate gradient H and mesh.tau:
% the local operators and the condensed system. solve keeps it in
% Data.newton_operators, where newton_residual (the residual at H) and
% the next Newton step (the linearization at H) take it instead of
% assembling it again.
pde.Data.iterate = struct("gradient", H);
pde.Data.newton_operators = [];
[local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, k);
operators = struct("gradient", H, "tau", mesh.tau, "local_operators", {local_operators}, ...
    "lhs_entries", lhs_entries, "rhs_entries", rhs_entries);
end
