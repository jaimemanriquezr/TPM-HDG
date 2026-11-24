function [solution, global_system] = solve(problem, mesh, k, options)
%SOLVE  Solve the variational problem on a mesh.
%
%   solution = solve(problem, mesh, k)
%   solution = solve(problem, mesh, Order = k, Cache = cache)
%   [solution, global_system] = solve(...)
%
%   solution is an FESolution; global_system is a struct with the global
%   (hybrid) matrix, right-hand side, free dofs, solution vector and the
%   residue of the global solve.
%
%   The path is chosen by content. With problem.Operator set (declarative
%   path), the operator is assembled, statically condensed, the boundary
%   conditions are applied and the local fields are recovered as FEFields.
%   With problem.Operator empty, the subclass supplies get_matrices,
%   get_assemble_indices, get_boundary_conditions and get_local_solutions.
%
%   Options:
%     Order  Polynomial degree; the positional k overrides it, and
%            problem.Order is used when neither is given.
%     Cache  Struct whose field indices (with lhs_rows, lhs_cols and
%            rhs_rows) skips recomputing the assemble indices, e.g. in
%            time loops.
%
%See also: FEProblem, FESolution, condense_global_system
arguments
    problem     (1,1)   FEProblem
    mesh            AbstractGrid
    k               double          = []
    options.Order (1,1) double = NaN
    % --- optionally cached quantities: precompute outside a time loop and
    %     pass in to avoid redundant work ---
    options.Cache             struct = struct.empty  % struct with any of: indices
end
% k resolves positional > name-value > problem property
if isempty(k)
    if ~isnan(options.Order)
        k = options.Order;
    else
        k = problem.Order;
    end
end
problem.Order = k;

solution = FESolution(problem, mesh, k);
problem.Mesh = mesh;
% bind the spaces to the solve mesh: a space enumerates its own dofs
if ~isempty(problem.Spaces)
    [problem.Spaces.Mesh] = deal(mesh);
end

if isempty(problem.Operator)
    %% ================================================================== %%
    %% Operator overrides (low-level path): get_matrices etc.
    %% ================================================================== %%
    if isempty(options.Cache) || ~isfield(options.Cache, "indices")
        [lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(problem, mesh);
    else
        [lhs_rows, lhs_cols, rhs_rows] = deal( ...
            options.Cache.indices.lhs_rows, ...
            options.Cache.indices.lhs_cols, ...
            options.Cache.indices.rhs_rows);
    end
    [local_operators, lhs_entries, rhs_entries] = get_matrices(problem, mesh);
    H   = sparse(lhs_rows(:), lhs_cols(:), lhs_entries(:));
    RHS = accumarray(rhs_rows, rhs_entries(:), [size(H, 1), 1]);

    [x_hat, dir_dof, g_hat] = get_boundary_conditions(problem, mesh);
    % TPM transfer of the essential datum (Cockburn & Solano): a class may
    % replace the Dirichlet rows of the transferred faces and free their
    % dofs (see transfer_boundary_rows); a no-op without a mapping
    if ismethod(problem, "apply_transfer_conditions")
        [H, RHS, dir_dof_T] = apply_transfer_conditions(problem, mesh, local_operators, H, RHS, dir_dof);
        x_hat(setdiff(dir_dof, dir_dof_T)) = 0;
        dir_dof = dir_dof_T;
    end
    RHS = RHS - H*x_hat + g_hat;

    [H, RHS] = problem.AdditionalConstraints(problem, mesh, H, RHS);

    freeDOF = setdiff((1:size(H, 1)).', dir_dof);
    x_hat(freeDOF) = H(freeDOF, freeDOF) \ RHS(freeDOF);

    if nargout > 1
        global_system = struct( ...
            "matrix",   H, ...
            "rhs",      RHS, ...
            "dof",      freeDOF, ...
            "solution", x_hat, ...
            "residue",  RHS(freeDOF) - H(freeDOF, freeDOF)*x_hat(freeDOF));
    end

    [local_solutions, trace_solutions] = get_local_solutions(problem, x_hat, local_operators, mesh);
    solution = solution.assign_unknowns(local_solutions, trace_solutions);

else
    %% ================================================================== %%
    %% Declarative path: assemble the Operator, condense, recover
    %% ================================================================== %%
    S = condense_global_system(problem, mesh);
    L = S.L;  H = S.H;  RHS = S.RHS;
    n_el = size(mesh.elements, 2);
    n_g  = numel(problem.GlobalSpaces);
    d_g  = S.d_g;  f_g = S.f_g;
    S12  = S.S12;  S1f = S.S1f;
    dims_loc = arrayfun(@(sp) sp.Dimension, problem.LocalSpaces);

    [x_hat, g_hat, H_hat, dof] = problem.BoundaryConditions.apply(problem);
    % uncondensed local modes (FESpace.UncondensedModes) append
    % element-supported unknowns after the global spaces; they carry no
    % boundary condition and are always free
    n_pad = L.n_dof - numel(x_hat);
    if n_pad > 0
        x_hat = [x_hat; zeros(n_pad, 1)];  g_hat = [g_hat; zeros(n_pad, 1)];
        H_hat = blkdiag(H_hat, sparse(n_pad, n_pad));
        dof.free = [dof.free(:).', L.n_dof - n_pad + (1:n_pad)];
    end

    H(dof.essential.weak, :)   = H_hat(dof.essential.weak, :);
    RHS(dof.essential.weak, :) = 0;
    RHS = RHS - H*x_hat + g_hat;
    x_hat(dof.free) = H(dof.free, dof.free) \ RHS(dof.free);

    % local recovery: x_loc = A11 \ A1f - sum_t (A11 \ A12_t) x_hat_t,
    % gathering each retained space's dofs by its footprint
    x_loc = reshape(S1f, [], 1, n_el);
    for t = 1:numel(S12)   % the global spaces, plus the uncondensed modes
        S12flat = reshape(permute(S12{t}, [1 2 7 5 3 4 6]), size(S12{t}, 1), d_g(t)*f_g(t), n_el);
        u_t = reshape(x_hat(L.footprint{t}), d_g(t)*f_g(t), 1, n_el);
        x_loc = x_loc - pagemtimes(S12flat, u_t);
    end

    U = get_field_functions(problem.Space, problem.Id.field);
    [U.mesh] = deal(mesh);
    x_loc = reshape(x_loc, [], n_el);
    field_dofs = mat2cell(x_loc, dims_loc(:), size(x_loc, 2));
    [U.dof] = field_dofs{:};

    U_hat = get_field_functions(problem.Space, problem.Id.hybrid);
    [U_hat.mesh] = deal(mesh);
    for t = 1:n_g
        U_hat(t).dof = x_hat(L.offset(t) + (1:L.count(t)));
    end

    solution.Unknowns = struct('field', U, 'hybrid', U_hat);
    if nargout > 1
        global_system = struct("matrix", H, "rhs", RHS, "dof", dof.free, ...
            "solution", x_hat, ...
            "residue", norm(H(dof.free, dof.free)*x_hat(dof.free) - RHS(dof.free)));
    end
end
end
