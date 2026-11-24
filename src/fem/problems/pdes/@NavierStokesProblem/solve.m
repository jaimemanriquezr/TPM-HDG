function [solution, global_system, history, u_star] = solve(problem, mesh, k, options)
%SOLVE  Navier-Stokes by Picard iteration (one time step when density > 0).
%
%   solution = solve(problem, mesh, k)
%   solution = solve(problem, mesh, Order = k, Cache = cache)
%   [solution, global_system, history, u_star] = solve(...)
%
%   The first step uses Data.convective_velocity of the problem as beta
%   (empty: a Stokes step); each following step is an Oseen problem with
%   beta the postprocessed velocity u* of the previous iterate, passed
%   through Data.convective_velocity. The iteration stops when the
%   relative velocity increment drops below tolerance or after
%   max_iterations steps.
%
%   solution and global_system are those of the last step. history holds
%   the relative increments ||u_{n+1} - u_n|| / ||u_n|| of the Oseen
%   steps (the first step has none). u_star is the postprocessed
%   (divergence-free) velocity of the last iterate.
%
%   Options:
%     Order           Polynomial degree; the positional k overrides it,
%                     and problem.Order is used when neither is given.
%     Cache           Struct whose field indices (lhs_rows, lhs_cols,
%                     rhs_rows) skips recomputing the assemble indices.
%     tolerance       Stopping tolerance on ||u_{n+1} - u_n|| / ||u_n||
%                     (default 1e-8).
%     max_iterations  Maximum number of steps (default 50).
%
%See also: OseenProblem, solve_navier_stokes
arguments
    problem     (1,1)   FEProblem
    mesh            SimplicialGrid
    k               double          = []
    options.Order (1,1) double = NaN
    % --- optionally cached quantities: precompute outside a time loop and
    %     pass in to avoid redundant work ---
    options.Cache             struct = struct.empty  % struct with any of: indices
    % --- Picard iteration (beta = u* of the previous iterate) ---
    options.tolerance      (1,1) double = 1e-8   % on ||u_{n+1} - u_n|| / ||u_n||
    options.max_iterations (1,1) double = 50
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

if isempty(options.Cache) || ~isfield(options.Cache, "indices")
    [lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(problem, mesh);
else
    [lhs_rows, lhs_cols, rhs_rows] = deal( ...
        options.Cache.indices.lhs_rows, ...
        options.Cache.indices.lhs_cols, ...
        options.Cache.indices.rhs_rows);
end

% Picard: get_matrices adds the Oseen term of Data.convective_velocity.
% The first step uses the problem's own convective velocity as the
% starting beta (empty: Stokes); then beta = u* of the previous iterate.
beta = problem.Data.convective_velocity;
u_old = [];
increment = Inf;
history = [];   % increments of the Oseen steps (the Stokes step has none)
n = 0;
% the boundary data do not depend on beta; the loop edits copies
[x_hat_bc, dir_dof_bc, g_hat] = get_boundary_conditions(problem, mesh);

while increment >= options.tolerance && n < options.max_iterations
    n = n + 1;
    problem.Data.convective_velocity = beta;

    [local_operators, lhs_entries, rhs_entries] = get_matrices(problem, mesh);
    H   = sparse(lhs_rows(:), lhs_cols(:), lhs_entries(:));
    RHS = accumarray(rhs_rows, rhs_entries(:));

    x_hat = x_hat_bc;
    dir_dof = dir_dof_bc;
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

    % the postprocessed u* of this iterate convects the next one
    beta = divergence_free_velocity(solution);
    [increment, u_old] = picard_increment(solution, u_old);
    if n > 1, history(end + 1) = increment; end %#ok<AGROW>
end
u_star = beta;
end

function [increment, u] = picard_increment(solution, u_old)
% relative change of the HDG velocity from the previous iterate (Inf first)
u = solution.Unknowns.field(2).dof;
increment = Inf;
if ~isempty(u_old)
    increment = norm(u - u_old, "fro")/norm(u_old, "fro");
end
end
