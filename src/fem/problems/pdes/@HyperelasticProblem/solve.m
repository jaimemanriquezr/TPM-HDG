function [solution, global_system, newton] = solve(problem, mesh, k, options)
%SOLVE  Newton's method for the hyperelastic HDG scheme.
%
%   solution = solve(problem, mesh, k)
%   [solution, global_system, newton] = solve(problem, mesh, Order = k, Cache = cache)
%
%   Each Newton step is one low-level solve (FEProblem.solve) with
%   Data.iterate set to the current (u, û): get_matrices returns the
%   linearization written for the new iterate, so the boundary
%   conditions and the static condensation are those of the linear class.
%   The initial iterate is Data.iterate when set, else the state's
%   (displacement, displacement_skeleton) when Data.state holds them, else
%   none (the first step is then the linear solve). The iteration stops
%   when |r| <= newton_tolerance |r_0| (newton_residual at the new
%   iterate; r_0 the linear residual of the zero iterate with the
%   Dirichlet data imposed), or when a residual below 1e-6 no longer
%   halves (the round-off floor of the linear solves); it warns
%   HyperelasticProblem:newtonNotConverged after newton_max_iterations
%   steps.
%
%   solution       FESolution of the last iterate
%   global_system  the global system of the last Newton step
%   newton         struct: iterations, residual (relative residual after
%                  each step, with the τ of the new iterate under the
%                  lagged adaptive stabilization), converged, stagnated
%                  (true when the stop was the round-off floor, not
%                  newton_tolerance), stop_reason ("tolerance",
%                  "stagnation" or "max_iterations"; converged is true for
%                  the first two), tau_max (the largest τ of the initial
%                  iterate and, lagged, of each new iterate)
%
%   With Data.stabilization = "adaptive" the mesh's τ is replaced by
%   adaptive_tau at the initial iterate (zero: the linear value), and,
%   under stabilization_update = "lagged", at every new iterate; τ is
%   not differentiated in the Jacobian.
arguments
    problem     (1,1)   HyperelasticProblem
    mesh            SimplicialGrid
    k               double          = []
    options.Order (1,1) double = NaN
    options.Cache             struct = struct.empty
end
if isempty(k)
    if ~isnan(options.Order)
        k = options.Order;
    else
        k = problem.Order;
    end
end
problem.Order = k;
data = problem.Data;

cache = options.Cache;
if isempty(cache) || ~isfield(cache, "indices")
    indices = struct();
    [indices.lhs_rows, indices.lhs_cols, indices.rhs_rows] = get_assemble_indices(problem, mesh, k);
    cache = struct("indices", indices);
end

N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = 4*d2 + (k > 0)*d1;

% initial iterate
if isempty(data.iterate) && data.density > 0 && ~isempty(data.state) ...
        && isfield(data.state, "solid") && has_kinematics(data.state.solid)
    problem.Data.iterate = struct( ...
        "displacement", data.state.solid.displacement, ...
        "displacement_skeleton", data.state.solid.displacement_skeleton);
end

% adaptive stabilization: τ of the initial iterate on the solve mesh
adaptive = isfield(data, "stabilization") && data.stabilization == "adaptive";
if adaptive
    mesh.tau = adaptive_tau(problem, mesh, k, problem.Data.iterate);
end
tau_max = max(mesh.tau, [], 'all');
% the iterate-independent operators, built once per τ (linear_operators)
problem.Data.newton_operators = linear_operators(problem, mesh, k);

[x_D, dir_dof] = get_boundary_conditions(problem, mesh);
x_ref = zeros(size(x_D));
x_ref(dir_dof) = x_D(dir_dof);
r_ref = norm(newton_residual(problem, mesh, k, zeros(d2peers + 3*d2, N_ELEMENTS), ...
    x_ref, cache.indices, nonlinear = false));
r_ref = max(r_ref, realmin);

history = zeros(1, 0);
converged = false;
stop_reason = "max_iterations";
for iteration = 1:data.newton_max_iterations
    [solution, global_system] = solve@LinearElasticityProblem(problem, mesh, k, Cache = cache);
    x_hat = global_system.solution;
    X = local_dofs(solution);
    problem.Data.iterate = struct( ...
        "displacement", X(d2peers + (1:2*d2), :), ...
        "displacement_skeleton", [reshape(x_hat(1:d1*N_FACES), d1, N_FACES);
                                  reshape(x_hat(d1*N_FACES + 1:end), d1, N_FACES)]);
    if adaptive && data.stabilization_update == "lagged"
        mesh.tau = adaptive_tau(problem, mesh, k, problem.Data.iterate);
        tau_max(end + 1) = max(mesh.tau, [], 'all'); %#ok<AGROW>
        problem.Data.newton_operators = linear_operators(problem, mesh, k);
    end
    history(end + 1) = norm(newton_residual(problem, mesh, k, X, x_hat, cache.indices))/r_ref; %#ok<AGROW>
    stagnated = iteration > 1 && history(end) <= 1e-6 && history(end) > history(end - 1)/2;
    if history(end) <= data.newton_tolerance || stagnated
        converged = true;
        stop_reason = "stagnation";
        if history(end) <= data.newton_tolerance, stop_reason = "tolerance"; end
        break
    end
end
if ~converged
    warning("HyperelasticProblem:newtonNotConverged", ...
        "Newton stopped after %d steps at relative residual %.3e.", ...
        data.newton_max_iterations, history(end));
end
newton = struct("iterations", iteration, "residual", history, "converged", converged, ...
    "stagnated", stop_reason == "stagnation", "stop_reason", stop_reason, "tau_max", tau_max);
solution = without_operators(solution);
end

function solution = without_operators(solution)
% the solution with Data.newton_operators removed from its problem
problem = solution.Problem;
problem.Data = rmfield(problem.Data, "newton_operators");
unknowns = solution.Unknowns;
solution = FESolution(problem, solution.Mesh, solution.Order);
solution.Unknowns = unknowns;
end

function tf = has_kinematics(s)
% the state holds a displacement and its skeleton as arrays
tf = all(isfield(s, ["displacement", "displacement_skeleton"])) ...
    && ~isempty(s.displacement) && ~isequal(s.displacement, 0) ...
    && ~isempty(s.displacement_skeleton) && ~isequal(s.displacement_skeleton, 0);
end

function X = local_dofs(solution)
% [σ; u; γ] dofs, one column per element
dofs = arrayfun(@(f) as_matrix(f.dof), solution.Unknowns.field, 'UniformOutput', false);
X = cat(1, dofs{:});
end

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end
