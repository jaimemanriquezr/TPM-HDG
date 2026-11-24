function [solution, global_system, newton] = solve(problem, mesh, k, options)
%SOLVE  Newton's method (with optional load stepping) for the
%gradient-based hyperelastic HDG scheme.
%
%   solution = solve(problem, mesh, k)
%   [solution, global_system, newton] = solve(problem, mesh, Order = k, Cache = cache)
%
%   Each Newton step is one low-level solve (FEProblem.solve) with
%   Data.iterate set to the current iterate: get_matrices returns the
%   linearization written for the new iterate, so (H, u) are condensed
%   onto û and the boundary conditions are those of the linear class.
%   The initial iterate is Data.iterate, else the state's (gradient,
%   displacement, displacement_skeleton) under Newmark, else none (the
%   first step is then the linear solve). τ is set on the mesh from
%   Data.stabilization at the initial iterate. For "adaptive" and
%   "cockburn_shen", stabilization_update = "lagged" recomputes it after
%   every step (not differentiated); "frozen" keeps it during a Newton
%   solve and, once converged, repeats the solve with τ of the solution
%   until τ no longer changes (at most 4 passes; not converged otherwise). The iteration stops
%   when |r| <= newton_tolerance |r_0| (newton_residual; r_0 the residual
%   of the zero iterate with the Dirichlet data imposed), or when a
%   residual below 1e-6 no longer halves (round-off floor); an iterate
%   with an inverted element (det F <= 0) ends it unconverged.
%   With load_steps = n > 1 (steady problems), the source and the
%   Dirichlet and Neumann data are scaled by the fractions 1/n, 2/n, ...,
%   each solved from the previous one; a fraction whose Newton fails is
%   halved (down to 1/(64 n)); component-wise specs in
%   Data.boundary_conditions are not scaled and raise an error.
%
%   newton  struct: iterations (total), residual (relative residuals of
%           the last fraction), converged, stagnated (true when the stop
%           was the round-off floor, not newton_tolerance), stop_reason
%           ("tolerance", "stagnation", "max_iterations" or "failure": an
%           inverted or non-finite iterate, a frozen τ that does not
%           settle, or load stepping short of the full load; converged is
%           true for the first two), tau_max (largest τ after each
%           step of the last fraction), fractions (the load fractions
%           solved), tau (the final τ, 3 x N) and marked (the elements
%           of T_h^(-,nb) under "cockburn_shen")
arguments
    problem     (1,1)   HyperelasticGradientProblem
    mesh            SimplicialGrid
    k               double          = []
    options.Order (1,1) double = NaN
    options.Cache             struct = struct.empty
end
if isempty(k)
    if ~isnan(options.Order), k = options.Order; else, k = problem.Order; end
end
problem.Order = k;
data = problem.Data;
cache = options.Cache;
if isempty(cache) || ~isfield(cache, "indices")
    indices = struct();
    [indices.lhs_rows, indices.lhs_cols, indices.rhs_rows] = get_assemble_indices(problem, mesh, k);
    cache = struct("indices", indices);
end

if isempty(data.iterate) && data.density > 0 && ~isempty(data.state) ...
        && isfield(data.state, "solid") && isfield(data.state.solid, "gradient") ...
        && ~isempty(data.state.solid.gradient) && ~isequal(data.state.solid.gradient, 0)
    s = data.state.solid;
    problem.Data.iterate = struct("gradient", s.gradient, "displacement", s.displacement, ...
        "displacement_skeleton", s.displacement_skeleton);
end

n_steps = 1;
if data.density == 0, n_steps = data.load_steps; end
if n_steps > 1
    if isfield(data, "boundary_conditions") && ~isempty(data.boundary_conditions)
        error("HyperelasticGradientProblem:loadSteps", "load_steps > 1 scales datum_dirichlet and datum_neumann only; Data.boundary_conditions is not supported.");
    end
    [solution, global_system, newton] = load_stepping(problem, mesh, k, cache, n_steps);
    return
end

% Newton
data = problem.Data;

N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
H0 = [];
if ~isempty(data.iterate), H0 = data.iterate.gradient; end
constants = [];
[mesh.tau, marked, constants] = stabilization_tau(problem, mesh, k, H0, constants);
tau_max = max(mesh.tau, [], 'all');

[x_D, dir_dof] = get_boundary_conditions(problem, mesh);
x_ref = zeros(size(x_D));
x_ref(dir_dof) = x_D(dir_dof);
% the assembly at the zero gradient serves r_0 and a first step from rest
% (newton_operators); after each step, the one at the new iterate serves
% its residual and the next step
problem.Data.newton_operators = newton_operators(problem, mesh, k, zeros(4*d2, N_ELEMENTS));
r_ref = max(norm(newton_residual(problem, mesh, k, zeros(6*d2, N_ELEMENTS), x_ref, cache.indices)), realmin);

adaptive = any(data.stabilization == ["adaptive", "cockburn_shen"]);
lagged = adaptive && data.stabilization_update == "lagged";
passes = 1;
if adaptive && data.stabilization_update == "frozen", passes = 4; end
history = zeros(1, 0);
iterations = 0;
for pass = 1:passes
    converged = false;
    stop_reason = "max_iterations";
    for iteration = 1:data.newton_max_iterations
        iterations = iterations + 1;
        [solution, global_system] = solve@FEProblem(problem, mesh, k, Cache = cache);
        x_hat = global_system.solution;
        X = [solution.Unknowns.field(1).dof; as_matrix(solution.Unknowns.field(2).dof)];
        problem.Data.iterate = struct("gradient", X(1:4*d2, :), "displacement", X(4*d2 + 1:end, :), ...
            "displacement_skeleton", [reshape(x_hat(1:d1*N_FACES), d1, N_FACES);
                                      reshape(x_hat(d1*N_FACES + 1:end), d1, N_FACES)]);
        try
            if lagged
                [mesh.tau, marked, constants] = stabilization_tau(problem, mesh, k, X(1:4*d2, :), constants);
                tau_max(end + 1) = max(mesh.tau, [], 'all'); %#ok<AGROW>
            end
            problem.Data.newton_operators = newton_operators(problem, mesh, k, X(1:4*d2, :));
            history(end + 1) = norm(newton_residual(problem, mesh, k, X, x_hat, cache.indices))/r_ref; %#ok<AGROW>
        catch err
            % an inverted element (det F <= 0) at the iterate ends the iteration
            if err.identifier ~= "HyperelasticProblem:inverted", rethrow(err); end
            history(end + 1) = NaN; %#ok<AGROW>
            stop_reason = "failure";
            break
        end
        if ~all(isfinite(X), 'all'), stop_reason = "failure"; break, end
        stagnated = iteration > 1 && history(end) <= 1e-6 && history(end) > history(end - 1)/2;
        if history(end) <= data.newton_tolerance || stagnated
            converged = true;
            stop_reason = "stagnation";
            if history(end) <= data.newton_tolerance, stop_reason = "tolerance"; end
            break
        end
    end
    if passes == 1 || ~converged, break, end
    % frozen τ: a further pass while τ of the solution differs from the τ used
    [tau_new, marked_new, constants] = stabilization_tau(problem, mesh, k, X(1:4*d2, :), constants);
    if max(abs(tau_new - mesh.tau), [], 'all') <= 1e-8*max(mesh.tau, [], 'all'), break, end
    if pass == passes
        converged = false;
        stop_reason = "failure";
        break
    end
    mesh.tau = tau_new;
    marked = marked_new;
    tau_max(end + 1) = max(mesh.tau, [], 'all'); %#ok<AGROW>
    r_ref = max(norm(newton_residual(problem, mesh, k, zeros(6*d2, N_ELEMENTS), x_ref, cache.indices)), realmin);
end
iteration = iterations;
if ~converged
    warning("HyperelasticGradientProblem:newtonNotConverged", ...
        "Newton stopped after %d steps at relative residual %.3e.", iteration, history(end));
end
newton = struct("iterations", iteration, "residual", history, "converged", converged, ...
    "stagnated", stop_reason == "stagnation", "stop_reason", stop_reason, "tau_max", tau_max, "iterate", problem.Data.iterate, "fractions", 1, ...
    "marked", marked, "tau", mesh.tau);
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

function A = as_matrix(A)
if iscell(A), A = cat(1, A{:}); end
end

function [solution, global_system, newton] = load_stepping(problem, mesh, k, cache, n_steps)
% the load fractions of solve, each a Newton solve from the previous one
full = problem.Data;
scale = @(h, s) {@(x, y) s*h{1}(x, y); @(x, y) s*h{2}(x, y)};
fraction = 0; step = 1/n_steps; total = 0; fractions = [];
warning_state = warning("off", "HyperelasticGradientProblem:newtonNotConverged");
restore = onCleanup(@() warning(warning_state));
while fraction < 1 - 1e-12
    target = min(1, fraction + step);
    trial = problem;
    trial.Data.load_steps = 1;
    trial.Data.source_term = scale(full.source_term, target);
    trial.Data.datum_dirichlet = scale(full.datum_dirichlet, target);
    trial.Data.datum_neumann = scale(full.datum_neumann, target);
    try
        [solution, global_system, newton] = solve(trial, mesh, k, Cache = cache);
        ok = newton.converged;
    catch err
        if ~startsWith(err.identifier, "HyperelasticProblem:"), rethrow(err); end
        ok = false;
    end
    if ok
        total = total + newton.iterations;
        fraction = target;
        fractions(end + 1) = target; %#ok<AGROW>
        problem.Data.iterate = newton.iterate;
    else
        step = step/2;
        if step < 1/(64*n_steps)
            break
        end
    end
end
delete(restore);
if fraction < 1 - 1e-12
    warning("HyperelasticGradientProblem:newtonNotConverged", ...
        "Load stepping stopped at the load fraction %.4f.", fraction);
    newton.converged = false;
    newton.stagnated = false;
    newton.stop_reason = "failure";
end
newton.iterations = total;
newton.fractions = fractions;
end
