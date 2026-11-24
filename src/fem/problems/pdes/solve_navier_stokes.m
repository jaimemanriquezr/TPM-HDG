function [solution, history, u_star] = solve_navier_stokes(problem, mesh, k, options)
%-------------------------------------------------------------------------
%Description: Steady Navier–Stokes by Picard iteration on the Oseen
%             solver (Solano & Vargas 2022, §6.4): the first solve is
%             Stokes (β = 0); then β is set from the previous iterate until
%             ||u_{n+1} - u_n|| / ||u_n|| < tol. With beta = "postprocessed"
%             (the paper's choice) β is the divergence-free u* of the
%             previous iterate (divergence_free_velocity); with "velocity"
%             it is the plain HDG velocity u_h.
%             problem is an OseenProblem or a NavierStokesProblem. An
%             OseenProblem's convective_velocity is ignored (the loop owns
%             β); a low-level NavierStokesProblem runs its own Picard loop
%             (NavierStokesProblem.solve), whose first step convects with
%             its Data.convective_velocity (empty = Stokes).
%             On the declarative path (problem.Operator set) pass
%             operator_of, a handle beta -> AffineOperator rebuilding the
%             Oseen operator (e.g. from StokesMatrices(..., beta)); the
%             boundary conditions on the problem are kept.
%Output:    - solution: the FESolution of the last Oseen solve
%           - history:  relative increments per iteration (1 x n_iter)
%           - u_star:   {u1*; u2*} of the last solve ("postprocessed" only)
%-------------------------------------------------------------------------
arguments
    problem (1,1) FEProblem
    mesh SimplicialGrid
    k (1,1) double = problem.Order
    options.beta (1,1) string {mustBeMember(options.beta, ["postprocessed", "velocity"])} = "postprocessed"
    options.tolerance (1,1) double = 1e-8
    options.max_iterations (1,1) double = 50
    options.verbose (1,1) logical = false
    options.operator_of function_handle = function_handle.empty   % declarative path
end
declarative = ~isempty(options.operator_of);
% a NavierStokesProblem runs its own Picard loop (NavierStokesProblem.solve),
% on whichever path the problem is set up for
if isa(problem, "NavierStokesProblem") && ~declarative
    assert(options.beta == "postprocessed", "solve_navier_stokes:beta", ...
        "NavierStokesProblem.solve always convects with the postprocessed u*");
    [solution, ~, history, u_star] = solve(problem, mesh, k, ...
        tolerance = options.tolerance, max_iterations = options.max_iterations);
    if options.verbose
        fprintf("  Picard %2d: relative increment %.3e\n", [1:numel(history); history]);
    end
    return
end
assert(declarative || isempty(problem.Operator), "solve_navier_stokes:operator_of", ...
    "a declarative problem needs operator_of(beta) to rebuild its Oseen operator");
problem = with_beta(problem, {});
solution = solve(problem, mesh, Order= k);
[beta, u_dof] = next_beta(solution, options.beta);

u_star = {};
history = zeros(1, options.max_iterations);
for n = 1:options.max_iterations
    problem = with_beta(problem, beta);
    solution = solve(problem, mesh, Order= k);
    [beta, u_new] = next_beta(solution, options.beta);
    history(n) = norm(u_new - u_dof, 'fro')/norm(u_dof, 'fro');
    if options.verbose
        fprintf("  Picard %2d: relative increment %.3e\n", n, history(n));
    end
    u_dof = u_new;
    if history(n) < options.tolerance
        history = history(1:n);
        if options.beta == "postprocessed"; u_star = beta; end
        return
    end
end
if options.beta == "postprocessed"; u_star = beta; end
warning("solve_navier_stokes:noConvergence", ...
    "Picard did not reach %.1e in %d iterations (last %.3e)", ...
    options.tolerance, options.max_iterations, history(end));

    function p = with_beta(p, b)
        if declarative
            p.Operator = options.operator_of(b);
        else
            p.Data.convective_velocity = b;
        end
    end
end

function [beta, u_dof] = next_beta(solution, choice)
u_dof = solution.Unknowns.field(2).dof;      % [u1; u2] per element
d2 = size(u_dof, 1)/2;
switch choice
    case "velocity"
        beta = {u_dof(1:d2, :); u_dof(d2+1:end, :)};
    case "postprocessed"
        beta = divergence_free_velocity(solution);
end
end
