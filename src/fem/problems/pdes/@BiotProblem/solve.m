function [solution, global_system] = solve(obj, mesh, k, options)
%SOLVE  Solve one Biot step on a mesh (declarative path).
%
%   [solution, global_system] = solve(obj, mesh, k)
%   [solution, global_system] = solve(obj, mesh, Order = k)
%
%   = solve(get_declarative_problem(obj, mesh, k), mesh, k). k resolves
%   positional > Order > obj.Order, as in FEProblem.solve. The fields are
%   solution.Unknowns.field(1..5) = (σ, η, γ, q, p) and
%   solution.Unknowns.hybrid(1) the trace, dof 3 d1 x n_faces with rows
%   [η̂_x; η̂_y; p̂].
%
%See also: BiotProblem, get_declarative_problem, FEProblem/solve
arguments
    obj  (1,1) BiotProblem
    mesh       SimplicialGrid
    k          double = []
    options.Order (1,1) double = NaN
end
if isempty(k)
    if ~isnan(options.Order)
        k = options.Order;
    else
        k = obj.Order;
    end
end
fe = get_declarative_problem(obj, mesh, k);
[solution, global_system] = solve(fe, mesh, k);
end
