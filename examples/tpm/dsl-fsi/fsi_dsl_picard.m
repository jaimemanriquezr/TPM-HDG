function [fluid_solution, solid_solution, n_iterations, increments] = fsi_dsl_picard(k, fluid_mesh, solid_mesh, interface, material, pd, interface_sources, options)
% Steady Navier-Stokes/elasticity FSI by Picard iteration on the coupled
% DSL problem of fsi_dsl_coupled_problem: the first solve is Stokes, then
% beta = u* of the previous iterate in the Oseen terms (the divergence-free
% postprocessed velocity in P_{k+1}, divergence_free_velocity), until
% ||u_h^{n+1} - u_h^n|| / ||u_h^{n+1}|| (Frobenius norm of the dofs) is
% below tolerance. pd.f_F must be the Navier-Stokes source (see
% fsi_dsl_convective_source).
%   n_iterations: number of Oseen solves after the Stokes one
%   increments:   relative increment of each Oseen solve (1 x n_iterations)
arguments
    k
    fluid_mesh
    solid_mesh
    interface
    material struct
    pd struct
    interface_sources cell
    options.tolerance (1,1) double = 1e-10
    options.max_iterations (1,1) double = 30
end
mu = material.viscosity;
solve_with = @(mesh, beta) solve(fsi_dsl_coupled_problem(k, mesh, solid_mesh, ...
    interface, material, pd, interface_sources, beta), k);

[fluid_solution, solid_solution] = solve_with(fluid_mesh, []);
beta = postprocessed_velocity(fluid_solution, fluid_mesh, k);
% tau rule: one constant on the whole fluid mesh, fixed over the Oseen
% iterations, tau = 1 + max|u*^Stokes|/(2 mu) (|u| >= |beta.n|), sampled at
% element vertices and centroids.
tau = 1 + max_speed(beta, fluid_mesh)/(2*mu);
fluid_mesh = add_stabilization(fluid_mesh, tau);

u_old = fluid_solution.Unknowns.field(2).dof;
increments = zeros(1, options.max_iterations);
for n = 1:options.max_iterations
    [fluid_solution, solid_solution] = solve_with(fluid_mesh, beta);
    beta = postprocessed_velocity(fluid_solution, fluid_mesh, k);
    u_new = fluid_solution.Unknowns.field(2).dof;
    increments(n) = norm(u_new - u_old, "fro")/norm(u_new, "fro");
    u_old = u_new;
    if increments(n) < options.tolerance
        break
    end
end
n_iterations = n;
increments = increments(1:n);
if increments(end) >= options.tolerance
    warning("fsi_dsl_picard:noConvergence", ...
        "Picard did not reach %.1e in %d iterations (last %.3e)", ...
        options.tolerance, options.max_iterations, increments(end));
end
end

%=========================================================================%
function beta = postprocessed_velocity(solution, mesh, k)
% u* of a fluid solution as an FEField on the vector P_{k+1} space
u_star = divergence_free_velocity(solution);
beta = FEField(FESpace(DubinerBasis(k + 1).to_vector(), mesh), [u_star{1}; u_star{2}]);
end

%=========================================================================%
function s = max_speed(u, mesh)
% max |u| over the reference vertices and centroid of every element
vertices = reference_element(2);
points = [vertices, mean(vertices, 2)];
values = evaluate(u, points, mesh);         % points x 1 x components x 1 x elements
s = max(vecnorm(values, 2, 3), [], "all");
end
