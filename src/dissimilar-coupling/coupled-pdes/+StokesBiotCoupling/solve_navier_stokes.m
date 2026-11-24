function [fluid_solution, biot_solution, state, info] = solve_navier_stokes(fluid_problem, ...
    biot_problem, fluid_mesh, biot_mesh, interface, sources, options)
%SOLVE_NAVIER_STOKES  One time step of the TPM-coupled Navier-Stokes-Biot problem.
%
%   [fs, bs, state, info] = StokesBiotCoupling.solve_navier_stokes( ...
%       fluid_problem, biot_problem, fluid_mesh, biot_mesh, interface, ...
%       sources, state = state, friction_bjs = g, mode = "iterated", ...
%       tolerance = 1e-6, max_iterations = 10, Cache = cache)
%
%   The step of StokesBiotCoupling.solve with the fluid's convective term
%   rho_c (beta . grad) u (rho_c = rho_F when the fluid density rho_F > 0,
%   else 1). The convective velocity beta is always the divergence-free
%   postprocessed velocity u* (divergence_free_velocity) of a fluid
%   solution, never the raw u_h; the loop owns beta, so a
%   convective_velocity in fluid_problem.Data is ignored, and
%   convective_stabilization is set to false on a copy of the problem:
%   the element and the transmission rows both use S = mu tau with tau
%   read from fluid_mesh.tau, which the caller sets (e.g. per face from a
%   bound of |u_F|; the mesh is not modified here).
%
%   Predictor: beta0 = fluid_step_history(time_scheme, time_step,
%   state.fluid).convective_velocity, i.e. u*^n (euler, or bdf2 without
%   state.fluid.u_star_previous) or 2 u*^n - u*^{n-1} (bdf2). Without
%   state.fluid.u_star (e.g. the state of StokesBiotCoupling
%   .initial_state) the step first solves the coupled problem without
%   convection (beta = 0, the same step) and takes beta0 = u* of it.
%
%   mode:
%     "semi-implicit"  one coupled solve with beta0;
%     "iterated"       Picard: solve with beta, beta <- u* of the
%                      solution, until the relative increment of the
%                      fluid velocity dofs ||u_j - u_{j-1}|| / ||u_j|| <
%                      tolerance, or max_iterations solves (warning
%                      StokesBiotCoupling:noConvergence). The first solve
%                      with a predictor has no previous iterate (its
%                      increment is NaN), so at least two solves are made;
%                      after a Stokes pre-solve the first increment is
%                      measured against it.
%   Every iterate solves the same step from the same input state (the
%   history of level n). The returned state is advanced once, by
%   StokesBiotCoupling.advance_state with the final solutions, and then
%     state.fluid.u_star_previous <- the input state.fluid.u_star ([]
%                                    when it has none, so that a bdf2 step
%                                    from it predicts with u* alone)
%     state.fluid.u_star          <- u* of the final fluid solution.
%   Start the state with u_star_previous missing or [] (as
%   StokesBiotCoupling.initial_state does): the first predictor is then
%   u*^0 (or the Stokes pre-solve), an Euler predictor under both schemes.
%
%   info:
%     iterations        number of coupled solves with convection
%     stokes_presolve   true when a coupled solve without convection
%                       preceded them (no predictor in the state)
%     increments        relative increments of the fluid velocity dofs,
%                       one per convective solve (NaN for the first one
%                       when there is no previous iterate)
%     max_bn_interface  max |beta . n| over the Sigma faces of the fluid
%                       mesh (fluid_mesh.intfaces)
%     max_bn            max |beta . n| over all faces of the fluid mesh
%     s_margin          min over element faces and face quadrature points
%                       of S - rho_c (beta . n_K)/2, S = mu fluid_mesh.tau
%                       of that element face, n_K the element's outward
%                       normal, at the points of the element's convective
%                       face term (warning StokesBiotCoupling:stabilization
%                       when <= 0)
%     (all three for the beta of the last solve)
%     u_star            {u1*; u2*} of the final fluid solution
%     divergence        max_K ||div u*||_{L2(K)} of it
%     Cache             the iterate-independent fluid element blocks
%                       (see StokesBiotCoupling.interface_problem), to
%                       pass back as Cache in the next step on the same
%                       fluid mesh and tau
%
%   The other options are those of StokesBiotCoupling.solve.
%
%See also: StokesBiotCoupling.solve, StokesBiotCoupling.interface_problem,
%          StokesBiotCoupling.advance_state, fluid_step_history,
%          divergence_free_velocity, solve_navier_stokes
arguments
    fluid_problem IncompressibleElasticityProblem
    biot_problem  BiotProblem
    fluid_mesh    SimplicialGrid
    biot_mesh     SimplicialGrid
    interface     AbstractCurve
    sources       struct = struct.empty
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.friction_bjs (1,1) double {mustBeNonnegative} = 0.0
    options.pore_pressure_mean (1,1) double = 0.0
    options.state struct = struct.empty
    options.mode (1,1) string {mustBeMember(options.mode, ["iterated", "semi-implicit"])} = "iterated"
    options.tolerance (1,1) double {mustBePositive} = 1e-6
    options.max_iterations (1,1) double {mustBeInteger, mustBePositive} = 10
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.vertex_tol        (1,1) double = 1e-8
    options.vertex_tol_relative double = []
    options.Cache struct = struct.empty
end
k = options.Order;
fd = fluid_problem.Data;
bd = biot_problem.Data;
fluid_state = struct.empty;
if ~isempty(options.state) && isfield(options.state, "fluid")
    fluid_state = options.state.fluid;
end
[~, history] = fluid_step_history(bd.time_scheme, bd.time_step, fluid_state);
rho_c = 1;
if fd.density > 0, rho_c = fd.density; end

fluid_problem.Data.convective_stabilization = false;
cache = options.Cache;
info = struct("iterations", 0, "stokes_presolve", false, "increments", zeros(1, 0));

% the predictor, or a solve without convection that provides it
beta = history.convective_velocity;
u_previous = [];
if isempty(beta)
    [fluid_solution, biot_solution] = coupled_solve({});
    beta = divergence_free_velocity(fluid_solution);
    u_previous = fluid_solution.Unknowns.field(2).dof;
    info.stokes_presolve = true;
end

n_solves = options.max_iterations;
if options.mode == "semi-implicit", n_solves = 1; end
converged = false;
for j = 1:n_solves
    beta_used = beta;
    [fluid_solution, biot_solution] = coupled_solve(beta_used);
    [beta, divergence] = divergence_free_velocity(fluid_solution);
    u = fluid_solution.Unknowns.field(2).dof;
    increment = NaN;
    if ~isempty(u_previous)
        increment = norm(u - u_previous, "fro")/max(norm(u, "fro"), realmin);
    end
    info.iterations = j;
    info.increments(j) = increment;
    u_previous = u;
    if options.mode == "iterated" && increment < options.tolerance
        converged = true;
        break
    end
end
if options.mode == "iterated" && ~converged
    warning("StokesBiotCoupling:noConvergence", ...
        "Picard did not reach %.1e in %d iterations (last increment %.3e).", ...
        options.tolerance, options.max_iterations, info.increments(end));
end

% the face check of the beta of the last solve
[info.max_bn_interface, info.max_bn, info.s_margin] = convective_margin(fluid_mesh, k, ...
    beta_used, fd.viscosity, rho_c);
if ~(info.s_margin > 0)
    warning("StokesBiotCoupling:stabilization", ...
        "S - rho_c beta . n_K/2 = %.3e <= 0 on some element face: raise fluid_mesh.tau " + ...
        "(max |beta . n| = %.3e).", info.s_margin, info.max_bn);
end
info.u_star = beta;
info.divergence = divergence;
info.Cache = cache;

% advance the state once, from the input state
u_star_old = [];
if ~isempty(fluid_state) && isfield(fluid_state, "u_star")
    u_star_old = fluid_state.u_star;
end
state = StokesBiotCoupling.advance_state(fluid_problem, biot_problem, ...
    fluid_solution, biot_solution, options.state);
state.fluid.u_star_previous = u_star_old;
state.fluid.u_star = info.u_star;

    function [fs, bs] = coupled_solve(b)
        % the coupled step from the input state with convective velocity b
        problem = fluid_problem;
        problem.Data.convective_velocity = b;
        [coupled, ~, ~, cache] = StokesBiotCoupling.interface_problem(problem, ...
            biot_problem, fluid_mesh, biot_mesh, interface, sources, Order = k, ...
            friction_bjs = options.friction_bjs, state = options.state, ...
            pore_pressure_mean = options.pore_pressure_mean, ...
            n_quadpoints = options.n_quadpoints, n_quadpoints_path = options.n_quadpoints_path, ...
            vertex_tol = options.vertex_tol, vertex_tol_relative = options.vertex_tol_relative, ...
            Cache = cache);
        [fs, bs] = solve(coupled, k);
    end
end

function [max_bn_interface, max_bn, s_margin] = convective_margin(mesh, k, beta, mu, rho_c)
% beta . n_K at the face quadrature points of convective_matrix (beta
% dofs of degree kb, Gauss rule of max(ceil((kb + 2k)/2), 1) points),
% n_K the outward unit normal of the element face
kb = find(arrayfun(@(d) nchoosek(d + 2, 2), 0:k+2) == size(beta{1}, 1)) - 1;
kq = max(ceil((kb + 2*k)/2), 1);
s = get_1d_reference_quadrature(kq);
N_NODES = numel(s);
quadrature_2d = [[0*s(:), 1 - s(:)];
                 [1 - s(:), 0*s(:)];
                 [1 - s(:), s(:)]];
pb = dubiner_basis(kb, quadrature_2d(:,1), quadrature_2d(:,2));
bx = pb*beta{1};  by = pb*beta{2};                     % 3nq x N
area = mesh.area(mesh.boundaries);
oriented_normals = mesh.normals .* kron(mesh.orientation, [1;1]);
on_sigma = ismember(mesh.boundaries, mesh.intfaces);
[max_bn_interface, max_bn, s_margin] = deal(0, 0, Inf);
for e = 1:3
    nodes_e = (e-1)*N_NODES + (1:N_NODES);
    bn = (oriented_normals(2*(e-1) + 1, :) .* bx(nodes_e, :) ...
        + oriented_normals(2*(e-1) + 2, :) .* by(nodes_e, :)) ./ area(e, :);   % nq x N
    max_bn = max(max_bn, max(abs(bn), [], "all"));
    if any(on_sigma(e, :))
        max_bn_interface = max(max_bn_interface, max(abs(bn(:, on_sigma(e, :))), [], "all"));
    end
    s_margin = min(s_margin, min(mu*mesh.tau(e, :) - rho_c*max(bn, [], 1)/2));
end
end
