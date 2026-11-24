function [fluid_solution, biot_solution, state, global_system] = solve(fluid_problem, ...
    biot_problem, fluid_mesh, biot_mesh, interface, sources, options)
%SOLVE  One time step of the TPM-coupled Stokes-Biot problem.
%
%   [fs, bs, state] = StokesBiotCoupling.solve(fluid_problem, biot_problem, ...
%       fluid_mesh, biot_mesh, interface, sources, state = state, ...
%       friction_bjs = g)
%
%   Builds the coupling with StokesBiotCoupling.interface_problem, solves
%   it, and returns the state advanced by the step
%   (StokesBiotCoupling.advance_state: the Biot part by
%   get_kinematic_state(biot_problem, bs, state), state.fluid from fs,
%   state.time and state.step once; always advanced, an empty input state
%   being the zero history at biot_problem.Data.time).
%   pore_pressure_mean: see StokesBiotCoupling.interface_problem.
%   Start a run from state = StokesBiotCoupling.initial_state(fp, bp, fm,
%   bm, Order = k) (or initial_state(biot_problem, biot_mesh, k) when the
%   fluid is quasi-static, density = 0) and set the sources and boundary
%   data at t_{n+1} before each call:
%     state = StokesBiotCoupling.initial_state(fp, bp, fm, bm, Order = k);
%     for n = 1:n_steps
%         data = mms.at(state.time + dt);    % sources at t_{n+1}
%         ...
%         [fs, bs, state] = StokesBiotCoupling.solve(fp, bp, fm, bm, ...
%             curve, data, state = state, friction_bjs = g);
%     end
%
%See also: StokesBiotCoupling.interface_problem, InterfaceProblem.stokes_biot
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
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.vertex_tol        (1,1) double = 1e-8
    options.vertex_tol_relative double = []
end
coupled = StokesBiotCoupling.interface_problem(fluid_problem, biot_problem, ...
    fluid_mesh, biot_mesh, interface, sources, Order = options.Order, ...
    friction_bjs = options.friction_bjs, state = options.state, ...
    pore_pressure_mean = options.pore_pressure_mean, ...
    n_quadpoints = options.n_quadpoints, n_quadpoints_path = options.n_quadpoints_path, ...
    vertex_tol = options.vertex_tol, vertex_tol_relative = options.vertex_tol_relative);
[fluid_solution, biot_solution, global_system] = solve(coupled, options.Order);
% state = get_kinematic_state(biot_problem, biot_solution, options.state);
state = StokesBiotCoupling.advance_state(fluid_problem, biot_problem, ...
    fluid_solution, biot_solution, options.state);
end
