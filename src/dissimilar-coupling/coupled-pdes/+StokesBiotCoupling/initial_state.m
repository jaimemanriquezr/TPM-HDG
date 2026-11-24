function state = initial_state(fluid_problem, biot_problem, fluid_mesh, biot_mesh, options)
%INITIAL_STATE  Time-stepping state of the coupled Stokes-Biot problem at t = Data.time.
%
%   state = StokesBiotCoupling.initial_state(fluid_problem, biot_problem, ...
%       fluid_mesh, biot_mesh, Order = k)
%
%   returns initial_state(biot_problem, biot_mesh, k) (state.time =
%   biot_problem.Data.time, state.step = 0, state.biot with the levels n
%   and no level n-1) with state.fluid
%     velocity           u⁰ = Π u(0) onto the elements (L2proj2D,
%                        2*d2 x n_el) of fluid_problem.Data.initial_velocity
%     velocity_skeleton  û⁰ = P_M u(0) onto the faces (L2projSkeleton1d)
%     velocity_previous  [] (no level n-1: a "bdf2" step from this state
%                        is an Euler step on both media, see
%                        fluid_step_history and biot_step_history)
%     stress, pressure   [] (not read by the step)
%   No acceleration problem is solved: the Euler and BDF2 steps of
%   StokesBiotCoupling need only the velocity. Advance the state with
%   StokesBiotCoupling.solve (or StokesBiotCoupling.advance_state).
%
%See also: StokesBiotCoupling.advance_state, StokesBiotCoupling.solve,
%          BiotProblem/initial_state, IncompressibleElasticityProblem/initial_state
arguments
    fluid_problem IncompressibleElasticityProblem
    biot_problem  BiotProblem
    fluid_mesh    SimplicialGrid
    biot_mesh     SimplicialGrid
    options.Order (1,1) double = fluid_problem.Order
end
k = options.Order;
state = initial_state(biot_problem, biot_mesh, k);
u0 = fluid_problem.Data.initial_velocity;
state.fluid = struct( ...
    "velocity",          L2proj2D(fluid_mesh, u0, k), ...
    "velocity_skeleton", L2projSkeleton1d(fluid_mesh, u0, k), ...
    "velocity_previous", [], ...
    "stress",            [], ...
    "pressure",          []);
end
