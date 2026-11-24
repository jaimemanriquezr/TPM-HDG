function state = advance_state(fluid_problem, biot_problem, fluid_solution, biot_solution, state)
%ADVANCE_STATE  The coupled Stokes-Biot state advanced by a solved step.
%
%   state = StokesBiotCoupling.advance_state(fluid_problem, biot_problem, ...
%       fluid_solution, biot_solution, state)
%
%   takes the state the step was solved from (see
%   StokesBiotCoupling.initial_state; an empty state is the zero history
%   at biot_problem.Data.time) and the two solutions of that step, and
%   returns it advanced:
%     - the Biot part, state.time and state.step by
%       get_kinematic_state(biot_problem, biot_solution, state)
%       (state.time += Data.time_step, state.step += 1, once);
%     - state.fluid:
%         velocity_previous  ← velocity (zeros when it is missing, [] or
%                              the scalar 0: a zero level n-1)
%         velocity           ← u¹              (unknowns.field(2))
%         velocity_skeleton  ← û¹              (unknowns.hybrid(1))
%         stress             ← σ¹              (unknowns.field(1), PEERS dofs)
%         pressure           ← -tr(σ¹)/2       (Dubiner dofs, d2 x n_el)
%       every other field of state.fluid (e.g. u_star, u_star_previous of
%       a Navier-Stokes step) is kept untouched.
%   The fluid is advanced whatever its density (with density 0 the step
%   does not read state.fluid). No acceleration is stored: the Euler and
%   BDF2 steps of StokesBiotCoupling read only velocity and
%   velocity_previous (fluid_step_history).
%
%See also: StokesBiotCoupling.solve, StokesBiotCoupling.initial_state,
%          fluid_step_history, BiotProblem/get_kinematic_state
arguments
    fluid_problem  IncompressibleElasticityProblem %#ok<INUSA>
    biot_problem   BiotProblem
    fluid_solution (1,1) FESolution
    biot_solution  (1,1) FESolution
    state          struct = struct.empty
end
fluid = struct();
if ~isempty(state) && isfield(state, "fluid") && ~isempty(state.fluid)
    fluid = state.fluid;
end
state = get_kinematic_state(biot_problem, biot_solution, state);

unknowns = fluid_solution.Unknowns;
u1 = unknowns.field(2).dof;
sigma1 = unknowns.field(1).dof;
u0 = [];
if isfield(fluid, "velocity"), u0 = fluid.velocity; end
if isempty(u0) || isequal(u0, 0)
    u0 = zeros(size(u1));
end
fluid.velocity_previous = u0;
fluid.velocity          = u1;
fluid.velocity_skeleton = unknowns.hybrid(1).dof;
fluid.stress            = sigma1;
fluid.pressure          = stress_pressure(fluid_solution.Mesh, fluid_solution.Order, sigma1);
state.fluid = fluid;
end

function p = stress_pressure(mesh, k, sigma)
% Dubiner dofs (d2 x N) of p = -tr(σ)/2, the element L² projection of the
% trace of the PEERS stress dofs sigma (d2peers x N), as
% IncompressibleElasticityProblem/get_kinematic_state computes it.
N_ELEMENTS = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
[~, M_tr_PEERS] = PEERS_matrix(mesh, k);
M = mass_matrix(mesh, k);
moments = pagemtimes(M_tr_PEERS(1:d2, :, :), reshape(sigma, [], 1, N_ELEMENTS));
p = -reshape(pagemldivide(M, moments), d2, N_ELEMENTS)/2;
end
