function [cdt, history] = fluid_step_history(time_scheme, time_step, fluid_state)
%FLUID_STEP_HISTORY  Effective time step and history combinations of a fluid step.
%
%   [cdt, history] = fluid_step_history(time_scheme, time_step, fluid_state)
%
%   The time derivative of the fluid momentum equation is replaced by
%       ∂_t u ≈ (u - w*)/cdt,
%   with
%       "euler":  cdt = Δt,       w* = uⁿ,
%       "bdf2":   cdt = 2Δt/3,    w* = (4uⁿ - uⁿ⁻¹)/3,
%   where uⁿ is the field velocity of fluid_state (state.fluid of
%   StokesBiotCoupling.initial_state / advance_state) and uⁿ⁻¹ its field
%   velocity_previous. Under "bdf2" a state whose velocity_previous is
%   missing or [] (no level n-1, e.g. the initial state) gives the Euler
%   step; the scalar 0 there is a zero level n-1. This is the rule of
%   biot_step_history, so the two media of a coupled step agree when both
%   states are advanced together.
%
%   history has the fields
%     velocity             w*, the volume dof array (2*d2 x n_el) of the
%                          inertia term; a missing field, [] or the scalar
%                          0 is zero, a combination of zeros is the scalar
%                          0, and an empty fluid_state is the zero history;
%     convective_velocity  the extrapolated convective velocity from the
%                          fields u_star (u*ⁿ) and u_star_previous
%                          (u*ⁿ⁻¹), each a 2 x 1 cell {u1*; u2*} of dof
%                          arrays as returned by divergence_free_velocity:
%                          u*ⁿ under "euler" or when u_star_previous is
%                          missing or empty, 2u*ⁿ - u*ⁿ⁻¹ under "bdf2"
%                          otherwise, and {} when the state has no u_star.
%                          It does not depend on the sentinel of velocity.
%
%   This is the only place where the scheme enters the fluid step:
%   IncompressibleElasticityMatrices receives cdt as its time_step (theta
%   = 1) and history.velocity as its state velocity.
%
%See also: biot_step_history, divergence_free_velocity,
%          StokesBiotCoupling.interface_problem
arguments
    time_scheme (1,1) string {mustBeMember(time_scheme, ["euler", "bdf2"])}
    time_step   (1,1) double {mustBePositive, mustBeFinite}
    fluid_state       struct = struct.empty
end
% the level n-1 exists unless velocity_previous is missing or [] (the
% scalar 0 is a zero level n-1)
is_bdf2 = time_scheme == "bdf2" && ~isempty(fluid_state) ...
    && isfield(fluid_state, "velocity_previous") ...
    && ~isempty(fluid_state.velocity_previous);
if is_bdf2
    cdt = 2*time_step/3;
else
    cdt = time_step;
end
history = struct();
current = field_of(fluid_state, "velocity");
if is_bdf2
    previous = field_of(fluid_state, "velocity_previous");
    if isequal(current, 0) && isequal(previous, 0)
        history.velocity = 0;
    else
        history.velocity = (4*current - previous)/3;
    end
else
    history.velocity = current;
end

% the predictor of the convective velocity, from its own levels
history.convective_velocity = {};
if ~isempty(fluid_state) && isfield(fluid_state, "u_star") && ~isempty(fluid_state.u_star)
    u_star = fluid_state.u_star;
    if time_scheme == "bdf2" && isfield(fluid_state, "u_star_previous") ...
            && ~isempty(fluid_state.u_star_previous)
        previous = fluid_state.u_star_previous;
        u_star = {2*u_star{1} - previous{1}; 2*u_star{2} - previous{2}};
    end
    history.convective_velocity = reshape(u_star, 2, 1);
end
end

function value = field_of(s, name)
% the field of the state, the scalar 0 when absent, empty or zero
value = 0;
if ~isempty(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end
