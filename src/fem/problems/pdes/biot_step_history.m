function [cdt, history] = biot_step_history(time_scheme, time_step, biot_state)
%BIOT_STEP_HISTORY  Effective time step and history combinations of a Biot step.
%
%   [cdt, history] = biot_step_history(time_scheme, time_step, biot_state)
%
%   Every time derivative of a Biot step is replaced by
%       ∂_t w ≈ (w - w*)/cdt,
%   with
%       "euler":  cdt = Δt,       w* = wⁿ,
%       "bdf2":   cdt = 2Δt/3,    w* = (4wⁿ - wⁿ⁻¹)/3,
%   where wⁿ is a field of biot_state (state.biot of
%   BiotProblem.initial_state / get_kinematic_state) and wⁿ⁻¹ its
%   "_previous" field. Under "bdf2" a state whose
%   displacement_skeleton_previous is missing or [] (no level n-1, e.g.
%   the initial state) gives the Euler step; the scalar 0 there is a zero
%   level n-1.
%
%   history has the fields displacement_skeleton, stress, rotation,
%   pressure, displacement and velocity (the last two the volume dof
%   arrays, 2*d2 x n_el, of the inertia term), each the combination w* of
%   the state's field of that name.
%   A missing field, [] or the scalar 0 is zero; a combination of zeros is
%   the scalar 0, and an empty biot_state is the zero history.
%
%   This is the only place where the scheme enters: BiotMatrices and
%   InterfaceProblem.stokes_biot receive cdt and history as numbers and
%   arrays.
%
%See also: BiotProblem, BiotMatrices
arguments
    time_scheme (1,1) string {mustBeMember(time_scheme, ["euler", "bdf2"])}
    time_step   (1,1) double {mustBePositive, mustBeFinite}
    biot_state        struct = struct.empty
end
names = ["displacement_skeleton", "stress", "rotation", "pressure", ...
    "displacement", "velocity"];
% the level n-1 exists unless displacement_skeleton_previous is missing or
% [] (the scalar 0 is a zero level n-1)
is_bdf2 = time_scheme == "bdf2" && ~isempty(biot_state) ...
    && isfield(biot_state, "displacement_skeleton_previous") ...
    && ~isempty(biot_state.displacement_skeleton_previous);
if is_bdf2
    cdt = 2*time_step/3;
else
    cdt = time_step;
end
history = struct();
for name = names
    current = field_of(biot_state, name);
    if is_bdf2
        previous = field_of(biot_state, name + "_previous");
        if isequal(current, 0) && isequal(previous, 0)
            history.(name) = 0;
        else
            history.(name) = (4*current - previous)/3;
        end
    else
        history.(name) = current;
    end
end
end

function value = field_of(s, name)
% the field of the state, the scalar 0 when absent, empty or zero
value = 0;
if ~isempty(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end
