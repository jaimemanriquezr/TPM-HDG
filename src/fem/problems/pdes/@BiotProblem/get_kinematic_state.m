function state = get_kinematic_state(obj, solution, state)
% Advances the time-stepping state by the solved step.
%
%   state = get_kinematic_state(obj, solution, state)
%
% takes the state the step was solved from (obj.Data.state; see
% initial_state) and the solution of that step, and returns the state
% advanced: state.time += Data.time_step, state.step += 1, and in
% state.biot every history field moved to its "_previous" field and
% replaced by the solved one:
%   stress                ← solution.Unknowns.field(1)
%   rotation              ← field(3)
%   pressure              ← field(5)
%   displacement_skeleton ← hybrid(1) rows 1:2*d1   (η̂_x; η̂_y)
%   displacement          ← field(2)   (volume η)
%   velocity              ← (field(2) − η*)/cdt, with (cdt, η*) =
%                           biot_step_history(time_scheme, time_step,
%                           state.biot) of the input state
%   pressure_skeleton     ← hybrid(1) rows 2*d1+1:3*d1   (p̂; not a history,
%                           no "_previous")
% A current field that is missing, [] or the scalar 0 moves to
% "_previous" as a zeros array (a zero level n-1, not the "no level n-1"
% sentinel of biot_step_history). An empty input state is the zero
% history at t = Data.time, step 0: the first step from it is a step from
% rest. The fluid and diagnostics fields are kept.
%
% The same call advances the Biot side of the coupled problem
% (StokesBiotCoupling.solve), with the side-2 solution.
%
%See also: BiotProblem/initial_state, biot_step_history
arguments
    obj      (1,1) BiotProblem
    solution (1,1) FESolution
    state          struct = struct.empty
end
if isempty(state)
    state = struct("time", obj.Data.time, "step", 0, "biot", struct(), ...
        "fluid", struct.empty, "diagnostics", struct.empty);
end
b = struct();
if isfield(state, "biot") && ~isempty(state.biot)
    b = state.biot;
end

unknowns = solution.Unknowns;
sigma1 = unknowns.field(1).dof;
gamma1 = unknowns.field(3).dof;
if iscell(gamma1), gamma1 = cat(1, gamma1{:}); end
p1     = unknowns.field(5).dof;
trace1 = unknowns.hybrid(1).dof;
d1 = size(trace1, 1)/3;
eta1_hat = trace1(1:2*d1, :);
p1_hat   = trace1(2*d1+1:end, :);

% the volume velocity of the step, with the (cdt, eta*) it was solved with
[cdt, h] = biot_step_history(obj.Data.time_scheme, obj.Data.time_step, b);
eta1 = unknowns.field(2).dof;
v1 = (eta1 - h.displacement)/cdt;

new = struct("displacement_skeleton", eta1_hat, "stress", sigma1, ...
    "rotation", gamma1, "pressure", p1, "displacement", eta1, "velocity", v1);
for name = string(fieldnames(new)).'
    b.(name + "_previous") = default_zeros(field_of(b, name), new.(name));
    b.(name) = new.(name);
end
b.pressure_skeleton = p1_hat;

state.biot = b;
state.time = state.time + obj.Data.time_step;
state.step = state.step + 1;
end

function value = field_of(s, name)
% a state.biot field, [] when absent
value = [];
if isfield(s, name), value = s.(name); end
end

function A = default_zeros(A, template)
% Empty or scalar-0 fields expand to a zero array shaped like the template.
if isempty(A) || isequal(A, 0)
    A = zeros(size(template));
end
end
