function state = translate_state(obj, mesh, state, shift, k)
% Newmark history of a solid whose reference configuration translates.
%
%   state = translate_state(obj, mesh, state, shift, k)   (k defaults to obj.Order)
%
% The displacement u(x, t) is measured from a reference that translates
% rigidly, x = X + d(t), on a fixed mesh. Its material derivative is
% D_t = d_t + d'(t).grad, and Newmark integrates D_t^2 u along the
% characteristics X = const: the history at x is the previous step's
% field at x - shift, shift = d(t_n) - d(t_n-1). Call this before each
% step with that shift; the step itself and get_kinematic_state are
% unchanged (state.solid.velocity and acceleration are D_t u, D_t^2 u).
%
% Each level (displacement, velocity, acceleration) of state.solid is
% replaced by its translate:
%   volume    the L2 projection of the translate, except for the
%             displacement: u_h + P_k[u*(. - shift) - u*], u* the P_k+1
%             local postprocessing of u_h from the solved stress and
%             rotation (postprocessed_field), whose inter-element jumps
%             are one order smaller than u_h's. Projecting u_h(. - shift)
%             itself leaves a defect of order |shift| h^k per step, which
%             the displacement sums over the steps (an O(h^k) error); the
%             velocity and acceleration enter the Newmark recursions with
%             alternating signs, so their defects do not accumulate, and
%             their stresses are recursions, not solves (a postprocessing
%             from them feeds the history back into itself and is unstable);
%   stress, rotation  the L2 projections of their translates onto their
%             own element spaces (PEERS, Dubiner);
%   skeleton  the face trace extended from x to x - shift along the
%             gradient c1 sigma - c2 tr(sigma) I + gamma of the level's
%             (untranslated) stress and rotation in an element owning the
%             face, as the TPM transfer extends a trace (shifted_skeleton).
% The translates are exact projections: the overlaps K ∩ (K' + shift)
% are clipped exactly, and where x - shift leaves the mesh the
% polynomial of x's element is extrapolated (shifted_projection). No
% translate depends on the translated skeleton: lifting the stress from
% it (lift_stress) feeds the skeleton back into its own extension and
% grows without bound.
% Empty levels stay empty. A zero shift returns the state unchanged.
% Newmark only (the onestep history u^n-2 would need its own shift).
arguments
    obj   (1,1) LinearElasticityProblem
    mesh        SimplicialGrid
    state (1,1) struct
    shift (2,1) double
    k     (1,1) double = obj.Order
end
if all(shift == 0)
    return
end
assert(obj.Data.time_scheme == "newmark", "LinearElasticityProblem:translateScheme", ...
    "translate_state supports time_scheme = ""newmark"" only.");
s = state.solid;
volume   = ["displacement", "velocity", "acceleration"];
skeleton = volume + "_skeleton";
stress   = ["stress", "stress_velocity", "stress_acceleration"];
rotation = ["rotation", "rotation_velocity", "rotation_acceleration"];
levels = find(arrayfun(@(name) isfield(s, name) && ~isempty(s.(name)) ...
    && ~isequal(s.(name), 0.0), volume));
if isempty(levels)
    return
end
d2 = nchoosek(k + 2, 2);
fields = @(names) arrayfun(@(name) s.(name), names, UniformOutput = false);
U_hat = shifted_skeleton(obj, mesh, k, fields(skeleton(levels)), ...
    fields(stress(levels)), fields(rotation(levels)), shift);
U = fields(volume(levels));
U = shifted_projection(mesh, k, cat(1, U{:}), shift, "dubiner");
U = mat2cell(U, 2*d2*ones(1, numel(levels)));
if levels(1) == 1
    U_star = postprocessed_field(obj, mesh, k, s.displacement, s.stress, s.rotation);
    increment = shifted_projection(mesh, k + 1, U_star, shift, "dubiner") - U_star;
    d2p = size(U_star, 1)/2;
    U{1} = s.displacement + increment([1:d2, d2p + (1:d2)], :);
end
Gamma = fields(rotation(levels));
Gamma = shifted_projection(mesh, k, cat(1, Gamma{:}), shift, "dubiner");
Sigma = fields(stress(levels));
Sigma = shifted_projection(mesh, k, cat(1, Sigma{:}), shift, "peers");
n_sigma = size(s.stress, 1);
for i = 1:numel(levels)
    L = levels(i);
    s.(volume(L)) = U{i};
    s.(rotation(L)) = Gamma((i - 1)*d2 + (1:d2), :);
    s.(stress(L)) = Sigma((i - 1)*n_sigma + (1:n_sigma), :);
    s.(skeleton(L)) = U_hat{i};
end
state.solid = s;
end
