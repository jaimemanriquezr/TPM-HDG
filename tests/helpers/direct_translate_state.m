function s = direct_translate_state(sp, old_mesh, new_mesh, k, old, shift, old_keys, new_keys, s)
%DIRECT_TRANSLATE_STATE  Volume Newmark history on a new carve from the old carve along the characteristics.
%
%   s = direct_translate_state(sp, old_mesh, new_mesh, k, old, shift, ...
%       old_keys, new_keys, s)
%
%   old is the untranslated state.solid on old_mesh (LinearElasticityProblem
%   sp, degree k); s is a state.solid on new_mesh whose skeleton levels are
%   kept (translate_state followed by remap_solid_state) and whose volume
%   levels are replaced. Both meshes are carved from one background grid,
%   their elements identified by the keys (centroids of the unsnapped
%   carve). For every element K of new_mesh and every level, the history is
%   the exact L2 projection onto K's space of the old piecewise field at
%   x - shift (LinearElasticityProblem.translated_field_between: the
%   overlaps with the translated old elements clipped exactly, the part of
%   K that no old element covers taken from the old element of largest
%   overlap, extrapolated):
%     velocity, acceleration, rotation levels   Dubiner P_k;
%     stress levels                              PEERS of K;
%     displacement   on an element present in both carves with the same
%                    vertices, u_K + P_K[u*(. - shift) - u*_K]; elsewhere
%                    P_K[u*(. - shift)]; u* the P_k+1 postprocessing
%                    (postprocessed_displacement) of the old state.
%   Translating on the old carve and then remapping extrapolates twice,
%   on the trailing side (x - shift outside the old carve) and on the
%   leading side (new elements); the direct projection reads the old field
%   where it is defined.
%
%See also: eulerian_solid_march, LinearElasticityProblem.translated_field_between,
%          LinearElasticityProblem.translate_state, remap_solid_state
arguments
    sp (1,1) LinearElasticityProblem
    old_mesh SimplicialGrid
    new_mesh SimplicialGrid
    k (1,1) double
    old (1,1) struct
    shift (2,1) double
    old_keys (:,2) double
    new_keys (:,2) double
    s (1,1) struct
end
between = @(U, degree_from, degree_to, basis, x) sp.translated_field_between(old_mesh, ...
    new_mesh, degree_from, degree_to, U, x, basis = basis);
has = @(name) isfield(old, name) && ~isempty(old.(name)) && ~isequal(old.(name), 0);

for name = ["velocity", "acceleration", "rotation", "rotation_velocity", "rotation_acceleration"]
    if ~has(name), continue, end
    D = old.(name);
    if iscell(D), D = cat(1, D{:}); end
    s.(name) = between(D, k, k, "dubiner", shift);
end
for name = ["stress", "stress_velocity", "stress_acceleration"]
    if ~has(name), continue, end
    s.(name) = between(old.(name), k, k, "peers", shift);
end

% displacement: the translate of u* added to u_K where K is unchanged
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
kept = present;
kept(present) = all(reshape(abs(new_mesh.coordinates(:, new_mesh.elements(:, present)) ...
    - old_mesh.coordinates(:, old_mesh.elements(:, from(present)))) < 1e-14, 6, []), 1).';
U_star = postprocessed_displacement(sp, old_mesh, k, old);
u = between(U_star, k + 1, k, "dubiner", shift);
if any(kept)
    % P_K u*_K: the unshifted projection, K covered by its old self
    own = between(U_star, k + 1, k, "dubiner", [0; 0]);
    u(:, kept) = old.displacement(:, from(kept)) + u(:, kept) - own(:, kept);
end
s.displacement = u;
end
