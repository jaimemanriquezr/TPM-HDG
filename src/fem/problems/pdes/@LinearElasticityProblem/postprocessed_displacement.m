function U_star = postprocessed_displacement(obj, mesh, k, s)
% Element-wise P_k+1 postprocessing of the HDG displacement.
%
%   U_star = postprocessed_displacement(obj, mesh, k, s)
%
% s is a solid state (the layout of state.solid) of one solution or of a
% linear combination of solutions (e.g. a Newmark predictor):
% s.displacement (2*d2 x N, Dubiner of degree k) with either s.stress and
% s.rotation (PEERS stress dofs; rotation dofs, a matrix or a cell of
% component blocks), or s.gradient (PEERS dofs of a displacement
% gradient). U_star (2*d2' x N, Dubiner of degree k + 1) solves on each
% element, per component c,
%   (grad u*_c, grad w)_K = (q_c, grad w)_K   for all w in P_k+1,
%   mean_K u*_c = mean_K u_c,
% with q the reconstructed displacement gradient: c1 sigma - c2 tr(sigma) I
% + gamma from the stress and the rotation, or s.gradient. With q of
% order k + 1, u* is of order k + 2 and its gradient of order k + 1.
arguments
    obj (1,1) LinearElasticityProblem
    mesh
    k (1,1) double
    s (1,1) struct
end
if isfield(s, "gradient")
    U_star = postprocessed_field(obj, mesh, k, s.displacement, s.gradient, []);
    return
end
gamma = s.rotation;
if iscell(gamma), gamma = cat(1, gamma{:}); end
U_star = postprocessed_field(obj, mesh, k, s.displacement, s.stress, gamma);
end
