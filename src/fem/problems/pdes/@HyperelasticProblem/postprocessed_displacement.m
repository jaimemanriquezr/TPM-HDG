function U_star = postprocessed_displacement(pde, mesh, k, s)
% Element-wise P_k+1 postprocessing of the hyperelastic HDG displacement.
%
%   U_star = postprocessed_displacement(pde, mesh, k, s)
%
% s is a solid state (the layout of state.solid or Data.iterate) of one
% solution or of a linear combination of solutions: s.displacement
% (2*d2 x N) and s.displacement_skeleton (2*(k+1) x F). The reconstructed
% gradient is the discrete gradient L_h of (u_h, û_h) in S_h,
%   (L_h, τ)_K = −(u_h, ∇·τ)_K + <û_h, τ n>_∂K   for every τ in S_h,
% the field at which the stress remainder is evaluated. With N ≡ 0, L_h
% is the S_h projection of the linear class's c1 σ − c2 tr(σ) I + γ, and
% since grad P_k+1 lies in S_h the result is that of
% LinearElasticityProblem.postprocessed_displacement. U_star
% (Dubiner of degree k + 1) solves, per element and component c,
%   (grad u*_c, grad w)_K = (L_h,c, grad w)_K,   mean_K u*_c = mean_K u_c.
% L_h is linear in (u_h, û_h), so a combination of states gives the same
% combination of reconstructions.
arguments
    pde  (1,1) HyperelasticProblem
    mesh       SimplicialGrid
    k    (1,1) double
    s    (1,1) struct
end
[R_u, R_hat] = HyperelasticProblem.discrete_gradient_operators(mesh, k);
N_ELEMENTS = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
u = reshape(s.displacement, 2*d2, 1, N_ELEMENTS);
u_K = element_traces(mesh, s.displacement_skeleton, k + 1);
l = pagemtimes(R_u, u) + pagemtimes(R_hat, u_K);
U_star = postprocessed_displacement@LinearElasticityProblem(pde, mesh, k, ...
    struct("displacement", s.displacement, "gradient", reshape(l, [], N_ELEMENTS)));
end
