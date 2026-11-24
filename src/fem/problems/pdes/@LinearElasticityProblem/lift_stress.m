function [sigma, gamma] = lift_stress(obj, mesh, k, u, u_hat)
% Stress and rotation of a given displacement and trace, element by element.
%
%   [sigma, gamma] = lift_stress(obj, mesh, k, u, u_hat)
%
% solves the constitutive and symmetry rows of the local problem with
% (u, û) given:
%   c1 (σ, τ) - c2 (tr σ I, τ) + (γ, τ) = <û, τ n>_∂K - (u, ∇·τ)   ∀τ
%                                (σ, η) = 0                          ∀η
% i.e. rows 1 and 3 of get_matrices' A1 [σ; u; γ] = A2 û + Af with the
% u column moved to the right. For a solved step (u, û) this returns the
% solved (σ, γ). Linear in (u, û).
%   u      volume dofs, 2*d2 x N (x block over y block)
%   u_hat  skeleton dofs, 2*d1 x F
%   sigma  PEERS stress dofs, d2peers x N (the layout of
%          solution.Unknowns.field(1).dof)
%   gamma  rotation dofs, d2 x N
arguments
    obj   (1,1) LinearElasticityProblem
    mesh        SimplicialGrid
    k     (1,1) double
    u           double
    u_hat       double
end
[A1, A2, Af] = static_local_operators(obj, mesh, k);
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = size(A1, 1) - 3*d2;
ix_s = 1:d2peers;
ix_u = d2peers + (1:2*d2);
ix_g = d2peers + 2*d2 + (1:d2);
ix = [ix_s, ix_g];

u_K = element_traces(mesh, u_hat, d1);
rhs = pagemtimes(A2(ix, :, :), u_K) + Af(ix, :, :) ...
    - pagemtimes(A1(ix, ix_u, :), reshape(u, 2*d2, 1, N_ELEMENTS));
x = reshape(pagemldivide(A1(ix, ix, :), rhs), [], N_ELEMENTS);
sigma = x(1:d2peers, :);
gamma = x(d2peers + 1:end, :);
end
