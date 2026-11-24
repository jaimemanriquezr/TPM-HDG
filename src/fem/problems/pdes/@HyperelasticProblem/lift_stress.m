function [sigma, gamma] = lift_stress(obj, mesh, k, u, u_hat)
% Stress and rotation of a given displacement and trace, element by element.
%
%   [sigma, gamma] = lift_stress(obj, mesh, k, u, u_hat)
%
% solves the constitutive and symmetry rows of the hyperelastic local
% problem with (u, û) given; they are linear in (σ, γ):
%   (ζ⁻¹σ, τ) + (γ, τ) = <û, τ n>_∂K − (u, ∇·τ) + (ζ⁻¹N(L_h), τ)   ∀τ
%                (σ, η) = (N(L_h), η)                               ∀η
% with L_h the discrete gradient of (u, û) (see get_matrices). For a
% solved step (u, û) this returns the solved (σ, γ), without the
% Kelvin–Voigt part when Data.solid_viscosity > 0. Layouts as in
% LinearElasticityProblem.lift_stress.
arguments
    obj   (1,1) HyperelasticProblem
    mesh        SimplicialGrid
    k     (1,1) double
    u           double
    u_hat       double
end
obj.Order = k;
obj.Data.density = 0;
obj.Data.iterate = struct.empty;
local_operators = get_matrices(obj, mesh, k);
[A1, A2, Af] = deal(local_operators{:});
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = size(A1, 1) - 3*d2;
ix_u = d2peers + (1:2*d2);
ix_c = [1:d2peers, d2peers + 2*d2 + (1:d2)];

u = reshape(u, 2*d2, 1, N_ELEMENTS);
u_K = element_traces(mesh, u_hat, d1);
[R_u, R_hat] = gradient_operators(mesh, k, A1, A2);
m = nonlinear_terms(obj, mesh, k, pagemtimes(R_u, u) + pagemtimes(R_hat, u_K));
rhs = pagemtimes(A2(ix_c, :, :), u_K) + Af(ix_c, :, :) + m ...
    - pagemtimes(A1(ix_c, ix_u, :), u);
x = reshape(pagemldivide(A1(ix_c, ix_c, :), rhs), [], N_ELEMENTS);
sigma = x(1:d2peers, :);
gamma = x(d2peers + 1:end, :);
end
