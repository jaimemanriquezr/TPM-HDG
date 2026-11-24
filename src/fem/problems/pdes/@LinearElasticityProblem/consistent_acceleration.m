function a = consistent_acceleration(obj, mesh, k, state)
% Volume acceleration satisfying the momentum equation of the local
% problem at the given state, element by element.
%
%   a = consistent_acceleration(obj, mesh, k, state)
%
% solves, with (σ, u, γ, û) read from state.solid (stress, displacement,
% rotation, displacement_skeleton),
%   ρ (a, v)_K = (f, v)_K + (∇·σ, v)_K - <S(u - û), v>_∂K - b₀ (u, v)_K   ∀v
% i.e. ρ [M 0; 0 M] a = V_f - V_s σ - V_U u + V_u û with the static
% momentum-row blocks of get_matrices. This is the identity every Newmark
% step leaves its acceleration in, so a Newmark march started from this a
% satisfies it at every step. f is obj.Data.source_term (set it to f(t)
% of the state), b₀ obj.Data.hooke_modulus; no boundary datum enters.
% When state.solid.stress is empty, (σ, γ) are lifted from (u, û)
% (lift_stress). Returns the 2*d2 x N volume dofs.
arguments
    obj   (1,1) LinearElasticityProblem
    mesh        SimplicialGrid
    k     (1,1) double
    state (1,1) struct
end
assert(obj.Data.density > 0, "LinearElasticityProblem:density", ...
    "consistent_acceleration requires density > 0.");
s = state.solid;
sigma = s.stress;
gamma = s.rotation;
if isempty(sigma)
    [sigma, gamma] = obj.lift_stress(mesh, k, s.displacement, s.displacement_skeleton);
end
[A1, A2, Af] = static_local_operators(obj, mesh, k);
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = size(A1, 1) - 3*d2;
ix_u = d2peers + (1:2*d2);

x = reshape([sigma; s.displacement; gamma], [], 1, N_ELEMENTS);
u_K = element_traces(mesh, s.displacement_skeleton, d1);
rhs = Af(ix_u, :, :) + pagemtimes(A2(ix_u, :, :), u_K) - pagemtimes(A1(ix_u, :, :), x);

M = mass_matrix(mesh, k);
O = zeros(d2, d2, N_ELEMENTS);
mass = obj.Data.density*[M, O; O, M];
a = reshape(pagemldivide(mass, rhs), 2*d2, N_ELEMENTS);
end
