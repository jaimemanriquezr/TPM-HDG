function [a, sigma, gamma, u_hat] = balanced_acceleration(obj, mesh, k, u, options)
% Acceleration of a given volume displacement with flux-balanced traces.
%
%   a = balanced_acceleration(obj, mesh, k, u)
%   [a, sigma, gamma, u_hat] = balanced_acceleration(obj, mesh, k, u, boundary_trace = g_hat)
%
% (σ, γ, û) = balanced_lift(obj, mesh, k, u, ...) (same options), then
% consistent_acceleration of the state (u, σ, γ, û):
%   ρ (a, v)_K = (f, v)_K + (∇·σ, v)_K - <S(u - û), v>_∂K - b₀ (u, v)_K   ∀v
% with f = Data.source_term and b₀ = Data.hooke_modulus. Then M a = f - K u
% with K the statically condensed stiffness acting on u: for a Newmark
% step's solved state (u¹ and its acceleration a¹, with the source and
% boundary data of the new time in obj) it returns a¹. Requires
% density > 0. Returns the 2*d2 x N volume dofs and the lifted fields.
arguments
    obj   (1,1) LinearElasticityProblem
    mesh        SimplicialGrid
    k     (1,1) double
    u           double
    options.boundary_trace double = []
end
[sigma, gamma, u_hat] = obj.balanced_lift(mesh, k, u, boundary_trace = options.boundary_trace);
state.solid = struct("displacement", u, "stress", sigma, "rotation", gamma, ...
    "displacement_skeleton", u_hat);
a = obj.consistent_acceleration(mesh, k, state);
end
