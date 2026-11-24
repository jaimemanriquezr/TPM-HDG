function [R_u, R_hat] = discrete_gradient_operators(mesh, k)
% The operators of the discrete displacement gradient L_h in S_h (PEERS),
%   (L_h, τ)_K = −(u_h, ∇·τ)_K + <û_h, τ n>_∂K   for every τ in S_h,
% whose dofs are l = R_u u_K + R_hat û_K per element: u_K the 2*d2
% displacement dofs, û_K the 6*(k+1) trace dofs of the element's three
% faces in the order of mesh.boundaries (x block over y block). The
% operators are geometric (no material data); they are those of the
% Newton step (gradient_operators).
%   R_u    d2peers x 2*d2 x N,   R_hat  d2peers x 6*(k+1) x N
arguments
    mesh    SimplicialGrid
    k (1,1) double
end
local_operators = get_matrices(LinearElasticityProblem(k), mesh, k);
[R_u, R_hat] = gradient_operators(mesh, k, local_operators{1}, local_operators{2});
end
