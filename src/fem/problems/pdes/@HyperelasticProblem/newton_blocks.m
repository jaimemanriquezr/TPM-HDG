function [lhs, rhs] = newton_blocks(pde, mesh, k)
%NEWTON_BLOCKS  The Newton step of get_matrices as declarative blocks.
%
%   [lhs, rhs] = newton_blocks(pde, mesh, k)
%
%   The terms that turn the declarative linear-elasticity operator
%   (LinearElasticityMatrices: local ids 1, 2, 3 = σ, u, γ, global id 4
%   = û, rows A1 X − A2 û = Af) into the Newton step of this class at the
%   iterate Data.iterate (see get_matrices): in the constitutive (test
%   id 1) and weak-symmetry (test id 3) rows,
%     lhs  −Dm R_u on u (trial 2) and −Dm R_û on û (trial 4),
%     rhs  m⁰ − Dm L⁰ (trial −1),
%   with m, Dm the moments of the stress remainder at the discrete
%   gradient L⁰ of the iterate (nonlinear_terms, gradient_operators),
%   the Kelvin–Voigt stress included when Data.solid_viscosity > 0 (its
%   history read from Data.state, which the caller sets to the state of
%   the step).
%   lhs and rhs are IndexedMatrix block lists; add both to the solid
%   operator. The hybrid blocks carry the face in dimension 7, the
%   x block of the face's 2(k+1) trace dofs over the y block.
arguments
    pde  (1,1) HyperelasticProblem
    mesh       SimplicialGrid
    k    (1,1) double
end
assert(pde.Data.stabilization == "constant", "HyperelasticProblem:newtonBlocks", ...
    "newton_blocks supports stabilization = ""constant"" only.");
iterate = pde.Data.iterate;
operators = linear_operators(pde, mesh, k);
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
d2peers = size(operators.R_u, 1);

u0 = reshape(iterate.displacement, 2*d2, 1, N_ELEMENTS);
u0_K = element_traces(mesh, iterate.displacement_skeleton, d1);
l0 = pagemtimes(operators.R_u, u0) + pagemtimes(operators.R_hat, u0_K);
[m, Dm] = nonlinear_terms(pde, mesh, k, l0, operators.viscous);
B_u = -pagemtimes(Dm, operators.R_u);           % (d2peers + d2) x 2 d2 x N
B_hat = -pagemtimes(Dm, operators.R_hat);       % (d2peers + d2) x 6 d1 x N
f = m - pagemtimes(Dm, l0);                     % (d2peers + d2) x 1 x N

volume = @(A) reshape(A, size(A, 1), size(A, 2), 1, 1, N_ELEMENTS);
% columns [x faces 1..3, y faces 1..3] -> 2 d1 columns per face, face in dim 7
hybrid = @(A) permute(reshape(permute(reshape(A, size(A, 1), d1, 3, 2, N_ELEMENTS), ...
    [1 2 4 3 5]), size(A, 1), 2*d1, 3, N_ELEMENTS), [1 2 5 6 4 7 3]);
rows = {1:d2peers, d2peers + (1:d2)};
lhs = IndexedArray(volume(B_u(rows{1}, :, :)), [2, 1]) ...
    + IndexedArray(hybrid(B_hat(rows{1}, :, :)), [4, 1]) ...
    + IndexedArray(volume(B_u(rows{2}, :, :)), [2, 3]) ...
    + IndexedArray(hybrid(B_hat(rows{2}, :, :)), [4, 3]);
rhs = IndexedArray(volume(f(rows{1}, :, :)), [-1, 1]) ...
    + IndexedArray(volume(f(rows{2}, :, :)), [-1, 3]);
end
