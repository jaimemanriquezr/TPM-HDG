function [A3, A4] = trace_operator(mesh, k, d2star)
% The trace rows of the linear class, Σ_K <σn − τ(u − û), μ>_∂K:
% A3 [σ; u; γ] + A4 û per element (A4 = <τ û, μ>_∂K), as in
% LinearElasticityProblem.get_matrices.
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
[~, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);
O = zeros(3*d1, d2, N_ELEMENTS);
o = zeros(3*d1, d2star, N_ELEMENTS);
A3 = [nxLP,    O, nyLP,    O, o, -tLP,    O, O;
         O, nxLP,    O, nyLP, o,    O, -tLP, O];
O = zeros(3*d1, 3*d1, N_ELEMENTS);
A4 = [tLL, O;
      O, tLL];
end
