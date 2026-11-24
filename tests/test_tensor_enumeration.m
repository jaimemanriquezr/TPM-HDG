% test_tensor_enumeration.m
% TensorFEBasis dof enumeration composes the wrapped scalar basis's ids
% with a support-aware component stride (scalar dof count over the mesh),
% giving component-major global numbering while the returned arrays keep
% per-entity columns [comp_1; comp_2; ...].
clc; clear;
mesh = StructuredSimplicialGrid(4, 4);
n_faces = size(mesh.faces, 2);
n_el = size(mesh.elements, 2);
n_fail = 0;

% vector trace (skeleton): stride = d1 * n_faces
VMk = LegendreBasis(1).to_vector(2);
g = get_global_dof(VMk, mesh);
stride = 2 * n_faces;
ok = isequal(size(g), [4, n_faces]) ...
    && isequal(g(:, 1), [1; 2; stride + 1; stride + 2]) ...
    && isequal(g(:, 2), [3; 4; stride + 3; stride + 4]);
if ~ok, fprintf("FAIL trace global enumeration\n"); n_fail = n_fail + 1; end

gs = g(:, 3);
ok = isequal(gs, [5; 6; stride + 5; stride + 6]);   % subset-consistent ids
if ~ok, fprintf("FAIL trace subset consistency\n"); n_fail = n_fail + 1; end

[r, c, r_rhs] = get_assembly_dof(VMk, mesh);
ok = size(r, 3) == n_el && size(r, 4) == 3 && size(r, 5) == 3 ...
    && isequal(size(r, 1:2), [4, 4]) && isequal(size(r_rhs, 2), 1);
f1 = mesh.boundaries(1, 1);   % first side of element 1
ok = ok && isequal(r(:, 1, 1, 1, 1), ...
    [2*(f1-1) + (1:2), stride + 2*(f1-1) + (1:2)].');
if ~ok, fprintf("FAIL trace local enumeration\n"); n_fail = n_fail + 1; end

% vector element support: stride = d * n_elements
VR0 = DubinerBasis(0).to_vector(2);
ge = get_global_dof(VR0, mesh);
ok = isequal(size(ge), [2, n_el]) && isequal(ge(:, 1), [1; n_el + 1]) ...
    && max(ge, [], 'all') == 2*n_el;
if ~ok, fprintf("FAIL element-support tensor enumeration\n"); n_fail = n_fail + 1; end

[re, ~, ~] = get_assembly_dof(VR0, mesh);
ok = ndims(re) == 3 && size(re, 3) == n_el && isequal(re(:, 1, 5), [5; n_el + 5]);
if ~ok, fprintf("FAIL element-support local enumeration\n"); n_fail = n_fail + 1; end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_tensor_enumeration:failed", "%d checks failed", n_fail);
end
