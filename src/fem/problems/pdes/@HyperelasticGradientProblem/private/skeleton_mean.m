function a_hat = skeleton_mean(mesh, k, a)
% face L² projection of the mean of the one-sided traces of the volume
% dofs a (2*d2 x N): on each face F, Σ_K <a_K, L>_F / n_F against the face
% mass, n_F = the number of elements sharing F
d1 = k + 1;
d2 = size(a, 1)/2;
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
% <P_j, L_i>_e and <L_j, L_i>_e per element and local face e (unit weight)
[~, ~, ~, LP, LL] = boundary_matrix(add_stabilization(mesh, 1), k);
incidence = sparse(1:3*N_ELEMENTS, mesh.boundaries(:), 1, 3*N_ELEMENTS, N_FACES);
face_mass = zeros(d1, d1, 3, N_ELEMENTS);
for e = 1:3
    block_e = (e-1)*d1 + (1:d1);
    face_mass(:, :, e, :) = reshape(LL(block_e, block_e, :), d1, d1, 1, N_ELEMENTS);
end
face_mass = reshape(reshape(face_mass, d1*d1, 3*N_ELEMENTS)*incidence, d1, d1, N_FACES);
a_hat = zeros(2*d1, N_FACES);
for c = 1:2
    a_c = reshape(a((c-1)*d2 + (1:d2), :), d2, 1, N_ELEMENTS);
    moments = reshape(pagemtimes(LP, a_c), d1, 3*N_ELEMENTS)*incidence;
    a_hat((c-1)*d1 + (1:d1), :) = reshape( ...
        pagemldivide(face_mass, reshape(moments, d1, 1, N_FACES)), d1, N_FACES);
end
end
