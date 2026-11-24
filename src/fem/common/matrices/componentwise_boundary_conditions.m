function [x_hat, dir_dof, g_hat] = componentwise_boundary_conditions(mesh, specs, k, n_dof_total)
% Assemble Dirichlet values and Neumann loads for a vector trace unknown
% from component-wise boundary_condition specs (one per boundary tag).
%
% DOF layout assumed: [u_hat_x (d1*N_F); u_hat_y (d1*N_F); ...extra],
% with d1 = k + 1 basis functions per face. Trailing PDE-specific DOF
% (e.g. the pressure-mean unknown) are untouched; n_dof_total sizes the
% output vectors.
%
%   x_hat(dir rows) = L2 face projection of the dirichlet data,
%   g_hat(neu rows) = <data, mu>_F traction moments of the neumann data,
% and dir_dof lists every constrained row (spec order, x before y).
arguments
    mesh        SimplicialGrid
    specs       (:,1) struct
    k           (1,1) double
    n_dof_total (1,1) double
end
d1 = k + 1;
N_DOF = d1*size(mesh.faces, 2);     % rows of one scalar trace component
[x_hat, g_hat] = deal(zeros(n_dof_total, 1));
dir_dof = zeros(0, 1);

for s = specs.'
    if ~isKey(mesh.boundary_tags, s.tag)
        error("componentwise_boundary_conditions:unknownTag", ...
            "Boundary tag ""%s"" is not registered on the mesh " + ...
            "(apply_tag it before solving).", s.tag);
    end
    faces = get_face_by_tag(mesh, s.tag);
    if isempty(faces), continue; end
    face_dof = (faces(:) - 1)*d1;
    face_dof = reshape((face_dof + (1:d1)).', [], 1);
    for c = 1:2
        rows = face_dof + (c - 1)*N_DOF;
        if s.type(c) == "dirichlet"
            x_hat(rows) = L2projSkeleton1d(mesh, s.data{c}, k, subset=faces);
            dir_dof = [dir_dof; rows]; %#ok<AGROW>
        else
            g_hat(rows) = test_faces(mesh, s.data{c}, k, subset=faces);
        end
    end
end
end
