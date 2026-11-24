function [x_hat, dir_dof, g_hat] = get_boundary_conditions(obj, mesh, gD, gN)
%-------------------------------------------------------------------------
%Description: Initialises the hybrid DOF vector x_hat and Neumann load
%             vector g_hat, then populates them with boundary data for
%             each vector component of the velocity trace û = (û_x, û_y):
%               x_hat(dir_dof) = L2-projection of g_D on Dirichlet faces
%               g_hat(neu_dof) = flux integral of q_N on Neumann faces
%             DOF layout: [û_x (N_DOF); û_y (N_DOF); ρ (Nelts)]
%-------------------------------------------------------------------------
arguments (Input)
    obj     (1,1)   FEProblem          % Stokes problem object
    mesh            SimplicialGrid      % triangular mesh
    gD      (2,1)   cell                = obj.Data.datum_dirichlet  % Dirichlet data {g_x, g_y}
    gN      (2,1)   cell                = obj.Data.datum_neumann    % Neumann data {q_x, q_y}
end
arguments (Output)
    x_hat   (:,1)   double  % hybrid DOF vector with Dirichlet data imposed (n_dof x 1)
    dir_dof (:,1)   double  % Dirichlet-constrained DOF indices (2*(k+1)*Ndir x 1)
    g_hat   (:,1)   double  % Neumann load vector (n_dof x 1)
end
k = obj.Order;
[x_hat, g_hat] = deal(zeros(get_ndof(obj, mesh, k), 1));
d1 = k + 1;
N_DOF = d1*size(mesh.faces, 2);  % # of DOF belonging to 1 scalar variable
%================================================================%
N_DIRDOF = size(mesh.dirfaces, 2);
dir_dof = (mesh.dirfaces(:) - 1)*d1;
dir_dof = reshape((dir_dof + (1:d1)).', d1*N_DIRDOF, 1);
x_hat(dir_dof)        = dirichlet_boundary_vector(mesh, gD{1}, k);
x_hat(dir_dof + N_DOF) = dirichlet_boundary_vector(mesh, gD{2}, k);
dir_dof = [dir_dof; dir_dof + N_DOF];
%================================================================%

%================================================================%
N_NEUDOF = size(mesh.neufaces, 2);
neuDOF = (mesh.neufaces(:) - 1)*d1;
neuDOF = reshape((neuDOF + (1:d1)).', d1*N_NEUDOF, 1);
g_hat(neuDOF)        = neumann_boundary_vector(mesh, gN{1}, k);
g_hat(neuDOF + N_DOF) = neumann_boundary_vector(mesh, gN{2}, k);
%================================================================%
end
