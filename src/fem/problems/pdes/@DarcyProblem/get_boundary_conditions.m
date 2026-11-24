function [x_hat, dir_dof, g_hat, neu_dof] = get_boundary_conditions(obj, mesh, gD, gN)
%-------------------------------------------------------------------------
%Description: Initialises the hybrid DOF vector x_hat and Neumann load
%             vector g_hat, then populates them with boundary data:
%               x_hat(dir_dof) = L2-projection of g_D on Dirichlet faces
%               g_hat(neu_dof) = flux integral of q_N on Neumann faces
%-------------------------------------------------------------------------
arguments (Input)
    obj     DarcyProblem
    mesh    SimplicialGrid      % Mesh
    gD      function_handle     = obj.Data.datum_dirichlet  % Dirichlet datum for pressure
    gN      function_handle     = obj.Data.datum_neumann    % Neumann datum for flux u·n
end
arguments (Output)
    x_hat   (:,1)   double  % hybrid DOF vector with Dirichlet data imposed (n_dof x 1)
    dir_dof (:,1)   double  % Dirichlet-constrained DOF indices ((k+1)*Ndir x 1)
    g_hat   (:,1)   double  % Neumann load vector (n_dof x 1)
    neu_dof (:,1)   double  % Neumann face DOF indices ((k+1)*Nneu x 1)
end
k = obj.Order;
[x_hat, g_hat] = deal(zeros(get_ndof(obj, mesh, k), 1));
d1 = k+1;
%========================================================%
Ndir = size(mesh.dirfaces, 2);
dir_dof = (mesh.dirfaces(:) - 1)*d1 + (1:d1);
dir_dof = reshape(dir_dof', d1*Ndir, 1);
x_hat(dir_dof) = dirichlet_boundary_vector(mesh, gD, k);
%========================================================%

%========================================================%
Nneu = size(mesh.neufaces, 2);
neu_dof = (mesh.neufaces(:) - 1)*d1 + (1:d1);
neu_dof = reshape(neu_dof',d1*Nneu,1);
g_hat(neu_dof) = neumann_boundary_vector(mesh, gN, k);
%========================================================%
end
