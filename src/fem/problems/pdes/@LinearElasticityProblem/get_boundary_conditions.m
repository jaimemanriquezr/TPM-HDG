function [x_hat, dir_dof, g_hat] = get_boundary_conditions(obj, mesh, gD, gN)
%-------------------------------------------------------------------------
%Description: Initialises the hybrid DOF vector x_hat and Neumann load
%             vector g_hat for the displacement trace û = (û_x, û_y).
%             Boundary conditions are specified component-wise per boundary
%             tag: obj.Data.boundary_conditions holds an array of
%             boundary_condition specs, each fixing û_c (dirichlet) or
%             loading the traction component (σn)_c (neumann) on the faces
%             of one tag. Without specs, full Dirichlet data gD is imposed
%             on the "Dirichlet" tag and full traction data gN on the
%             "Neumann" tag.
%             DOF layout: [û_x (N_DOF); û_y (N_DOF)]
%-------------------------------------------------------------------------
arguments (Input)
    obj     (1,1)   FEProblem          % linear elasticity problem object
    mesh            SimplicialGrid      % triangular mesh
    gD      (2,1)   cell                = obj.Data.datum_dirichlet  % Dirichlet data {g_x, g_y}
    gN      (2,1)   cell                = obj.Data.datum_neumann    % Neumann data {q_x, q_y}
end
arguments (Output)
    x_hat   (:,1)   double  % hybrid DOF vector with Dirichlet data imposed (n_dof x 1)
    dir_dof (:,1)   double  % Dirichlet-constrained DOF indices
    g_hat   (:,1)   double  % Neumann load vector (n_dof x 1)
end
k = obj.Order;
if isfield(obj.Data, "boundary_conditions") && ~isempty(obj.Data.boundary_conditions)
    specs = obj.Data.boundary_conditions;
else
    specs = [ ...
        boundary_condition("Dirichlet", x_data = gD{1}, y_data = gD{2});
        boundary_condition("Neumann", x = "neumann", y = "neumann", ...
            x_data = gN{1}, y_data = gN{2})];
end
[x_hat, dir_dof, g_hat] = componentwise_boundary_conditions( ...
    mesh, specs, k, get_ndof(obj, mesh, k));
end
