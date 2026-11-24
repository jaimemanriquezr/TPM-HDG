function [local_solutions, trace_solutions] = get_local_solutions(obj, x_hat, local_operators, mesh)
%-------------------------------------------------------------------------
%Description: Recovers the local volume solutions (u,p) per element (flux u
%             and pressure p) by solving the element-local problem
%               A1*[u;p] = A2*x_hat_K + Af
%             where x_hat_K are the hybrid DOFs on the boundary of element K.
%-------------------------------------------------------------------------
arguments (Input)
    obj             (1,1)   FEProblem      % Darcy problem object
    x_hat           (:,1)   double          % global hybrid DOF vector (n_dof x 1)
    local_operators (3,1)   cell            % {A1, A2, Af} from get_matrices
    mesh            SimplicialGrid          % triangular mesh
end
arguments (Output)
    local_solutions (:,1)   cell    % {{u1_h; u2_h}, p_h}: volume flux and pressure per element
    trace_solutions (:,1)   cell    % {x_hat}: pressure trace reshaped by face (k+1 x Nfaces)
end

d1 = size(local_operators{2},2)/3;
d2 = size(local_operators{2},1)/3;
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);

faces = mesh.boundaries(:);
% index u_hat D.O.F. by element
x_hat = reshape(x_hat, d1, N_FACES);
u_aux = reshape(x_hat(:, faces), [3*d1, N_ELEMENTS]);
[A1, A2, Af] = deal(local_operators{:});
% page-wise over elements
rhs = pagemtimes(A2, reshape(u_aux, [], 1, N_ELEMENTS)) + Af;
sol = reshape(pagemldivide(A1, rhs), [], N_ELEMENTS);

u1_h = sol(0 + (1:d2),:);
u2_h = sol(d2 + (1:d2),:);
p_h = sol(2*d2 + (1:d2),:);

local_solutions = {{u1_h; u2_h}, p_h};
trace_solutions = {x_hat};
end
