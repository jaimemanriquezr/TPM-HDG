function [local_solutions, trace_solutions] = get_local_solutions(obj, x_hat, local_operators, mesh)
%-------------------------------------------------------------------------
%Description: Recovers the local volume solutions (L,u,p) per element by
%             solving the element-local problem
%               A1*[L;u;p;ξ] = A2*[û_K; ρ_K] + Af
%             where û_K and ρ_K are the hybrid DOFs on element K.
%-------------------------------------------------------------------------
arguments (Input)
    obj             (1,1)   FEProblem      % Stokes problem object
    x_hat           (:,1)   double          % global hybrid DOF vector (n_dof x 1)
    local_operators (3,1)   cell            % {A1, A2, Af} from get_matrices
    mesh            SimplicialGrid          % triangular mesh
end
arguments (Output)
    local_solutions (:,1)   cell    % {L_h (2x2 cell), {u1_h;u2_h}, p_h}: volume fields per element
    trace_solutions (:,1)   cell    % {{ux_hat;uy_hat}, rho}: velocity trace (k+1 x Nfaces) and pressure mean (1 x Nelts)
end
d1 = (size(local_operators{2},2) - 1)/6;
d2 = (size(local_operators{2},1) - 1)/7;
N_FACES = size(mesh.faces, 2);
N_ELEMENTS = size(mesh.elements, 2);
N_DOF = d1*N_FACES;  % # of DOFs belonging to one scalar variable

u_hat = x_hat(1:(2*N_DOF));
rho = x_hat((1:N_ELEMENTS) + 2*N_DOF).';

ux_hat = reshape(u_hat(1:N_DOF), d1, N_FACES);
uy_hat = reshape(u_hat(N_DOF+1:end), d1, N_FACES);

u_aux = [reshape(ux_hat(:, mesh.boundaries), [3*d1, N_ELEMENTS]);
         reshape(uy_hat(:, mesh.boundaries), [3*d1, N_ELEMENTS])
         rho];
     
[A1, A2, Af] = deal(local_operators{:});
% page-wise over elements
rhs = pagemtimes(A2, reshape(u_aux, [], 1, N_ELEMENTS)) + Af;
sol = reshape(pagemldivide(A1, rhs), [], N_ELEMENTS);

blockd2 = @(x) d2*(x-1) + (1:d2);
L11_h = sol(blockd2(1),:);
L21_h = sol(blockd2(2),:);
L12_h = sol(blockd2(3),:);
L22_h = sol(blockd2(4),:);
u1_h = sol(4*d2 + blockd2(1),:);
u2_h = sol(4*d2 + blockd2(2),:);
p_h = sol(4*d2 + 2*d2 + (1:d2),:);

local_solutions = {{L11_h, L12_h; L21_h, L22_h}, {u1_h; u2_h}, p_h};
trace_solutions = {{ux_hat; uy_hat}, rho};
end