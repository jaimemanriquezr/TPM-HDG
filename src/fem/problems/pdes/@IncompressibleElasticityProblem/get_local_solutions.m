function [local_solutions, trace_solutions] = get_local_solutions(obj, x_hat, local_operators, mesh)
%-------------------------------------------------------------------------
%Description: Recovers the local volume solutions (σ,u,γ) per element by
%             solving the element-local problem
%               A1*[σ;u;γ] = A2*[û_K; ρ_K] + Af
%             where û_K and ρ_K are the hybrid DOFs on element K.
%-------------------------------------------------------------------------
arguments (Input)
    obj             (1,1)   FEProblem      % incompressible elasticity problem object
    x_hat           (:,1)   double          % global hybrid DOF vector (n_dof x 1)
    local_operators (3,1)   cell            % {A1, A2, Af} from get_matrices
    mesh            SimplicialGrid          % triangular mesh
end
arguments (Output)
    local_solutions (:,1)   cell    % {S_h (PEERS cell), {u1_h;u2_h}, G_h}: stress, velocity, rotation per element
    trace_solutions (:,1)   cell    % {{ux_hat;uy_hat}, rho}: velocity trace (k+1 x Nfaces) and pressure mean (1 x Nelts)
end
d1 = (size(local_operators{2}, 2) - 1) / 6;  % A2 has 6*d1+1 = 2*3*(k+1)+1 columns (+1 for ρ)
d2 = d1*(d1+1)/2;                             % C(k+2,2) = (k+1)(k+2)/2 = d1*(d1+1)/2
N_FACES = size(mesh.faces, 2);
N_ELEMENTS = size(mesh.elements, 2);
N_DOF = d1*N_FACES;

u_hat = x_hat(1:2*N_DOF);
ux_hat = reshape(u_hat(1:N_DOF), d1, N_FACES);
uy_hat = reshape(u_hat(N_DOF+1:end), d1, N_FACES);
rho = reshape(x_hat((1:N_ELEMENTS) + 2*N_DOF), 1, N_ELEMENTS);

u_aux = [reshape(ux_hat(:, mesh.boundaries), 3*d1, N_ELEMENTS);
         reshape(uy_hat(:, mesh.boundaries), 3*d1, N_ELEMENTS);
         rho];
     
[A1, A2, Af] = deal(local_operators{:});
% page-wise over elements
rhs = pagemtimes(A2, reshape(u_aux, [], 1, N_ELEMENTS)) + Af;
sol = reshape(pagemldivide(A1, rhs), [], N_ELEMENTS);

S_h = cell(4, 1);
for i = 1:4
    S_h{i} = sol(d2*(i-1) + (1:d2), :);
end

if d1 > 1
    S_h{5} = sol(4*d2 + (1:d1), :);
    sigma_dim = 4*d2 + d1;
else
    sigma_dim = 4*d2;
end

u_h = {sol(sigma_dim + (1:d2),:);
       sol(sigma_dim + d2 + (1:d2),:)};
G_h = sol(sigma_dim + 2*d2 + (1:d2), :);

local_solutions = {S_h, u_h, G_h};
trace_solutions = {{ux_hat; uy_hat}, rho};
end