function [local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, polynomial_degree)
%-------------------------------------------------------------------------
%Description: Computes the global and local matrices for the Darcy problem
%             Find (u,p) \in V^2 x V s.t.:
%               u + κ*∇p = 0
%                    ∇·u = f
% whose HDG discretization is
%               (κ⁻¹u, v) - (p, ∇·v) = -<p̂,Gn>
%                   (∇·u,q) + <Su,q> = <Sp̂,q> + (f,q)
%                〚(u + S(p - p̂))·n〛= 0
% and, in matrix form,
%                      A1*[u;p] = A2*p̂ + Af     (PDE_h)
%              A3*[u;p] - tLL*p̂ = 0             (NUMFLUX NORMAL TRACE)
%-------------------------------------------------------------------------
arguments (Input)
    pde               (1,1)   FEProblem      % HDG problem object
    mesh              SimplicialGrid          % triangular mesh
    polynomial_degree (1,1)   double          = pde.Order  % polynomial degree
end
arguments (Output)
    local_operators (3,1)   cell    % {A1, A2, Af}: local-to-global operators per element
    lhs_entries     (:,:,:) double  % LHS Schur complement (3*(k+1) x 3*(k+1) x Nelts)
    rhs_entries     (:,:)   double  % RHS Schur complement (3*(k+1) x Nelts)
end
k = polynomial_degree;
data = pde.Data;

if isa(data.diffusion, 'double')
    kappa = data.diffusion;
    det_kappa = kappa(1)*kappa(4) - kappa(2)*kappa(3);
    kappa_inverse = {@(x,y) kappa(4)/det_kappa*x.^0, @(x,y) -kappa(3)/det_kappa*x.^0;
        @(x,y) -kappa(2)/det_kappa*x.^0, @(x,y) kappa(1)/det_kappa*x.^0};
elseif isa(data.diffusion, 'cell')
    kappa = data.diffusion;
    det_kappa = @(x,y) kappa{1}(x,y).*kappa{4}(x,y) - kappa{2}(x,y).*kappa{3}(x,y);
    kappa_inverse = {@(x,y) kappa{4}(x,y)./det_kappa(x,y).*x.^0, @(x,y) -kappa{3}(x,y)./det_kappa(x,y).*x.^0;
        @(x,y) -kappa{2}(x,y)./det_kappa(x,y).*x.^0, @(x,y) kappa{1}(x,y)./det_kappa(x,y).*x.^0};
end
f = data.source;
%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1 = [M, -C;
%       C^t, tPP]
% A2 = [nLP^t;
%       -tLP^t]
% Af = [0;
%       Fb];
% A3 = [nLP, tLP]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1*[q;u] - A2*û = Af - LOCAL PROBLEM
% A3*[q;u] + tLL*û = 0 - FLUX BALANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MASS MATRIX
M = variable_mass_matrix(mesh, k, kappa_inverse);
M = cell2mat(M);

% CONVECTION MATRICES
[Cx, Cy] = convection_matrix(mesh,k);
C = cat(2, Cx, Cy);
C_T = permute(C, [2 1 3]);

% EDGE MATRICES (edge/volume + edge/edge)
[tauPP, nxLP, nyLP, tauLP, tauLL] = boundary_matrix(mesh, k);
nLP = cat(2, nxLP, nyLP);
nLP_T = permute(nLP, [2 1 3]);
tauLP_T = permute(tauLP, [2 1 3]);

% LOAD VECTOR
Fb = test_elements(mesh, f, k);
O = zeros(size(Fb));

% LOCAL MATRICES
A1 = [M -C_T;
    C tauPP];
A2 = [-nLP_T;
    tauLP_T];
Af = [O;
    O;
    Fb];
A3 = [nLP, tauLP];

% COMPUTING GLOBAL SYSTEM AND STORING LOCAL SOLVERS
local_operators = {A1, A2, Af};
% page-wise over elements: one LU of A1 per element serves both right-hand sides
X = pagemldivide(A1, cat(2, A2, Af));
lhs_entries = pagemtimes(A3, X(:, 1:end-1, :)) - tauLL;
rhs_entries = -pagemtimes(A3, X(:, end, :));
end
