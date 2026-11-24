function [local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, polynomial_degree)
%-------------------------------------------------------------------------
%Description: Computes the global and local matrices for the Stokes problem
%             Find (L,U,p) \in V^(2x2)xV^2xV
%                  L - ∇U = 0
%            -∇·(ν*L - p*Id) = f             (PDE)         
%                   ∇·U = 0                  
%                    ∫p = 0                  (PRS)
% whose HDG discretization is
%                    (L, G) + (u, ∇·G) = <û,Gn>
%       -ν*(∇·L,V) + ν*<SU,V> + (∇p,V) = <Sû,V> + (f,V)
%                       (U,∇q) + (ξ,q) = <û·n,q>
%                                (p,η) = (ρ,η)
%        〚(ν*L - p*Id - νS(U - û))·n〛= 0
%                               〚û·n〛= 0
% and, in matrix form,
%              A1*[L;U;p;ξ] = A2*[û;ρ] + Af     (PDE_h)
%              A3*[L;U;p;ξ] + tLL*û = 0         (NUMFLUX NORMAL TRACE) 
%              [nxLP^t nyLP^t]*û = 0            (VELOCITY NORMAL TRACE)      
%-------------------------------------------------------------------------
arguments (Input)
    pde               (1,1)   FEProblem      % HDG problem object
    mesh              SimplicialGrid          % triangular mesh
    polynomial_degree (1,1)   double          = pde.Order  % polynomial degree
end
arguments (Output)
    local_operators (3,1)   cell    % {A1, A2, Af}: local-to-global operators per element
    lhs_entries     (:,:,:) double  % LHS Schur complement (2*3*(k+1)+1 x 2*3*(k+1)+1 x Nelts)
    rhs_entries     (:,:)   double  % RHS Schur complement (2*3*(k+1)+1 x Nelts)
end
k = polynomial_degree;
data = pde.Data;

f = data.source_term;
nu = data.viscosity;
% Oseen (Solano & Vargas 2022, eq. 3.4): a convective velocity beta adds
% -(u x beta, grad v) + <(uhat x beta) n, v> to the momentum row, -(uhat x
% beta) n to the numerical flux, and scales the stabilization to
% S = nu*tau with tau = max|beta.n|/(2 nu) + 1 (their (3.5)). Absent, the
% Stokes path below is untouched (S = nu).
is_oseen = isfield(data, "convective_velocity") && ~isempty(data.convective_velocity);
S = nu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRELIMINARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USEFUL QUANTITIES
N_ELEMENTS = size(mesh.elements, 2);
d1 = k+1;             % dim P_k(e)
d2 = nchoosek(k+2,2); % dim P_k(K)

% VOLUME INTEGRALS
M = mass_matrix(mesh, k);              % M_ij = \int_K P_i * P_j
[Cx, Cy] = convection_matrix(mesh, k); % C_ij = \int_K P_i * dP_j
CxT = permute(Cx,[2 1 3]); CyT = permute(Cy,[2 1 3]);
[tPP, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);

% BOUNDARY INTEGRALS
nxPL = permute(nxLP,[2 1 3]);
nyPL = permute(nyLP,[2 1 3]);
tPL = permute(tLP,[2 1 3]);

% CONVECTIVE INTEGRALS (Oseen only)
if is_oseen
    [Cb, BbPL, BbLL, max_bn] = convective_matrix(mesh, k, data.convective_velocity, skip_faces = mesh.neufaces);
    S = nu*(max_bn/(2*nu) + 1);
    BbLP = permute(BbPL, [2 1 3]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (L, G)_K + (u, div(G))_K = <\uhat,Gn>_dK
% -nu*(div(L),v)_K + nu*<tu,v>_dK + (grad(p),v)_K = <t\uhat,v>_dK + (f,v)_K
%  (u,grad(q))_K + (xi,q)_K = <\uhat*n,q>_dK
%                  (p,eta)_K = (\rho,eta)_K

O = zeros(d2,d2,N_ELEMENTS);
o = zeros(d2,3*d1,N_ELEMENTS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% EQ 1: (L, G) + (u, div(G) = <\uhat,Gn> %%%%%%%%%%%%%%%%%%%%%
% (L, G)_K
G_L = [M, O, O, O;
       O, M, O, O;
       O, O, M, O;
       O, O, O, M];
% (u, div(G))_K
G_U = [CxT, O ;
        O ,CxT;
       CyT, O ;
        O, CyT];
% <\uhat,Gn>_dK
G_u = [nxPL, o;
      o, nxPL;
      nyPL, o;
      o, nyPL];

% ----------------------------------------------------------------------- %
G_p = [O; 
       O; 
       O; 
       O];
G_xi = [O(:,1,:); 
       O(:,1,:); 
       O(:,1,:); 
       O(:,1,:)];
G_rho = [O(:,1,:);
         O(:,1,:);
         O(:,1,:);
         O(:,1,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% EQ 2: -nu*(div(L),V) + <nu*tU,V> + (grad(p),V) = <nu*t\uhat,V> + (f,V) %%
% -nu*(div(L), V)_K
V_L = -nu*[Cx, O, Cy, O;
           O,  Cx, O, Cy];
% S*<tU,V>_dK
V_U = [S*tPP, O;
          O, S*tPP];
% (grad(p), V)_K
V_p = [Cx;
       Cy];
% <St\uhat,V>_dK
V_u = [S*tPL, o;
          o, S*tPL];
% -(u x beta, grad(V))_K + <(\uhat x beta) n, V>_dK on the left, so the
% face term moves to the A2 (right-hand) side with its sign flipped
if is_oseen
    V_U = V_U - [Cb, O; O, Cb];
    V_u = V_u - [BbLP, o; o, BbLP];
end

% (f,v)_K
F_V = [test_elements(mesh, f{1}, k); 
       test_elements(mesh, f{2}, k)];

% ----------------------------------------------------------------------- %
V_xi = [O(:,1,:); 
        O(:,1,:)];
V_rho = [O(:,1,:); 
        O(:,1,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%% EQ 3: (u,grad(q)) + (xi,q) = <\uhat*n,q> %%%%%%%%%%%%%%%%%%%
% (u, grad(q))_K
q_U = [CxT, CyT];
% (xi, q)_K
q_xi = M(:,1,:);
% <\uhat*n,q>_dK      
q_u = [nxPL, nyPL];

% ----------------------------------------------------------------------- %
q_L = [O, O, O, O];
q_p = O;
q_rho = O(:,1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% EQ 4: (p,eta) = (rho,eta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (p,eta)_K
e_P = M(1,:,:);
% (\rho,eta)_K
e_rho = M(1,1,:);

% ----------------------------------------------------------------------- %
e_L = [O(1,:,:), O(1,:,:), O(1,:,:), O(1,:,:)];
e_U = [O(1,:,:), O(1,:,:)];
e_xi = O(1,1,:);
e_u = [o(1,:,:), o(1,:,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = [G_L, G_U, G_p, G_xi;
      V_L, V_U, V_p, V_xi;
      q_L, q_U, q_p, q_xi;
      e_L, e_U, e_P, e_xi];

A2 = [G_u, G_rho;
      V_u, V_rho;
      q_u, q_rho;
      e_u, e_rho];

Af = [zeros(4*d2,1,N_ELEMENTS);  % G
     F_V;                   % V
     zeros(d2,1,N_ELEMENTS);     % q
     zeros(1,1,N_ELEMENTS)];     % eta

local_operators = {A1, A2, Af};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O = zeros(3*d1,d2,N_ELEMENTS);
% \sum_K <nu*\Lhat*n - \phat*n, mu>_dK = 0
% nu*\Lhat*n - \phat*n = nu*L*n - p*n - S(u - û) [- (û x beta) n]
A3 = [nu*nxLP, O, nu*nyLP, O, -S*tLP, O, -nxLP, O(:,1,:);
      O, nu*nxLP, O, nu*nyLP, O, -S*tLP, -nyLP, O(:,1,:)];

% lhs_entries*[\uhat;rho] = rhs_entries
% First 3*2*d1 equations: [[L*n - p*n, mu]] = 0
% 3*2*d1+1-th equation: \int_dK \uhat*n  = 0

O = zeros(3*d1,3*d1,N_ELEMENTS);
vtLL = [S*tLL, O, O(:,1,:);
         O, S*tLL, O(:,1,:)];
if is_oseen
    vtLL = vtLL - [BbLL, O, O(:,1,:);
                   O, BbLL, O(:,1,:)];
end
nPL = [nxPL(1,:,:), nyPL(1,:,:), O(1,1,:)];
% page-wise over elements: one LU of A1 per element serves both right-hand sides
X = pagemldivide(A1, cat(2, A2, Af));
lhs_entries = cat(1, pagemtimes(A3, X(:, 1:end-1, :)) + vtLL, nPL);
rhs_entries = cat(1, -pagemtimes(A3, X(:, end, :)), zeros(1, 1, N_ELEMENTS));
end