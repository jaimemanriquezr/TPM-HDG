function [local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, polynomial_degree)
%-------------------------------------------------------------------------
%Description: Computes the local operators and the global (Schur complement)
% matrices for the incompressible elasticity (Stokes-type) problem in
% stress-velocity form, with the PEERS element for the stress:
%       Find (σ,u,γ) ∈ S^(2x2) x V^2 x A
%  (1/2μ)(σ - 1/2 tr(σ)Id) - ∇u + γ = 0          (CAUCHY)
%                              -∇·σ = f          (PDE)
%                              (σ,η) = 0         (SYMMETRY)
% whose HDG discretization, with the hybrid unknowns (û,ρ) (velocity trace
% and per-element constant ρ) and a per-element constant multiplier ξ, is
%  c1*(σ,τ) - c2*(tr(σ)Id,τ) + (u,∇·τ) + (γ,τ) + (ξ,tr(τ)) = <û,τn>
%                        -(∇·σ,v) + <Su,v> = <Sû,v> + (f,v)
%                                      (σ,η) = 0
%                                (tr(σ),z) = (ρ,z)        (z constant)
% with c1 = 1/(2μ), c2 = 1/(4μ) and S = μ. For density > 0 the momentum row
% gains the θ-scheme mass and history terms; with data.convective_velocity
% set it gains the convective terms and S changes accordingly (see the
% body). In matrix form,
%              A1*[σ;u;γ;ξ] = A2*[û;ρ] + Af     (PDE_h)
%              A3*[σ;u;γ;ξ] + S*tLL*û = 0      (NUMFLUX NORMAL TRACE)
% where the last row of the global system is ∫_∂K û·n = 0 on each element.
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
is_time_dependent = (data.density > 0);
dt = data.time_step;

f = data.source_term;
mu = data.viscosity;
% Navier-Stokes/Oseen (Solano & Vargas 2022, eq. 3.4, kept out of the
% stress): a convective velocity beta adds rho_c[-(u x beta, grad v) +
% <(uhat x beta) n, v>] to the momentum row and -rho_c (uhat x beta) n to
% the numerical flux, with rho_c = density when transient and 1 when
% steady (kinematic form), and scales the stabilization to S = mu*tau,
% tau = rho_c max|beta.n|/(2 mu) + 1. Without it the Stokes operator is unchanged.
is_oseen = isfield(data, "convective_velocity") && ~isempty(data.convective_velocity);
S = mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRELIMINARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USEFUL QUANTITIES
N_ELEMENTS = size(mesh.elements, 2);
d1 = k+1;   % dim P_k(e)
d2 = nchoosek(k+2,2); % dim P_k(K)

if k == 0
    d2star = 0;
else
    d2star = nchoosek(k + 1, 1);
end
d2peers = 4*d2 + d2star;

c1 = 1/(2*mu);
c2 = 1/(4*mu);

% VOLUME INTEGRALS
M = mass_matrix(mesh, k);              % M_ij = \int_K P_i * P_j
[M_PEERS, M_tr_PEERS, M_A_PEERS] = PEERS_matrix(mesh, k);
[Cx, Cy] = convection_matrix(mesh, k); % C_ij = \int_K P_i * dP_j
CxT = permute(Cx,[2 1 3]); CyT = permute(Cy,[2 1 3]);
[tPP, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);


% BOUNDARY INTEGRALS
nxPL = permute(nxLP,[2 1 3]);
nyPL = permute(nyLP,[2 1 3]);
tPL = permute(tLP,[2 1 3]);

% CONVECTIVE INTEGRALS (Oseen only)
if is_oseen
    rho_c = 1;
    if is_time_dependent; rho_c = data.density; end
    % do-nothing on every face with a natural (traction) component: the
    % "Neumann" tag, or the componentwise specs of data.boundary_conditions
    natural_faces = mesh.neufaces(:);
    if isfield(data, "boundary_conditions")
        for s = data.boundary_conditions(:).'
            if any(s.type == "neumann") && isKey(mesh.boundary_tags, s.tag)
                natural_faces = [natural_faces; reshape(get_face_by_tag(mesh, s.tag), [], 1)]; %#ok<AGROW>
            end
        end
    end
    [Cb, BbPL, BbLL, max_bn] = convective_matrix(mesh, k, data.convective_velocity, skip_faces = natural_faces);
    Cb = rho_c*Cb;  BbPL = rho_c*BbPL;  BbLL = rho_c*BbLL;
    % convective_stabilization = false leaves S = mu*(mesh tau): the caller
    % has put tau(beta) on the mesh so that transmission rows built from
    % mesh.tau (InterfaceProblem) match the local problem
    if ~isfield(data, "convective_stabilization") || data.convective_stabilization
        S = mu + rho_c*max_bn/2;
    end
    BbLP = permute(BbPL, [2 1 3]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c1*(σ,t)_K - c2*(tr(σ)Id,t)_K + (U,div(t))_K + (g,t)_K = <u,tn>_dK
%                              -(div(σ),v)_K + <SU,v>_dK = <Su,v>_dK + (f,v)_K
%                                             (σ,e)_K = 0
% where c1 = 1/(2*mu), c2 = 1/(4*mu), so c1*σ - c2*tr(σ)Id = (1/2*mu)*dev(σ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% EQ 1: k1*(σ,t) - k2*(tr(σ)Id,t) + (u,∇·t) + (g,t) + (xi, tr(t))= <û,tn> %%%%%%%%
% k1*(σ,τ)_K + k2*(tr(σ)Id,τ)_K
t_s = c1*M_PEERS - c2*M_tr_PEERS;

% (u, div(t))_K
O = zeros(d2, d2,N_ELEMENTS);
o = zeros(d2star, d2,N_ELEMENTS);
t_U = [CxT, O ;
        O, CxT;
       CyT, O ;
        O, CyT;
        o, o];

% (g, t)_K
t_g = M_A_PEERS;

% (xi, t)_K
t_xi = M_tr_PEERS(:, 1, :);

% <\uhat,tn>_dK
O = zeros(d2, 3*d1,N_ELEMENTS);
o = zeros(d2star, 3*d1,N_ELEMENTS);
t_u = [nxPL, O;
      O, nxPL;
      nyPL, O;
      O, nyPL;
      o, o];
% ----------------------------------------------------------------------- %
t_rho   = zeros(d2peers, 1, N_ELEMENTS);
t_f     = zeros(d2peers,1,N_ELEMENTS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% EQ 2: -(div(σ),v) + <SU,v> = <Su,v> + (f,v) %%%%%%%%%%%%%
O = zeros(d2, d2,N_ELEMENTS);
o = zeros(d2, 3*d1,N_ELEMENTS);
% -(div(σ), V)_K
V_s = -permute(t_U, [2 1 3]);
% <StU,V>_dK
V_U = [S*tPP, O; O, S*tPP];
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
V_f = [test_elements(mesh, f{1}, k); 
       test_elements(mesh, f{2}, k)];
% ----------------------------------------------------------------------- %
V_xi    = zeros(2*d2, 1, N_ELEMENTS);
V_g     = zeros(2*d2, d2, N_ELEMENTS); 
V_rho   = zeros(2*d2, 1, N_ELEMENTS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% EQ 3: (σ,e) = 0 %%%%%%%%%%%%%
% (σ, e)_K
e_s = permute(t_g, [2 1 3]); 

% ----------------------------------------------------------------------- %
e_U     = zeros(d2, 2*d2, N_ELEMENTS);
e_g     = zeros(d2, d2, N_ELEMENTS);
e_u     = zeros(d2, 2*3*d1, N_ELEMENTS);
e_xi    = zeros(d2, 1, N_ELEMENTS);
e_rho   = zeros(d2, 1, N_ELEMENTS);
e_f     = zeros(d2,1,N_ELEMENTS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% EQ 4: (tr(σ),z) = (ρ,z) %%%%%%%%%%%%%
% (tr(σ),z)_K
z_s     = permute(t_xi, [2 1 3]); 
z_rho   = M_PEERS(1, 1, :);
% ----------------------------------------------------------------------- %
z_U     = zeros(1, 2*d2, N_ELEMENTS);
z_g     = zeros(1, d2, N_ELEMENTS);
z_u     = zeros(1, 2*3*d1, N_ELEMENTS);
z_xi    = zeros(1, 1, N_ELEMENTS);
z_f     = zeros(1, 1, N_ELEMENTS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_time_dependent
    mass = data.density*[M, O; O, M];

    % theta-scheme on the momentum equation: theta = 1 is Backward Euler;
    % theta = 1/2 is Crank-Nicolson, written as BE with dt -> theta*dt plus
    % the history term ((1-theta)/theta) M a0, where a0 = u_t(t^n) absorbs
    % the old spatial residual (M a^n = -K u^n + f^n) so the source stays
    % f^{n+1}. theta comes from data.theta / data.time_scheme (get_theta).
    theta = pde.get_theta();
    % u0 and a0 are dof arrays or zero ([] or the scalar 0), read
    % independently from data.state.fluid
    h = fluid_history(data);

    V_U = V_U + mass/(theta*dt);
    if ~is_zero(h.velocity)
        u_0 = reshape(h.velocity, 2*d2, 1, []);
        V_f = V_f + pagemtimes(mass, u_0)/(theta*dt);
    end
    if theta < 1 && ~is_zero(h.acceleration)
        a_0 = reshape(h.acceleration, 2*d2, 1, []);
        V_f = V_f + (1 - theta)/theta * pagemtimes(mass, a_0);
    end
end

A1 = [t_s, t_U, t_g, t_xi;
      V_s, V_U, V_g, V_xi;
      e_s, e_U, e_g, e_xi;
      z_s, z_U, z_g, z_xi];

A2 = [t_u, t_rho;
      V_u, V_rho;
      e_u, e_rho;
      z_u, z_rho];      

Af = [t_f;  
      V_f;                   
      e_f;
      z_f];

local_operators = {A1, A2, Af};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O = zeros(3*d1,d2,N_ELEMENTS);
o = zeros(3*d1,d2star,N_ELEMENTS);
% \sum_K <\sigmahat*n, mu>_dK = 0
% \sigmahat*n = \sigma*n - S*tau*(u - \uhat) [- (\uhat x beta) n]
A3 = [nxLP,    O, nyLP,    O, o,  -S*tLP,    O,  O, O(:,1,:);
         O, nxLP,    O, nyLP, o,     O, -S*tLP,  O, O(:,1,:)];

% lhs_entries*[\uhat;rho] = rhs_entries
% First 3*2*d1 equations: [[sigma*n, mu]] = 0
% 3*2*d1+1-th equation: \int_dK \uhat*n  = 0

O = zeros(3*d1,3*d1,N_ELEMENTS);
bigtLL = [S*tLL, O, O(:,1,:);
          O, S*tLL, O(:,1,:)];
if is_oseen
    bigtLL = bigtLL - [BbLL, O, O(:,1,:);
                       O, BbLL, O(:,1,:)];
end
nPL = [nxPL(1,:,:), nyPL(1,:,:), O(1,1,:)];
% page-wise over elements: one LU of A1 per element serves both right-hand sides
X = pagemldivide(A1, cat(2, A2, Af));
lhs_entries = cat(1, pagemtimes(A3, X(:, 1:end-1, :)) + bigtLL, nPL);
rhs_entries = cat(1, -pagemtimes(A3, X(:, end, :)), zeros(1, 1, N_ELEMENTS));
end

function h = fluid_history(data)
    % The volume histories u0, a0 from data.state.fluid ([] when absent;
    % without a state the history is zero).
    names = ["velocity", "acceleration"];
    source = [];
    if isfield(data, "state") && ~isempty(data.state)
        source = data.state.fluid;
    end
    for i = 1:numel(names)
        h.(names(i)) = [];
        if ~isempty(source) && isfield(source, names(i)) && isnumeric(source.(names(i)))
            h.(names(i)) = source.(names(i));
        end
    end
end

function z = is_zero(A)
    % [] and the scalar 0 stand for a zero history
    z = isempty(A) || isequal(A, 0.0);
end