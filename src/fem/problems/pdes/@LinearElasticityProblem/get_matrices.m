function [local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, polynomial_degree)
%-------------------------------------------------------------------------
%Description: Computes the global and local matrices 
% for the linear elasticity problem:
%       Find (σ,u) ∈ S^(2x2) x V^2
%                  ζ⁻¹σ = ε(u)          (CAUCHY)
%                  -∇·σ = f             (PDE)         
% which is equivalent to:
%       Find (σ,u,γ) \in S^(2x2) x V^2 x A^(2x2)
%  (1/2μ)(σ - λ/(2μ + 2λ)tr(σ)Id) - ∇u + γ = 0          (CAUCHY)
%                                     -∇·σ = f             (PDE)   
% whose HDG discretization is
%       Find (σ,u,γ) ∈ S^(2x2) x V^2 x A^(2x2)
%  c1*(σ,τ) - c2*(tr(σ)Id,τ) + (u,∇·τ) + (γ,τ) = <û,τn>
%                            -(∇·σ,v) + <Su,v> = <Sû,v> + (f,v)
%                                        (σ,η) = 0
% and, in matrix form,
%              A1*[σ;u;γ] = A2*û + Af     (PDE_h)
%              A3*[σ;u;γ] + SLL*û = 0     (NUMFLUX NORMAL TRACE)           
%-------------------------------------------------------------------------
arguments (Input)
    pde               (1,1)   FEProblem      % HDG problem object
    mesh              SimplicialGrid          % triangular mesh
    polynomial_degree (1,1)   double          = pde.Order  % polynomial degree
end
arguments (Output)
    local_operators (3,1)   cell    % {A1, A2, Af}: local-to-global operators per element
    lhs_entries     (:,:,:) double  % LHS Schur complement (2*3*(k+1) x 2*3*(k+1) x Nelts)
    rhs_entries     (:,:)   double  % RHS Schur complement (2*3*(k+1) x Nelts)
end
k = polynomial_degree;
data = pde.Data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRELIMINARIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USEFUL QUANTITIES
N_ELEMENTS = size(mesh.elements, 2);
d1 = k+1;   % dim P_k(e)
d2 = nchoosek(k+2,2); % dim P_k(K)
if k == 0
    d2star = 0;
else
    d2star = nchoosek(k+1,1);
end
d2peers = 4*d2 + d2star;


% VOLUME INTEGRALS
M = mass_matrix(mesh, k);
[M_PEERS, M_tr_PEERS, M_A_PEERS] = PEERS_matrix(mesh, k);
[Cx, Cy] = convection_matrix(mesh, k); % C_ij = \int_K P_i * dP_j
CxT = permute(Cx, [2 1 3]); 
CyT = permute(Cy, [2 1 3]);

% BOUNDARY INTEGRALS
[tPP, nxLP, nyLP, tLP, tLL] = boundary_matrix(mesh, k);
nxPL = permute(nxLP, [2 1 3]);
nyPL = permute(nyLP, [2 1 3]);
tPL = permute(tLP, [2 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_time_dependent = (data.density > 0);

beta0 = data.hooke_modulus; % elastic coefficient
dt = data.time_step;

f = data.source_term;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c1*(σ,t)_K - c2*(tr(σ)Id,t)_K + (U,div(t))_K + (g,t)_K = <u,tn>_dK
%                              -(div(σ),v)_K + <SU,v>_dK = <Su,v>_dK + (f,v)_K
%                                             (σ,e)_K = 0
% where c1 = (1/2*mu), c2 = (1/2*mu)*(lambda/(4*lambda + 2*mu))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% EQ 1: c1*(σ,t) - c2*(tr(σ)Id,t) + (u,∇·t) + (g,t) = <û,tn> %%%%%%%%
% c1*(σ,τ)_K - c2*(tr(σ)Id,τ)_K
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(data.young_modulus, data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(4*mu*(lambda + mu));
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

% <\uhat,tn>_dK
O = zeros(d2, 3*d1,N_ELEMENTS);
o = zeros(d2star, 3*d1,N_ELEMENTS);
t_u = [nxPL, O;
      O, nxPL;
      nyPL, O;
      O, nyPL;
      o, o];
% ----------------------------------------------------------------------- %
t_f = zeros(d2peers,1,N_ELEMENTS);
% optional right-hand side of the constitutive row: data.constitutive_source
% (d2peers x N_ELEMENTS) holds, per element, the moments (r, t)_K against
% the PEERS test functions; absent or empty, the row is homogeneous
if isfield(data, "constitutive_source") && ~isempty(data.constitutive_source)
    t_f = t_f + reshape(data.constitutive_source, d2peers, 1, N_ELEMENTS);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% EQ 2: -(div(σ),v) + <SU,v> = <Su,v> + (f,v) %%%%%%%%%%%%%
O = zeros(d2, d2,N_ELEMENTS);
o = zeros(d2, 3*d1,N_ELEMENTS);
% -(div(σ), V)_K
V_s = -permute(t_U, [2 1 3]);
% b0*(U, V)_K + <tU,V>_dK
V_U = beta0*[M, O; O, M] + [tPP, O; O, tPP];
% <t\uhat,V>_dK
V_u = [tPL, o;
       o, tPL];

% (f,v)_K
V_f = [test_elements(mesh, f{1}, k); 
       test_elements(mesh, f{2}, k)];
% ----------------------------------------------------------------------- %
V_g = [O;
       O]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% EQ 3: (σ,e) = 0 %%%%%%%%%%%%%
O = zeros(d2, d2,N_ELEMENTS);
o = zeros(d2, 3*d1,N_ELEMENTS);
% (σ, e)_K
e_s = permute(t_g, [2 1 3]); 

% ----------------------------------------------------------------------- %
e_U = [O, O];
e_g = O;
e_u = [o, o];
e_f = zeros(d2,1,N_ELEMENTS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if is_time_dependent
    mass = data.density*[M, O; O, M];
    % each history is a dof array or zero ([] or the scalar 0), read
    % independently from data.state.solid
    history = zeros(2*d2, N_ELEMENTS);
    h = solid_history(data, history);

    if isfield(data, "time_scheme") && data.time_scheme == "newmark"
        % Newmark-β in displacement form: the acceleration
        %   a^{n+1} = (u^{n+1} - u^n - Δt v^n - Δt²(1/2-β) a^n) / (β Δt²)
        % is eliminated from M a^{n+1} + K u^{n+1} = f^{n+1}.
        beta = data.newmark_beta;
        V_U = V_U + mass/(beta*dt^2);
        u_0  = reshape(h.displacement, 2*d2, 1, []);
        du_0 = reshape(h.velocity,     2*d2, 1, []);
        V_f = V_f + pagemtimes(mass, u_0 + dt*du_0)/(beta*dt^2);
        a_0 = reshape(h.acceleration, 2*d2, 1, []);
        V_f = V_f + (1/2 - beta)/beta * pagemtimes(mass, a_0);
    elseif onestep_start(data)
        V_U = V_U + 2*mass/dt^2;
        du_0 = reshape(h.velocity, 2*d2, 1, []);
        u_0 = reshape(h.displacement, 2*d2, 1, []); 
        V_f = V_f + 2*pagemtimes(mass, u_0)/dt^2 + 2*pagemtimes(mass, du_0)/dt;
    else
        
        V_U = V_U + mass/dt^2;
        u_00 = reshape(h.displacement_previous, 2*d2, 1, []); 
        u_0 = reshape(h.displacement, 2*d2, 1, []); 
        V_f = V_f + pagemtimes(mass, 2*u_0 - u_00)/dt^2;
    end
end

A1 = [t_s, t_U, t_g;
      V_s, V_U, V_g;
      e_s, e_U, e_g];

A2 = [t_u;
      V_u;
      e_u];      
  
Af = [t_f;  
      V_f;                   
      e_f];

local_operators = {A1, A2, Af};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O = zeros(3*d1, d2, N_ELEMENTS);
o = zeros(3*d1, d2star, N_ELEMENTS);
% \sum_K <\sigmahat*n, mu>_dK = 0
% \sigmahat*n = \sigma*n - tau*(u - \uhat)
A3 = [nxLP,    O, nyLP,    O, o,  -tLP,    O,  O;
         O, nxLP,    O, nyLP, o,     O, -tLP,  O];

% lhs_entries*[\uhat] = rhs_entries
% First 3*2*d1 equations: [[sigma*n,mu]] = 0

O = zeros(3*d1, 3*d1, N_ELEMENTS);
bigtLL = [tLL, O;
         O, tLL];
% page-wise over elements: one LU of A1 per element serves both right-hand sides
X = pagemldivide(A1, cat(2, A2, Af));
lhs_entries = pagemtimes(A3, X(:, 1:end-1, :)) + bigtLL;
rhs_entries = -pagemtimes(A3, X(:, end, :));
end

function h = solid_history(data, template)
    % The volume histories u⁰, v⁰, a⁰, u⁻¹ as arrays shaped like the
    % template, from data.state.solid; without a state the history is zero.
    names = ["displacement", "velocity", "acceleration", "displacement_previous"];
    source = [];
    if isfield(data, "state") && ~isempty(data.state)
        source = data.state.solid;
    end
    for i = 1:numel(names)
        value = [];
        if ~isempty(source) && isfield(source, names(i)) && ~iscell(source.(names(i)))
            value = source.(names(i));
        end
        h.(names(i)) = default_zeros(value, template);
    end
end

function tf = onestep_start(data)
    % The first one-step step (Taylor start): a state whose
    % displacement_previous is unset.
    tf = isfield(data, "state") && ~isempty(data.state) && isfield(data.state, "solid") ...
        && ~isempty(data.state.solid) && (~isfield(data.state.solid, "displacement_previous") ...
        || isempty(data.state.solid.displacement_previous));
end

function A = default_zeros(A, template)
    % Empty or scalar-0 histories expand to a zero array shaped like the template.
    if isempty(A) || isequal(A, 0.0)
        A = zeros(size(template));
    end
end
