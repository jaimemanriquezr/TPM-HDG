function [lhs, flux, rhs] = BiotMatrices(k, mesh, E, nu, kappa, options)
% HDG quasi-static Biot (total stress, displacement, rotation, Darcy flux,
% pore pressure) as IndexedMatrix blocks, one implicit time step:
%   local  spaces (ids 1-5): sigma = PEERSBasis(k) (dS), total stress;
%       eta = Pk.to_vector() (2 d2), displacement; gamma = Pk (d2),
%       rotation; q = Pk.to_vector() (2 d2), Darcy flux; p = Pk (d2),
%       pore pressure
%   global space  (id 6):    (eta_hat_x, eta_hat_y, p_hat) =
%       LegendreBasis(k).to_vector(3), per face [eta_hat_x; eta_hat_y; p_hat]
%
% With c_b = alpha/(2(lambda + mu)), cdt the effective time step and n the
% outward normal of K, the element equations are
%   (A sigma, tau) + c_b (p, tr tau) + (eta, div tau) + (gamma, tau) - <eta_hat, tau n> = 0
%   -(div sigma, v) + b0 (eta, v) + <tau_s (eta - eta_hat), v>                         = (f, v)
%   (sigma, zeta)                                                                      = 0
%   (kappa^-1 q, w) - (p, div w) + <p_hat, w.n>                                        = 0
%   (1/cdt)[(s0 + 2 alpha c_b)(p, r) + c_b (tr sigma, r)] + (div q, r)
%       + <tau_p (p - p_hat), r>  = (g, r) + (1/cdt)[(s0 + 2 alpha c_b)(p*, r) + c_b (tr sigma*, r)]
% and the trace rows <sigma n - tau_s (eta - eta_hat), mu> (components 1-2)
% and -<q.n + tau_p (p - p_hat), r_hat> (component 3), with tau_s =
% stabilization_solid * mesh.tau and tau_p = stabilization_darcy *
% ||kappa|| * mesh.tau (||.||
% the spectral norm: tau_p carries the units of the permeability, so the
% Darcy block is invariant under kappa -> c kappa up to the factor c).
%
% With density rho_P > 0 the inertia rho_P d_tt eta is the backward
% difference applied twice, v = (eta - eta*)/cdt, a = (v - v*)/cdt: the
% momentum row gains rho_P/cdt^2 (eta, v) (added to b0) and the history
% rho_P (eta*/cdt^2 + v*/cdt, v) on its right-hand side.
%
% lhs holds the blocks tested by ids 1-5, flux those tested by id 6, rhs
% the [-1, 2] and [-1, 5] history blocks (empty when the history is zero). The volume
% sources (f, g) are not added here. Hybrid blocks carry the face index in
% dim 7 (trials) / dim 6 (tests).
%
% The elastic part is LinearElasticityMatrices (ids 1-3, trace components
% 1-2) and the Darcy part DarcyMatrices with its tau scaled (ids 4-5,
% trace component 3); with biot_alpha = 0 and storage = 0 the blocks are
% their direct sum.
%
% options.history carries the time-scheme combinations w*: fields stress
% (dS x n_el), pressure (d2 x n_el), and, read only when density > 0,
% displacement and velocity (2 d2 x n_el); [], the scalar 0 or a missing
% field mean zero.
arguments
    k     (1,1) double
    mesh        SimplicialGrid
    E     (1,1) double
    nu    (1,1) double
    kappa (2,2) double = eye(2)
    options.biot_alpha          (1,1) double = 1.0
    options.storage             (1,1) double = 0.0
    options.effective_time_step (1,1) double = 1.0
    options.history             struct = struct.empty
    options.hooke_modulus       (1,1) double = 0.0
    options.stabilization_darcy (1,1) double = 1.0
    options.stabilization_solid (1,1) double = 1.0
    options.density             (1,1) double {mustBeNonnegative} = 0.0
end
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
c_b = options.biot_alpha/(2*(lambda + mu));
cdt = options.effective_time_step;
s0 = options.storage;

S_k = PEERSBasis(k);
P_k = DubinerBasis(k);
d1 = k + 1;
d2 = P_k.Dimension;

%% elastic part: ids 1-3 kept, the trace (id 4) becomes components 1-2 of id 6
rho = options.density;
mesh_s = mesh;
mesh_s.tau = options.stabilization_solid * mesh.tau;
[e_lhs, e_flux] = LinearElasticityMatrices(k, mesh_s, E, nu, ...
    hooke_modulus = options.hooke_modulus + rho/cdt^2);
%% Darcy part: ids 1-2 -> 4-5, the trace (id 3) becomes component 3 of id 6
mesh_p = mesh;
mesh_p.tau = BiotProblem.darcy_tau_factor(options.stabilization_darcy, kappa) * mesh.tau;
[d_lhs, d_flux] = DarcyMatrices(k, mesh_p, kappa);

elastic_ids = [1 2 3 6];
darcy_ids   = [4 5 6];
lhs  = [remap(e_lhs.Arrays,  elastic_ids, 4, 0, d1); remap(d_lhs.Arrays,  darcy_ids, 3, 2*d1, 0)];
flux = [remap(e_flux.Arrays, elastic_ids, 4, 0, d1); remap(d_flux.Arrays, darcy_ids, 3, 2*d1, 0)];

%% coupling: c_b (p, tr tau), (1/cdt)[(s0 + 2 alpha c_b)(p, r) + c_b (tr sigma, r)]
M_PS = trace_mass(S_k, k, mesh);                                              % (p_j, tr tau_i)
M_SP = pagetranspose(M_PS);                                                   % (tr sigma_j, r_i)
M_PP = pages(mass_matrix(mesh, P_k, P_k));
m_p = (s0 + 2*options.biot_alpha*c_b)/cdt;
if c_b ~= 0
    lhs = [lhs; {IndexedArray(c_b*M_PS, [5 1]); IndexedArray((c_b/cdt)*M_SP, [1 5])}];
end
if m_p ~= 0
    lhs = [lhs; {IndexedArray(m_p*M_PP, [5 5])}];
end
lhs  = sum_blocks(lhs);
flux = sum_blocks(flux);

%% history
n_el = size(mesh.elements, 2);
sigma_star = history_field(options.history, "stress");
p_star = history_field(options.history, "pressure");
rhs = [];
V = zeros(d2, 1, 1, 1, n_el);
if c_b ~= 0 && ~isempty(sigma_star)
    V = V + (c_b/cdt)*pagemtimes(M_SP, reshape(sigma_star, S_k.Dimension, 1, 1, 1, n_el));
end
if m_p ~= 0 && ~isempty(p_star)
    V = V + m_p*pagemtimes(M_PP, reshape(p_star, d2, 1, 1, 1, n_el));
end
if any(V(:))
    rhs = IndexedMatrix(IndexedArray(V, [-1, 5]));
end
if rho > 0
    % inertia history rho (eta*/cdt^2 + v*/cdt, v)
    eta_star = history_field(options.history, "displacement");
    v_star = history_field(options.history, "velocity");
    W = zeros(2*d2, n_el);
    if ~isempty(eta_star), W = W + eta_star/cdt^2; end
    if ~isempty(v_star),   W = W + v_star/cdt;     end
    if any(W(:))
        M_UU = pagekron(eye(2), M_PP);
        H = rho*pagemtimes(M_UU, reshape(W, 2*d2, 1, 1, 1, n_el));
        if isempty(rhs)
            rhs = IndexedMatrix(IndexedArray(H, [-1, 2]));
        else
            rhs = rhs + IndexedMatrix(IndexedArray(H, [-1, 2]));
        end
    end
end
end

function out = remap(arrays, new_ids, trace_id, before, after)
% Relabel the blocks with new_ids(old id) and pad the trace rows/columns
% (old id trace_id) with `before` zero rows/columns ahead and `after`
% behind, so that they sit in the 3-component trace layout.
out = cell(numel(arrays), 1);
for i = 1:numel(arrays)
    A = double(arrays{i});
    id = arrays{i}.id;
    if id(1) == trace_id, A = pad(A, 2, before, after); end
    if id(2) == trace_id, A = pad(A, 1, before, after); end
    new_id = id;
    new_id(id > 0) = new_ids(id(id > 0));
    out{i} = IndexedArray(A, new_id);
end
end

function M = trace_mass(S_k, k, mesh)
% (p_j, tr tau_i) on every element, [dS, d2, 1, 1, n_el]: tau_i the PEERS
% functions (bubbles included), p_j the Dubiner functions of degree k
[x_local, w] = get_2d_reference_quadrature(k + 1);
vol = permute(mesh.volume, [1 3 4 5 2]);
S = S_k.evaluate(x_local, mesh);
tr_S = S(:, :, 1, 1, :) + S(:, :, 2, 2, :);
P = DubinerBasis.evaluate_basis(k, x_local(1, :), x_local(2, :));
M = vol .* pagemtimes(w .* tr_S, 'transpose', P, 'none');
end

function A = pad(A, dim, before, after)
sz = size(A, 1:max(ndims(A), dim));
z_before = sz;  z_before(dim) = before;
z_after  = sz;  z_after(dim)  = after;
A = cat(dim, zeros(z_before), A, zeros(z_after));
end

function M = sum_blocks(blocks)
M = IndexedMatrix(blocks{1});
for i = 2:numel(blocks)
    M = M + IndexedMatrix(blocks{i});
end
end

function w = history_field(history, name)
% a history dof array, or [] when absent or zero
w = [];
if ~isempty(history) && isfield(history, name)
    w = history.(name);
    if isequal(w, 0) || ~any(w(:))
        w = [];
    end
end
end
