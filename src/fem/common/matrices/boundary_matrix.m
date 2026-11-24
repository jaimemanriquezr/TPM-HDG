function varargout = boundary_matrix(mesh, varargin)
% BOUNDARY_MATRIX Element boundary (face) matrices, two call forms.
%   boundary = boundary_matrix(mesh, trial, test, boundary_weights)
%       AbstractFEBasis form (DSL path).
%   [tauPP, nxLP, nyLP, tauLP, tauLL] = boundary_matrix(mesh, k, ...
%       basis_2d=..., basis_1d=...)
%       function-handle form (low-level path, k = degree).
if nargin >= 2 && isa(varargin{1}, 'double')
    [varargout{1:max(nargout, 1)}] = handle_boundary_matrix(mesh, varargin{:});
else
    varargout{1} = basis_boundary_matrix(mesh, varargin{:});
end
end

function boundary = basis_boundary_matrix(mesh, trial, test, boundary_weights)
    arguments
        mesh SimplicialGrid
        trial AbstractFEBasis {mustBeSameElement(trial, mesh)} = DubinerBasis(0);
        test AbstractFEBasis {mustBeSameElement(test, mesh)} = LegendreBasis(0);
        boundary_weights {mustBeBoundaryWeight} = ones(1, 1, 3, size(mesh.elements, 2));
    end
    if isequal(size(boundary_weights, 1:2), [3 size(mesh.elements, 2)])
        boundary_weights = reshape(boundary_weights, [1 1 1 1 size(boundary_weights)]);
    else
        boundary_weights = reshape(boundary_weights, [1 1 size(boundary_weights)]);
    end
    bases = [trial, test];
    [s, w] = get_1d_reference_quadrature(max([bases.Order]) + 1);

    p_trial = evaluate_on_boundary(trial, s, mesh);
    p_test = evaluate_on_boundary(test, s, mesh);

    area_boundary = reshape(mesh.area(mesh.boundaries), 1, 1, 1, 1, 3, []);
    area_weights = area_boundary .* boundary_weights;
    boundary = pagemtimes(w .* p_test, 'transpose', area_weights .* p_trial, 'none');
    boundary = reshape(sum(boundary, 3:4), size(boundary, 1), size(boundary, 2), 3, []);
end

function p = evaluate_on_boundary(basis, s, mesh)
    switch basis.GeometricInputDimension
        case 1
            s_local = {s; 1-s};
            p = basis.evaluate(s_local);
            p = arrayfun(@(mu) p{mu}, mesh.permutation, 'UniformOutput', false);
        case 2
            x_local = {[0*s; 1 - s]; [1 - s; 0*s]; [1 - s; s]};
            p = basis.evaluate(x_local);
    end
    n_points = size(p{1}, 1);
    dof_dim = size(p{1}, 2);
    out_sz_1 = size(p{1}, 3);
    out_sz_2 = size(p{1}, 4);
    p = reshape(cat(5, p{:}), n_points, dof_dim, out_sz_1, out_sz_2, 3, []);
end

function [tauPP, nxLP, nyLP, tauLP, tauLL] = handle_boundary_matrix(mesh, k, options)
arguments
    mesh SimplicialGrid
    k double
    options.basis_2d function_handle = @dubiner_basis
    options.basis_1d function_handle = @legendre_basis
end
%-------------------------------------------------------------------------
%Description: Computes the element-boundary matrices, with P the volume
%             basis, L the face basis, tau = mesh.tau and n the normal:
%             tauPP_ij = <tau P_j, P_i>_{dK}
%             nxLP_ij  = <n_x P_j, L_i>_{dK}
%             nyLP_ij  = <n_y P_j, L_i>_{dK}
%             tauLP_ij = <tau P_j, L_i>_{dK}
%             tauLL_ij = <tau L_j, L_i>_{dK}
%Input:     - mesh: SimplicialGrid (with stabilization, mesh.tau)
%           - k: polynomial degree
%Output:    - tauPP: matrix (d2 x d2 x Nelts)
%           - nxLP: matrix (3*d1 x d2 x Nelts)
%           - nyLP: matrix (3*d1 x d2 x Nelts)
%           - tauLP: matrix (3*d1 x d2 x Nelts)
%           - tauLL: matrix (3*d1 x 3*d1 x Nelts)
%-------------------------------------------------------------------------
[s, w] = get_1d_reference_quadrature(k);
N_ELEMENTS = size(mesh.elements,2);
N_NODES = length(w);

quadrature_2d = [[0*s(:), 1 - s(:)];
                 [1 - s(:), 0*s(:)];
                 [1 - s(:), s(:)]];
p = options.basis_2d(k, quadrature_2d(:,1), quadrature_2d(:,2));
d2 = size(p, 2);

quadrature_1d = [s(:);
                 1 - s(:)];
l = options.basis_1d(k, quadrature_1d(:,1));
d1 = size(l, 2);

[tauLP, nxLP, nyLP] = deal(zeros(3*d1, d2*N_ELEMENTS));
tauLL = zeros(3*d1, 3*d1, N_ELEMENTS);
tauPP = zeros(d2, d2*N_ELEMENTS);

% BOUNDARY MATRIX (FACE/FACE): REFERENCE ELEMENT
ll = (w.*l(1:N_NODES, :))' * l(1:N_NODES,:);

area_tau = mesh.area(mesh.boundaries) .* mesh.tau;
oriented_normals = mesh.normals .* kron(mesh.orientation, [1;1]);
for e = 1:3
    Nx = oriented_normals(2*(e-1) + 1, :);
    Ny = oriented_normals(2*(e-1) + 2, :);
    block_e = (e-1)*d1 + (1:d1);
    area_tau_e = area_tau(e, :);

    %%% BOUNDARY MATRIX (FACE/FACE): PHYSICAL ELEMENT
    tau_ll = kron(area_tau_e, ll);
    tauLL(block_e, block_e, :) = reshape(tau_ll, [d1, d1, N_ELEMENTS]);

    %%% BOUNDARY MATRIX (VOLUME/VOLUME)
    wpp = (w.*p(block_e, :))'*(p(block_e, :)); % REFERENCE ELEMENT
    tauPP = tauPP + kron(area_tau_e, wpp); % PHYSICAL ELEMENT

    for mu = 1:2 % loop over permutations
        is_valid_edge = (mesh.permutation(e, :) == mu);

        block_mu = (mu - 1)*N_NODES + (1:N_NODES);
        %%% BOUNDARY MATRICES (FACE/VOLUME)
        % ON REFERENCE ELEMENT
        wlp = (w.*l(block_mu, :))'*p(block_e, :);

        % ON PHYSICAL ELEMENT
        nxmu = Nx .* is_valid_edge;
        nxLP(block_e, :) = nxLP(block_e, :) + kron(nxmu, wlp);

        nymu = Ny .* is_valid_edge;
        nyLP(block_e, :) = nyLP(block_e, :) + kron(nymu, wlp);

        taumu = area_tau(e, :) .* is_valid_edge;
        tauLP(block_e, :) = tauLP(block_e, :) + kron(taumu, wlp);
    end
end

% RESHAPING
nxLP = reshape(nxLP, [3*d1, d2, N_ELEMENTS]);
nyLP = reshape(nyLP, [3*d1, d2, N_ELEMENTS]);
tauLP = reshape(tauLP, [3*d1, d2, N_ELEMENTS]);
tauPP = reshape(tauPP, [d2, d2, N_ELEMENTS]);
end
