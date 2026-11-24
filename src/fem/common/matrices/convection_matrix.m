function varargout = convection_matrix(mesh, varargin)
% CONVECTION_MATRIX Element convection matrices, two call forms.
%   convection = convection_matrix(mesh, trial, test, ...)  - AbstractFEBasis
%       form (DSL path); returns a struct with fields x, y.
%   [Cx, Cy] = convection_matrix(mesh, k, basis=...)  - function-handle
%       form (low-level path, k = degree).
if nargin >= 2 && isa(varargin{1}, 'double')
    [varargout{1:max(nargout, 1)}] = handle_convection_matrix(mesh, varargin{:});
else
    varargout{1} = basis_convection_matrix(mesh, varargin{:});
end
end

function convection = basis_convection_matrix(mesh, trial, test, options)
    arguments
        mesh SimplicialGrid
        trial AbstractFEBasis {mustBeSameElement(trial, mesh)} = DubinerBasis(0);
        test AbstractFEBasis {mustBeCompatible(test, trial)} = trial

        options.test_operator FEOperator = IdentityOperator();
    end
    bases = [trial, test];
    [x_local, w] = get_quadrature(bases);
    [p_trial_x, p_trial_y] = trial.evaluate_derivatives(1, x_local);
    p_test = options.test_operator.evaluate(test, x_local, mesh);

    % [ASSUMPTION] only scalar bases
    p_trial_x = squeeze(p_trial_x);
    p_trial_y = squeeze(p_trial_y);
    p_test = squeeze(p_test);

    ops = [options.test_operator];
    if ~any([bases.IsPiolaTransformed, ops.is_spatial_dependent])
        w_p_test = w .* p_test(:, :, 1);
        local_conv_x = 1/2 * w_p_test' * p_trial_x(:, :, 1);
        local_conv_y = 1/2 * w_p_test' * p_trial_y(:, :, 1);
    else
        w_p_test = w .* p_test;
        local_conv_x = 1/2 * pagemtimes(w_p_test, 'transpose', p_trial_x, 'none');
        local_conv_y = 1/2 * pagemtimes(w_p_test, 'transpose', p_trial_y, 'none');
    end

    % detB*B^-T
    B = get_affine_map(mesh, output="array");
    axx = B(2, 2, :); axy = -B(2, 1, :);
    ayx = -B(1, 2, :); ayy = B(1, 1, :);

    convection.x = axx.*local_conv_x + axy.*local_conv_y;
    convection.y = ayx.*local_conv_x + ayy.*local_conv_y;
end

function [Cx, Cy] = handle_convection_matrix(mesh, k, options)
arguments
    mesh SimplicialGrid
    k double
    options.basis function_handle = @dubiner_basis
end
%-------------------------------------------------------------------------
%Description: Computes the convection matrices
%             Cx_ij = int_{K} P_i*Px_j
%             Cy_ij = int_{K} P_i*Py_j
%Input:     - T: SimplicialGrid
%           - k: polynomial degree
%Output:    - Cx: x-direction convection matrix (d2 x d2 x Nelts)
%           - Cy: y-direction convection matrix (d2 x d2 x Nelts)
%-------------------------------------------------------------------------

% QUADRATURE
N_ELEMENTS = size(mesh.elements, 2);
% QUADRATURE
[~, rs, w, B] = quadpoints_elements_2d(mesh, k, output="cell");

% COLLOCATION
[p, px, py] = options.basis(k, rs{:});
d2 = size(p,2);

% CONVECTION MATRICES (reference element)
wp = w .* p;
convhatx = 1/2*wp'*px;
convhaty = 1/2*wp'*py;

% detB*B^-T
axx = B{2, 2}(:); axy = -B{2, 1}(:);
ayx = -B{1, 2}(:); ayy = B{1,1}(:);

% CONVECTION MATRICES (physical element)
convx = kron(axx', convhatx) + kron(axy', convhaty);
convy = kron(ayx', convhatx) + kron(ayy', convhaty);

% RESHAPING
Cx = reshape(convx, [d2, d2, N_ELEMENTS]);
Cy = reshape(convy, [d2, d2, N_ELEMENTS]);
end
