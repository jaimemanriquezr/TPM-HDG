function mass = mass_matrix(mesh, varargin)
% MASS_MATRIX Element mass matrix, two call forms sharing one name.
%   mass_matrix(mesh, trial, test, ...)  - AbstractFEBasis form (DSL path)
%   mass_matrix(mesh, k, Vh_trial, Vh_test) - function-handle form
%                                          (low-level path, k = degree)
if nargin >= 2 && isa(varargin{1}, 'double')
    mass = handle_mass_matrix(mesh, varargin{:});
else
    mass = basis_mass_matrix(mesh, varargin{:});
end
end

function mass = basis_mass_matrix(mesh, trial, test, options)
    arguments
        mesh AbstractGrid
        trial AbstractFEBasis {mustBeSameElement(trial, mesh)} = DubinerBasis(0);
        test AbstractFEBasis {mustBeCompatible(test, trial)} = trial

        options.trial_operator FEOperator = IdentityOperator();
        options.test_operator FEOperator = IdentityOperator();
    end
    bases = [trial, test];
    [x_local, w] = get_quadrature(bases);
    p_trial = options.trial_operator.evaluate(trial, x_local, mesh);
    p_test = options.test_operator.evaluate(test, x_local, mesh);

    ops = [options.trial_operator, options.test_operator];
    if ~any([bases.IsPiolaTransformed, ops.is_spatial_dependent])
        p_trial = p_trial(:, :, :, :, 1);
        p_test = p_test(:, :, :, :, 1);
    end
    local_mass = pagemtimes(w .* p_test, 'transpose', p_trial, 'none');
    local_mass = reshape(sum(local_mass, 3:4), size(local_mass, 1), size(local_mass, 2), []);

    el_dim = test.GeometricInputDimension;
    switch el_dim
        case 2
            el_measure = reshape(mesh.volume, 1, 1, []);
        case 1
            el_measure = reshape(mesh.area, 1, 1, []);
    end
    mass = el_measure .* local_mass;
end

function M = handle_mass_matrix(mesh, k, Vh_trial, Vh_test)
arguments
    mesh AbstractGrid
    k double
    Vh_trial function_handle = @dubiner_basis
    Vh_test function_handle = Vh_trial
end

n_elements = size(mesh.elements, 2);
[xhat, w] = get_2d_reference_quadrature(k + 2, mesh.ElementType);

% Mass matrix computation
vol = mesh.volume;

p_trial = Vh_trial(k, xhat);
trial_dim = size(p_trial, 2);

p_test = Vh_test(k, xhat);
test_dim = size(p_test, 2);

mass = zeros(test_dim, n_elements*trial_dim);
for i = 1:length(w)
    p_dot = pagemtimes(p_trial(i, :, :, :), 'transpose', ...
                       p_test(i, :, :, :), 'none');
    mass = mass + kron(vol, w(i)*sum(p_dot, 3:4));
end
M = reshape(mass, [test_dim, trial_dim, n_elements]);
end
