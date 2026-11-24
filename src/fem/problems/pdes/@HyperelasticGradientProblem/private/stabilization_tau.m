function [tau, marked, constants] = stabilization_tau(obj, mesh, k, H, constants)
% The face stabilization τ (3 x N, local face e of element K) of
% Data.stabilization (see HyperelasticGradientProblem); H are the dofs
% of the iterate's gradient H_h (4*d2 x N, [] = 0), read by "adaptive"
% and "cockburn_shen". marked (1 x N) flags the elements given
% C_KF,o/h_F under "cockburn_shen" (false otherwise). constants are the
% C_KF,o of cockburn_shen_constants (optional, [] = computed here when
% needed), returned for the next call on the same mesh and k.
if nargin < 5, constants = []; end
data = obj.Data;
marked = false(1, size(mesh.elements, 2));
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(data.young_modulus, data.poisson_ratio);
N_ELEMENTS = size(mesh.elements, 2);
switch data.stabilization
    case "constant"
        tau = data.stabilization_constant*mu/data.characteristic_length*ones(3, N_ELEMENTS);
    case "mesh_size"
        tau = data.stabilization_constant*mu./mesh.area(mesh.boundaries);
    case "adaptive"
        if isempty(H) || ~any(H, 'all')
            % the largest eigenvalue of Hooke's law on full gradients
            tau = (2*mu + 2*lambda)/(2*data.characteristic_length)*ones(3, N_ELEMENTS);
            return
        end
        d2 = nchoosek(k + 2, 2);
        s = get_1d_reference_quadrature(2*k + 2);
        s = s(:); nf = numel(s);
        p = dubiner_basis(k, [0*s; 1 - s; 1 - s], [1 - s; 0*s; s]);   % the three faces
        G = pagemtimes(p, reshape(H, d2, 4, N_ELEMENTS));               % 3nf x 4 x N
        G = reshape(permute(G, [1 3 2]), [], 4);
        [~, dP] = HyperelasticGradientProblem.first_piola(G, lambda, mu, data.material);
        D = permute(dP, [2 3 1]);
        rho = max(pageeig((D + pagetranspose(D))/2), [], 1);
        rho = reshape(max(reshape(rho, nf, 3, N_ELEMENTS), [], 1), 3, N_ELEMENTS);
        rho_face = accumarray(mesh.boundaries(:), rho(:), [size(mesh.faces, 2), 1], @max);
        tau = reshape(rho_face(mesh.boundaries), 3, N_ELEMENTS)/(2*data.characteristic_length);
    case "cockburn_shen"
        % Cockburn–Shen: τ_KF = C_KF,o/h_F on the elements of T_h^(-,nb)
        % (the elements where A(F_h) = ∂P/∂F has a negative eigenvalue at
        % a quadrature point, and every element sharing a vertex with
        % one), τ = 1 elsewhere; both times stabilization_scale
        tau = ones(3, N_ELEMENTS);
        if ~isempty(H) && any(H, 'all')
            d2 = nchoosek(k + 2, 2);
            [xq, ~] = get_2d_reference_quadrature(2*k + 2);
            p = dubiner_basis(k, xq(1, :).', xq(2, :).');
            G = pagemtimes(p, reshape(H, d2, 4, N_ELEMENTS));
            n_q = size(p, 1);
            G = reshape(permute(G, [1 3 2]), [], 4);
            [~, dP] = HyperelasticGradientProblem.first_piola(G, lambda, mu, data.material);
            D = permute(dP, [2 3 1]);
            lambda_min = reshape(min(pageeig((D + pagetranspose(D))/2), [], 1), n_q, N_ELEMENTS);
            negative = any(lambda_min < -1e-10*(2*mu + 2*lambda), 1);
            touched = false(1, size(mesh.coordinates, 2));
            touched(mesh.elements(:, negative)) = true;
            marked = any(touched(mesh.elements), 1);
            if isempty(constants), constants = cockburn_shen_constants(mesh, k); end
            h_F = mesh.area(mesh.boundaries);
            tau(:, marked) = constants(:, marked)./h_F(:, marked);
        end
        tau = data.stabilization_scale*tau;
end
end
