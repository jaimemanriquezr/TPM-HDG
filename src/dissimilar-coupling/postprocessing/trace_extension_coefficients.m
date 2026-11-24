function ext = trace_extension_coefficients(problem, mesh, partition, sigma_dof, gamma_dof, options)
%-------------------------------------------------------------------------
%Description: Face-basis coefficients of the TPM path extension
%
%               Ext(u) := int_0^1 grad u(x(xi)) . m dxi = T_h(u_hat) - u_hat,
%
%             the correction transfer_trace_to_interface adds to the
%             skeleton trace to reach Sigma, returned as Legendre
%             coefficients on the interface faces so it can be added to a
%             right-hand side.
%
%             The T2 kinematic row compares a path-extended trace at step
%             n+1 with the raw face trace of the history h0 at step n. On a
%             retreated mesh this leaves a (1/cdt)*Ext term, an O(delta/dt)
%             datum error that does not vanish under refinement. Adding Ext
%             to h0 transports the history to Sigma and removes it.
%
%             grad u = c1 sigma - c2 tr(sigma) I + gamma is reconstructed
%             as in transfer_trace_to_interface, then L2-projected face by
%             face onto the skeleton's Legendre basis.
%
%             Takes dofs, not a solution: the T2 history needs the Newmark
%             combination u^0 + c_v v^0 + c_a a^0, which the caller forms
%             in (sigma, gamma) since Ext is linear in them.
%
%See also: transfer_trace_to_interface, FluidStructureCoupling.solve
%-------------------------------------------------------------------------
arguments (Input)
    problem     (1,1)   LinearElasticityProblem
    mesh                SimplicialGrid
    partition   (1,1)   struct
    sigma_dof   (:,:)   double  % PEERS stress dofs of the field to extend
    gamma_dof   (:,:)   double  % rotation dofs of the same field
    options.Order double {mustBeScalarOrEmpty} = []
end
arguments (Output)
    ext (:,:) double    % (2*d1) x n_patches, matching solve.m's h0 layout
end
k = options.Order;
if isempty(k), k = problem.Order; end
d1 = k + 1;

[lambda, mu] = LinearElasticityProblem.get_lame_parameters( ...
    problem.Data.young_modulus, problem.Data.poisson_ratio);
c1 = 1/(2*mu);
c2 = lambda/(2*mu*(2*lambda + 2*mu));

s_quad  = size(partition.quadrature.local, 2);
[~, wS]  = lgwt(s_quad, 0, 1);

n_of_patches = size(partition.faces, 3);
% One column per patch; a face may carry several patches. The solve
% expands h0 to one column per patch and integrates column e over patch
% e's quadrature points only, so each fit below uses patch e's own points.

G = get_affine_map(mesh, partition.elements, inverse = true);

ext = zeros(2*d1, n_of_patches);
for e = 1:n_of_patches
    K = partition.elements(e);

    path_integral = path_gradient_integral(k, partition, e, G(:, :, e), ...
        sigma_dof(:, K), gamma_dof(:, K), c1, c2);                % s x 2

    % L2 projection onto the face's Legendre basis. The patch length is a
    % common factor of both sides and cancels, so it is omitted.
    l = legendre_basis(k, partition.quadrature.local(1, :, e));    % s x d1
    if isfield(partition, "collapsed") && partition.collapsed(e) && partition.area(e) == 0
        % the single point paired with a collapsed face (get_partitions):
        % every node is the same point, so fit the minimum-norm
        % coefficients that reproduce the extension there
        coeffs = pinv(l(1, :)) * path_integral(1, :);              % d1 x 2
    else
        coeffs = (l' * (wS .* l)) \ (l' * (wS .* path_integral));  % d1 x 2
    end
    ext(:, e) = [coeffs(:, 1); coeffs(:, 2)];
end
end
