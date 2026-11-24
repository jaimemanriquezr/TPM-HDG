function k_I = get_inverse_tensor(k, x, y)
%-------------------------------------------------------------------------
%Description: Evaluates the inverse of the 2×2 diffusion tensor κ at
%             quadrature points (x,y). Accepts κ as either a (2×2) cell
%             of function handles or a constant (2×2) double matrix.
%-------------------------------------------------------------------------
arguments (Input)
    k               % diffusion tensor κ: (2×2) cell of functions or (2×2) double
    x   (:,:)   double  % x-coordinates of quadrature points
    y   (:,:)   double  % y-coordinates of quadrature points
end
arguments (Output)
    k_I (2,2)   cell    % κ⁻¹ as (2×2) cell of double arrays matching size of x
end
if isa(k, 'cell')
    k_11 = k{1,1}(x, y); % dims: (p_gauss * cs_gauss) x n_of_patches
    k_21 = k{2,1}(x, y); % dims: (p_gauss * cs_gauss) x n_of_patches
    k_12 = k{1,2}(x, y); % dims: (p_gauss * cs_gauss) x n_of_patches
    k_22 = k{2,2}(x, y); % dims: (p_gauss * cs_gauss) x n_of_patches
elseif isa(k, 'double')
    k_11 = k(1,1) * ones(size(x));
    k_21 = k(2,1) * ones(size(x));
    k_12 = k(1,2) * ones(size(x));
    k_22 = k(2,2) * ones(size(x));
end

% computing kappa^{-1}
det_k = k_11.*k_22 - k_12.*k_21;
kI_11 = k_22./det_k; kI_12 = -k_12./det_k;
kI_21 = -k_21./det_k; kI_22 = k_11./det_k;

k_I = {kI_11, kI_12; kI_21, kI_22};
end
