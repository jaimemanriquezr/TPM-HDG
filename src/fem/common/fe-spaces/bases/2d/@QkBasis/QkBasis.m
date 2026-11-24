classdef QkBasis < ScalarFEBasis
    % Discontinuous tensor-product basis of degree up to k in each variable
    % on the reference square [0, 1]^2, Dimension = (k+1)^2.
    %
    %   basis = QkBasis(k)
    %
    % The functions are products of the orthonormal Legendre polynomials
    % L_i of legendre_eval: phi_(i,j)(r, s) = L_i(r) * L_j(s), with
    % column index i + (k+1)*(j-1) (r fastest). With this ordering the
    % reference mass matrix is kron(I, I) = I, and on an axis-aligned
    % rectangle every element matrix of a separable operator is a
    % Kronecker product kron(A_s, A_r).
    %
    % [P, dPx, dPy, dPxx, dPyy, dPxy] = QkBasis.evaluate_basis(k, r, s)
    % returns the values and first and second derivatives at (r, s).
    %
    % Supported on CartesianGrid; LegendreBasis(k) is the matching trace
    % space on its faces.
    %
    % See also: ScalarFEBasis, DubinerBasis, LegendreBasis, CartesianGrid, legendre_eval
    properties (Constant, Hidden)
        ElementType = "Quad"
        GeometricInputDimension = 2
        GeometricOutputDimension = 1
        IsDiscontinuous = true
        IsPiolaTransformed = false
    end

    methods
        function obj = QkBasis(k)
            arguments
                k double = [];
            end
            obj.Order = k;
        end

        function E = compute_l2_error(obj, dof, mesh, options)
            % L2 norm of (dof field - exact_solution) on a CartesianGrid;
            % dof is Dimension-by-n_el (as in DubinerBasis.compute_l2_error)
            arguments
                obj QkBasis
                dof double
                mesh CartesianGrid
                options.exact_solution function_handle = @(x,y) 0*(x.*y);
            end
            [xy, rs, w] = quadpoints_elements_2d(mesh, obj.Order + 1, output="array");
            P = obj.evaluate(rs) * dof;
            P_ex = reshape(options.exact_solution(xy(1,:,:), xy(2,:,:)), size(P));
            E = sqrt(sum(mesh.volume .* (w'*(P - P_ex).^2)));
        end
    end

    methods (Access = protected)
        function dim_space = return_Dimension(obj)
            if isempty(obj.Order)
                dim_space = 0;
            else
                dim_space = (obj.Order + 1)^2;
            end
        end
    end

    methods (Static)
        function [P, dPx, dPy, dPxx, dPyy, dPxy] = evaluate_basis(k, r, s)
            arguments
                k (1, 1) double
                r double
                s double = r(2, :);
            end
            if nargin < 3
                r = r(1, :);
            end
            [Lr, dLr, ddLr] = legendre_1d(k, r(:));
            [Ls, dLs, ddLs] = legendre_1d(k, s(:));

            % row-wise Kronecker product: column i + (k+1)*(j-1) is a(:,i).*b(:,j)
            tp = @(a, b) reshape(a .* permute(b, [1, 3, 2]), size(a, 1), []);
            P = tp(Lr, Ls);
            dPx = tp(dLr, Ls);
            dPy = tp(Lr, dLs);
            dPxx = tp(ddLr, Ls);
            dPyy = tp(Lr, ddLs);
            dPxy = tp(dLr, dLs);
        end
    end
end

function [L, dL, ddL] = legendre_1d(k, s)
% The k+1 Legendre polynomials orthonormal on [0, 1] (the normalization of
% legendre_eval), with their first and second derivatives, at the column s.
t = 2*s - 1;
[L, dL, ddL] = deal(zeros(numel(s), k + 1));
for n = 0:k
    [p, dp, ddp] = jacobi_eval(n, 0, 0, t);
    c = sqrt(2*n + 1);
    L(:, n+1) = c*p;
    dL(:, n+1) = 2*c*dp;
    ddL(:, n+1) = 4*c*ddp;
end
end
