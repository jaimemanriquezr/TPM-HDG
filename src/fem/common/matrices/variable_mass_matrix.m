function M = variable_mass_matrix(mesh, k, m, options)
arguments
    mesh SimplicialGrid
    k double
    m cell = {@(x,y) 1*(x.*y).^0}
    options.basis function_handle = @dubiner_basis
end
%-------------------------------------------------------------------------
%Description: Computes the variable mass matrix
%               M_ij = int_{K} m(x)*P_i*P_j
%Input:     - T: SimplicialGrid
%           - k: polynomial degree
%           - m: vectorized 2D function
%Output:    - M: element-wise mass matrix (d2 x d2 x Nelts)
%-------------------------------------------------------------------------
N_ELEMENTS = size(mesh.elements, 2);
[xy, rs, w] = quadpoints_elements_2d(mesh, k, output="cell");
p = options.basis(k, rs{:});
d2 = size(p, 2);

% MASS MATRIX
M = cell(size(m));
for i = 1:size(m, 1)
    for j = 1:size(m, 2)
        mass = zeros(d2, d2*N_ELEMENTS);
        m_ij = (mesh.volume).*m{i,j}(xy{:});
        for q = 1:length(w)
            mass = mass + kron(m_ij(q,:), w(q)*p(q,:)'*p(q,:));
        end
        M{i, j} = reshape(mass, [d2, d2, N_ELEMENTS]);
    end
end
end