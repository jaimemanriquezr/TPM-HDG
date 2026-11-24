function [lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(problem, mesh, k)
%-------------------------------------------------------------------------
%Description: Computes sparse assembly index arrays for the global HDG system.
%             Given lhs_entries and rhs_entries from get_matrices, the global
%             system is assembled as:
%                 H   = sparse(lhs_rows(:), lhs_cols(:), lhs_entries(:))
%                 RHS = accumarray(rhs_rows, rhs_entries(:))
%-------------------------------------------------------------------------
arguments (Input)
    problem (1,1)   FEProblem      % HDG problem object
    mesh            SimplicialGrid  % triangular mesh
    k       (1,1)   double          = problem.Order  % polynomial degree
end
arguments (Output)
    lhs_rows    (:,:,:) double  % row    index array (3*(k+1) x 3*(k+1) x Nelts)
    lhs_cols    (:,:,:) double  % column index array (3*(k+1) x 3*(k+1) x Nelts)
    rhs_rows    (:,1)   double  % RHS    index array (3*(k+1)*Nelts x 1)
end
d1 = k+1;
N_ELEMENTS = size(mesh.elements, 2);

% D.O.F. for edge "l" are (d1*(l-1), ..., d1*l) = d1*(l-1) + (1,...,d1)
face_DOF = (mesh.boundaries(:) - 1)*d1 + (1:d1);
% D.O.F. indexed by element (columns) and edge (row blocks)
face_DOF = reshape(face_DOF.', 3*d1, N_ELEMENTS);

[J, I] = meshgrid(1:3*d1);
lhs_rows = reshape(face_DOF(I(:), :), 3*d1, 3*d1, N_ELEMENTS);
lhs_cols = reshape(face_DOF(J(:), :), 3*d1, 3*d1, N_ELEMENTS);
rhs_rows = reshape(face_DOF, 3*d1*N_ELEMENTS, 1);
end
