function [lhs_rows, lhs_cols, rhs_rows] = get_assemble_indices(obj, mesh, k)
%-------------------------------------------------------------------------
%Description: Computes sparse assembly index arrays for the global HDG system.
%             Global hybrid unknowns are [û_x; û_y; ρ]: vector velocity traces
%             (2*(k+1)*Nfaces DOFs) and per-element pressure means (Nelts DOFs).
%             Given lhs_entries and rhs_entries from get_matrices:
%                 H   = sparse(lhs_rows(:), lhs_cols(:), lhs_entries(:))
%                 RHS = accumarray(rhs_rows, rhs_entries(:))
%-------------------------------------------------------------------------
arguments (Input)
    obj     (1,1)   FEProblem      % HDG problem object
    mesh            SimplicialGrid  % triangular mesh
    k       (1,1)   double          = obj.Order  % polynomial degree
end
arguments (Output)
    lhs_rows    (:,:,:) double  % row    index array (2*3*(k+1)+1 x 2*3*(k+1)+1 x Nelts)
    lhs_cols    (:,:,:) double  % column index array (2*3*(k+1)+1 x 2*3*(k+1)+1 x Nelts)
    rhs_rows    (:,1)   double  % RHS    index array ((2*3*(k+1)+1)*Nelts x 1)
end
d1 = k+1;
N_DOF = d1*size(mesh.faces, 2); % # of DOF belonging to 1 scalar variable
N_ELEMENTS = size(mesh.elements, 2);

face_DOF = (mesh.boundaries(:) - 1)*d1 + (1:d1);
face_DOF = reshape(face_DOF.', 3*d1, N_ELEMENTS);
element_DOF = 1:N_ELEMENTS;
all_DOF = [face_DOF;
           face_DOF + N_DOF;
           element_DOF + 2*N_DOF];
[J, I] = meshgrid(1:(2*3*d1 + 1));
lhs_rows = reshape(all_DOF(I(:), :), 2*3*d1 + 1, 2*3*d1 + 1, N_ELEMENTS);
lhs_cols = reshape(all_DOF(J(:), :), 2*3*d1 + 1, 2*3*d1 + 1, N_ELEMENTS);
rhs_rows = reshape(all_DOF, (2*3*d1 + 1)*N_ELEMENTS, 1);
end
