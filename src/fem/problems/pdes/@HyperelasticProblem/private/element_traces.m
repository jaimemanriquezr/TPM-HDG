function u_K = element_traces(mesh, u_hat, d1)
% the trace dofs of each element's three faces, as get_local_solutions
% gathers them: 6*d1 x 1 x N (u_hat is 2*d1 x F, x block over y block)
N_ELEMENTS = size(mesh.elements, 2);
ux_hat = u_hat(1:d1, :);
uy_hat = u_hat(d1+1:end, :);
u_K = reshape([reshape(ux_hat(:, mesh.boundaries), [3*d1, N_ELEMENTS]);
               reshape(uy_hat(:, mesh.boundaries), [3*d1, N_ELEMENTS])], [], 1, N_ELEMENTS);
end
