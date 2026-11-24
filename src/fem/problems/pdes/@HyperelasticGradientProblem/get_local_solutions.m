function [local_solutions, trace_solutions] = get_local_solutions(obj, x_hat, local_operators, mesh) %#ok<INUSL>
% Recovers (H, u) per element from A1 [H; u] = A2 û_K + Af.
d1 = size(local_operators{2}, 2)/6;
d2 = d1*(d1 + 1)/2;
N_FACES = size(mesh.faces, 2);
N_ELEMENTS = size(mesh.elements, 2);
N_DOF = d1*N_FACES;
u_hat = [reshape(x_hat(1:N_DOF), d1, N_FACES); reshape(x_hat(N_DOF + 1:2*N_DOF), d1, N_FACES)];
u_K = element_traces(mesh, u_hat, d1);
[A1, A2, Af] = deal(local_operators{:});
sol = reshape(pagemldivide(A1, pagemtimes(A2, u_K) + Af), [], N_ELEMENTS);
local_solutions = {sol(1:4*d2, :), {sol(4*d2 + (1:d2), :); sol(5*d2 + (1:d2), :)}};
trace_solutions = {{u_hat(1:d1, :); u_hat(d1 + 1:end, :)}};
end
