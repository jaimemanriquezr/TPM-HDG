function operators = linear_operators(pde, mesh, k)
% The operators of the Newton step that do not depend on the iterate, at
% mesh.tau: the linear local operators {A1, A2, Af} of
% LinearElasticityProblem.get_matrices (inertial terms and Data.state
% included), the gradient operators R_u, R_hat (gradient_operators) and
% the trace operators A3, A4 (trace_operator) and, with a solid
% viscosity, the data of the strain rate (viscous_rate; struct.empty
% without one). solve keeps them in
% Data.newton_operators during its Newton loop; they are taken from there
% when built for the current mesh.tau.
if isfield(pde.Data, "newton_operators") && ~isempty(pde.Data.newton_operators) ...
        && isequal(pde.Data.newton_operators.tau, mesh.tau)
    operators = pde.Data.newton_operators;
    return
end
pde.Data.iterate = struct.empty;
pde.Data.newton_operators = [];
local_operators = get_matrices(pde, mesh, k);
[R_u, R_hat] = gradient_operators(mesh, k, local_operators{1}, local_operators{2});
d2 = nchoosek(k + 2, 2);
[A3, A4] = trace_operator(mesh, k, size(local_operators{1}, 1) - 7*d2);
operators = struct("tau", mesh.tau, "local_operators", {local_operators}, ...
    "R_u", R_u, "R_hat", R_hat, "A3", A3, "A4", A4, ...
    "viscous", viscous_rate(pde, mesh, k, R_u, R_hat));
end
