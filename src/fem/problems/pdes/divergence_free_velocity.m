function [u_star, divergence] = divergence_free_velocity(solution)
%-------------------------------------------------------------------------
%Description: The divergence-free postprocessed velocity u* in P_{k+1}
%             (postprocess_velocity) of a Stokes/Oseen solution on either
%             fluid formulation: the velocity gradient is field 1 of a
%             StokesProblem and is reconstructed from the PEERS stress and
%             rotation of an IncompressibleElasticityProblem
%             (gradient_from_stress).
%Output:    - u_star: {u1*; u2*} dofs in dubiner_basis(k+1), n x Nelts
%           - divergence: max_K ||div u*||_{L2(K)}
%-------------------------------------------------------------------------
arguments
    solution (1,1) FESolution
end
mesh = solution.Mesh;
k = solution.Order;
problem = solution.Problem;
fields = solution.Unknowns.field;
d2 = nchoosek(k + 2, 2);

u = fields(2).dof;
u_h = {u(1:d2, :); u(d2+1:end, :)};
d1 = k + 1;
u_hat = reshape(solution.Unknowns.hybrid(1).dof, 2*d1, []);   % per-face [x; y], flattened on the declarative path
u_hat = {u_hat(1:d1, :); u_hat(d1+1:end, :)};

% the velocity gradient: field 1 of the gradient form (TensorFEBasis, entity
% layout [L11; L21; L12; L22]) or reconstructed from the PEERS stress and
% the rotation of the stress form; on the declarative path the class is a
% plain FEProblem, so dispatch on the basis of field 1 and read the
% viscosity from the Data struct
if isa(problem, "StokesProblem") || isa(fields(1).basis, "TensorFEBasis")
    L = fields(1).dof;
    L_h = {L(1:d2, :),        L(2*d2+1:3*d2, :);
           L(d2+1:2*d2, :),   L(3*d2+1:4*d2, :)};
elseif isa(problem, "IncompressibleElasticityProblem") || isa(fields(1).basis, "PEERSBasis")
    assert(isfield(problem.Data, "viscosity"), "divergence_free_velocity:viscosity", ...
        "the stress reconstruction needs problem.Data.viscosity");
    L_h = gradient_from_stress(mesh, k, fields(1).dof, fields(3).dof, problem.Data.viscosity);
else
    error("divergence_free_velocity:problem", "no velocity gradient for a %s", class(problem));
end
[u_star, divergence] = postprocess_velocity(mesh, k, u_hat, u_h, L_h);
end
