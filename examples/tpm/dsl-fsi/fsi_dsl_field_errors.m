function E = fsi_dsl_field_errors(solution, mesh, sigma, u, gamma)
% L2 errors of sigma_h, u_h, gamma_h and u_hat, and ||P_M u - u_hat||_h
% (P_M the L2 projection onto the trace space). On face fields
% compute_l2_error is h-weighted (sum_e h_e ||.||_e^2).
fields = solution.Unknowns.field;
u_hat = solution.Unknowns.hybrid(1);
P_u = Datum(u).project(mesh, u_hat.basis);
E = [compute_l2_error(fields(1), mesh, Datum(sigma));
     compute_l2_error(fields(2), mesh, Datum(u));
     compute_l2_error(fields(3), mesh, Datum(gamma));
     compute_l2_error(u_hat, mesh, Datum(u));
     compute_l2_error(FEField(FESpace(u_hat.basis, mesh), u_hat.dof - P_u.dof), mesh)];
end
