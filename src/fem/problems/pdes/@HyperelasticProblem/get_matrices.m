function [local_operators, lhs_entries, rhs_entries] = get_matrices(pde, mesh, polynomial_degree)
%-------------------------------------------------------------------------
%Description: Newton step of the hyperelastic HDG scheme at the iterate
%             Data.iterate = (u⁰, û⁰) (a struct with fields displacement,
%             2*d2 x N, and displacement_skeleton, 2*d1 x F, the layout
%             of state.solid), written for the new iterate (σ, u, γ, û) so
%             that it has the shape of the linear class. The local
%             residual is
%               A1 X − A2 û − Af − [m(L_h) in the constitutive and symmetry rows],
%             m = [(ζ⁻¹N(L_h), τ); (N(L_h), η)] the moments of the stress
%             remainder N(L) = P(I + L) − Cε(L) at the discrete gradient
%             L_h = R_u u + R_û û (nonlinear_terms, gradient_operators),
%             with the Kelvin–Voigt stress at the Newmark rate of L_h when
%             Data.solid_viscosity > 0 (viscous_rate);
%             its linearization, in those rows,
%               (A1 − Dm R_u) X = (A2 + Dm R_û) û + Af + m⁰ − Dm L⁰,
%             replaces A1, A2, Af, while the momentum rows and the trace
%             rows are those of LinearElasticityProblem.get_matrices
%             (inertial terms included). Returns the local operators (so
%             get_local_solutions recovers the new iterate) and the
%             statically condensed hybrid system. Without an iterate the
%             outputs are the linear ones. The iterate-independent
%             operators come from linear_operators (cached by solve in
%             Data.newton_operators during the Newton loop).
%-------------------------------------------------------------------------
arguments (Input)
    pde               (1,1)   HyperelasticProblem
    mesh              SimplicialGrid
    polynomial_degree (1,1)   double = pde.Order
end
k = polynomial_degree;
iterate = [];
if isfield(pde.Data, "iterate"), iterate = pde.Data.iterate; end
cached = isfield(pde.Data, "newton_operators") && ~isempty(pde.Data.newton_operators);
if isempty(iterate) && ~cached
    [local_operators, lhs_entries, rhs_entries] = get_matrices@LinearElasticityProblem(pde, mesh, k);
    return
end
operators = linear_operators(pde, mesh, k);
[A1, A2, Af] = deal(operators.local_operators{:});
J = A1;
if ~isempty(iterate)
    N_ELEMENTS = size(mesh.elements, 2);
    d1 = k + 1;
    d2 = nchoosek(k + 2, 2);
    d2peers = size(A1, 1) - 3*d2;
    ix_u = d2peers + (1:2*d2);
    ix_c = [1:d2peers, d2peers + 2*d2 + (1:d2)];   % constitutive and symmetry rows

    u0 = reshape(iterate.displacement, 2*d2, 1, N_ELEMENTS);
    u0_K = element_traces(mesh, iterate.displacement_skeleton, d1);
    l0 = pagemtimes(operators.R_u, u0) + pagemtimes(operators.R_hat, u0_K);
    [m, Dm] = nonlinear_terms(pde, mesh, k, l0, operators.viscous);

    J(ix_c, ix_u, :) = J(ix_c, ix_u, :) - pagemtimes(Dm, operators.R_u);
    A2(ix_c, :, :) = A2(ix_c, :, :) + pagemtimes(Dm, operators.R_hat);
    Af(ix_c, :, :) = Af(ix_c, :, :) + m - pagemtimes(Dm, l0);
end

local_operators = {J, A2, Af};
X = pagemldivide(J, cat(2, A2, Af));
lhs_entries = pagemtimes(operators.A3, X(:, 1:end-1, :)) + operators.A4;
rhs_entries = -pagemtimes(operators.A3, X(:, end, :));
end
