function coupled = fsi_dsl_coupled_problem(k, fluid_mesh, solid_mesh, interface, material, pd, interface_sources, beta)
% The no-slip steady FSI problem (stress-velocity PEERS Stokes/Oseen fluid,
% PEERS linear elasticity solid) written in the form language, coupled by
% TPM on the InterfaceProblem path. Test functions stand on the left:
% vectors pair as v'*u, tensors as sum(T.*S), scalars stand outside.
%
%   material: struct with viscosity, young_modulus, poisson_ratio, time_step
%   pd, interface_sources: from FluidStructureCoupling.load_problems
%   beta: convective velocity (empty = Stokes). An element-wise vector
%       field (Datum, or an FEField such as the previous Picard velocity);
%       adds the Oseen terms of Solano-Vargas (2022, eq. 3.4) to the
%       momentum row and the numerical flux:
%           -(u x beta, grad v)_K + <(beta.n) u_hat, v>_dK
%           -<(beta.n) u_hat, v_hat>_dK
%       The fluid stabilization is mu*(mesh tau): with convection the
%       caller puts tau = 1 + max|beta.n|/(2 mu) on fluid_mesh
%       (add_stabilization), so the transmission rows, which read mesh.tau,
%       match the local problem.
arguments
    k
    fluid_mesh
    solid_mesh
    interface
    material struct
    pd struct
    interface_sources cell
    beta = []
end
J  = [0, 1; -1, 0];
Pk = DubinerBasis(k);
dx = Measure(preset = "dx");
dS = Measure(preset = "dS");
n_ = ExteriorNormal();
tau_ = Stabilization();

% --- fluid: incompressible elasticity (PEERS stress / velocity) ---------
mu = material.viscosity;
c1 = 1/(2*mu);  c2 = 1/(4*mu);
fluid = FEProblem(Order = k, ...
    Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
              FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
fluid.Mesh = fluid_mesh;
fluid.Data.viscosity = mu;   % read by divergence_free_velocity
[sg, u, gm, xi, u_hat, rho] = get_trial_functions(fluid.Space);
[T, v, e, z, v_hat, s_rho]  = get_test_functions(fluid.Space);
local = ...
      c1*sum(T.*sg)*dx - c2*trace(T)*trace(sg)*dx ...
    + divergence(T)'*u*dx + sum(T.*(gm*J))*dx + trace(T)*xi*dx ...
    - sum(T.*(u_hat*n_.'))*dS ...
    - v'*divergence(sg)*dx + mu*v'*(tau_*u)*dS ...
    - mu*v'*(tau_*u_hat)*dS + v'*Datum(pd.f_F)*dx ...
    + sum((e*J).*sg)*dx ...
    + z*trace(sg)*dx - z*rho*dx;
flux = ...
      v_hat'*(sg*n_)*dS - mu*v_hat'*(tau_*u)*dS ...
    + mu*v_hat'*(tau_*u_hat)*dS ...
    + s_rho*(u_hat'*n_)*dS;
if ~isempty(beta)
    bn_u_hat = dot(beta, n_)*u_hat;
    local = local - sum(gradient(v).*(u*beta.'))*dx + dot(bn_u_hat, v)*dS;
    flux  = flux - dot(bn_u_hat, v_hat)*dS;
end
fluid.Operator = AffineOperator(local) + flux;
fluid.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_F), ...
    "Dirichlet", projection_order = k);

% --- solid: linear elasticity (PEERS stress / displacement) -------------
E_S = material.young_modulus;  nu_S = material.poisson_ratio;
lambda = E_S*nu_S/((1 + nu_S)*(1 - 2*nu_S));
mu_S   = E_S/(2*(1 + nu_S));
c1_S = 1/(2*mu_S);  c2_S = lambda/(4*mu_S*(lambda + mu_S));
solid = FEProblem(Order = k, ...
    Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
              FESpace(LegendreBasis(k).to_vector(2))]);
solid.Mesh = solid_mesh;
[sg, u, gm, u_hat] = get_trial_functions(solid.Space);
[T, v, e, v_hat]   = get_test_functions(solid.Space);
solid.Operator = AffineOperator( ...
      c1_S*sum(T.*sg)*dx - c2_S*trace(T)*trace(sg)*dx ...
    + divergence(T)'*u*dx + sum(T.*(gm*J))*dx - sum(T.*(u_hat*n_.'))*dS ...
    - v'*divergence(sg)*dx + v'*(tau_*u)*dS ...
    - v'*(tau_*u_hat)*dS + v'*Datum(pd.f_S)*dx ...
    + sum((e*J).*sg)*dx);
solid.Operator = solid.Operator + ( ...
      v_hat'*(sg*n_)*dS - v_hat'*(tau_*u)*dS + v_hat'*(tau_*u_hat)*dS);
if ~isempty(get_face_by_tag(solid_mesh, "Dirichlet"))
    solid.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_S), ...
        "Dirichlet", projection_order = k);
else
    solid.BoundaryConditions = StrongEssentialBoundaryConditions();   % enclosed solid
end

% --- coupling -----------------------------------------------------------
coupled = InterfaceProblem.fsi(fluid, solid, interface, ...
    viscosity = mu, young_modulus = E_S, poisson_ratio = nu_S, ...
    domain_measure = pd.measure_F, ...
    scheme = "steady", time_step = material.time_step, pressure_mean = "pinned", ...
    source_neumann = interface_sources{1}, ...
    source_dirichlet = interface_sources{2});
end
