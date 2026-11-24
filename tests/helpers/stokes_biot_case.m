function out = stokes_biot_case(geometry, k, h, options)
%STOKES_BIOT_CASE  One manufactured Stokes-Biot run for the tests.
%
%   out = stokes_biot_case(geometry, k, h, ...)
%
%   geometry  "flat": unit square, Sigma: y = 1/2, fluid above, Biot below,
%             fitted_split with sizes (0.8 h, h) (fluid, Biot);
%             "circle": centre (1/2, 1/2), R = 1/4, fluid outside, Biot
%             inside, interface_meshes(curve, "fitted", 0.8 h, gap = gap,
%             h_solid = h).
%   The Biot mesh is the coarser one at Sigma (StokesBiotCoupling
%   .interface_problem warns when it is finer).
%   Every outer boundary is Dirichlet (u_F for the fluid; eta and p for
%   the Biot side). The exact solution is
%   StokesBiotCoupling.get_manufactured_solution with nu_p the outward
%   normal of the Biot domain, and its options are the name-value options
%   below that carry the same names (defaults as there).
%
%   Other options:
%     gap                  retreat of the fitted meshes ([] = h^2)
%     time_scheme          "euler" | "bdf2" (the Biot Data.time_scheme)
%     time_step            Delta t
%     n_steps              number of steps (1)
%     rate                 "from_rest" (default): one step from the zero
%                          state, d_t w = w / time_step; "exact": start
%                          from initial_state(biot_problem, biot_mesh, k)
%                          at t = 0 and march n_steps with the exact d_t
%     stabilization_darcy  tau_p / (||kappa|| tau)
%     stabilization_solid  tau_s / tau (tau the Biot mesh.tau)
%     pore_pressure_mean   the mean of p over the Biot mesh, used only
%                          when the pore-pressure level is free
%                          (InterfaceProblem.stokes_biot); "exact" = the
%                          mean of the exact p_P over the Biot mesh
%
%   out.errors     L2 errors at the last step, fields sigma_F, u_F,
%                  sigma_P, u_P, q_P, p_P
%   out.h          h
%   out.residue    ||H x - b|| / ||b|| on the free dofs of the last solve
%   out.solutions  struct: fluid, biot (the FESolutions of the last
%                  step), state, global_system (of the last step), meshes
%                  {fluid_mesh, biot_mesh}, data (mms.at at the last time)
%
%   Only the tests use it.
arguments
    geometry (1,1) string {mustBeMember(geometry, ["flat", "circle"])}
    k        (1,1) double
    h        (1,1) double {mustBePositive}
    options.stream_function (1,1) string = "sin(pi*x)*sin(pi*y)"
    options.pressure_fluid  (1,1) string = "sin(pi*x)*exp(y - 0.5)"
    options.displacement    (2,1) string = ["x*sin(y)"; "y*cos(x)"]
    options.pressure_biot   (1,1) string = "cos(pi*x)*sin(pi*y)"
    options.viscosity       (1,1) double = 1.0
    options.young_modulus   (1,1) double = 1.0
    options.poisson_ratio   (1,1) double = 0.25
    options.biot_alpha      (1,1) double = 1.0
    options.storage         (1,1) double = 0.0
    options.diffusion       (2,2) double = eye(2)
    options.hooke_modulus   (1,1) double = 0.0
    options.friction_bjs    (1,1) double = 0.0
    options.time_step       (1,1) double = 1.0
    options.rate (1,1) string {mustBeMember(options.rate, ["exact", "from_rest"])} = "from_rest"
    options.gap double {mustBeScalarOrEmpty} = []
    options.time_scheme (1,1) string {mustBeMember(options.time_scheme, ["euler", "bdf2"])} = "euler"
    options.n_steps (1,1) double {mustBeInteger, mustBePositive} = 1
    options.stabilization_darcy (1,1) double {mustBePositive} = 1.0
    options.stabilization_solid (1,1) double {mustBePositive} = 1.0
    options.pore_pressure_mean = 0.0
end
assert(options.rate == "exact" || options.n_steps == 1, ...
    "stokes_biot_case:fromRest", ...
    "rate = ""from_rest"" is one step from the zero state; use ""exact"" to march");
gap = options.gap;  if isempty(gap), gap = h^2; end
dt = options.time_step;

%% geometry: (fluid, Biot) meshes and the outward normal of the Biot domain
switch geometry
    case "flat"
        curve = FlatCurveHorizontal(0.5);
        [fluid_mesh, biot_mesh] = fitted_split(curve, [0.8*h, h], gap);
        normal_p = ["0"; "1"];
    case "circle"
        curve = CircularCurve([0.5; 0.5], 0.25);
        [fluid_mesh, biot_mesh] = interface_meshes(curve, "fitted", 0.8*h, gap = gap, ...
            h_solid = h);
        r = "sqrt((x - 0.5)^2 + (y - 0.5)^2)";
        normal_p = ["(x - 0.5)/" + r; "(y - 0.5)/" + r];
end

%% exact solution
mms = StokesBiotCoupling.get_manufactured_solution(normal_p, ...
    stream_function = options.stream_function, ...
    pressure_fluid = options.pressure_fluid, ...
    displacement = options.displacement, ...
    pressure_biot = options.pressure_biot, ...
    viscosity = options.viscosity, ...
    young_modulus = options.young_modulus, ...
    poisson_ratio = options.poisson_ratio, ...
    biot_alpha = options.biot_alpha, ...
    storage = options.storage, ...
    diffusion = options.diffusion, ...
    hooke_modulus = options.hooke_modulus, ...
    friction_bjs = options.friction_bjs, ...
    time_step = dt, ...
    rate = options.rate);

%% problems (data at t = 0; refreshed at t_{n+1} in the loop)
data = mms.at(0);
fluid_problem = IncompressibleElasticityProblem(k, ...
    viscosity = options.viscosity, density = 0.0, time_step = dt, ...
    source_term = data.f_F, datum_dirichlet = data.u_F, ...
    pressure_mean = "pinned");
biot_problem = BiotProblem(k, ...
    young_modulus = options.young_modulus, ...
    poisson_ratio = options.poisson_ratio, ...
    biot_alpha = options.biot_alpha, ...
    storage = options.storage, ...
    diffusion = options.diffusion, ...
    hooke_modulus = options.hooke_modulus, ...
    stabilization_darcy = options.stabilization_darcy, ...
    stabilization_solid = options.stabilization_solid, ...
    time = 0.0, time_step = dt, time_scheme = options.time_scheme, ...
    source_term = data.f_P, source_mass = data.g_P, ...
    initial_displacement = data.u_P, initial_pressure = data.p_P, ...
    datum_dirichlet = {data.u_P{1}; data.u_P{2}; data.p_P});

%% time loop
if options.rate == "exact"
    state = initial_state(biot_problem, biot_mesh, k);
    t = state.time;
else
    state = struct.empty;     % the zero history at Data.time
    t = 0.0;
end
for n = 1:options.n_steps
    t = t + dt;
    data = mms.at(t);
    fluid_problem.Data.source_term = data.f_F;
    fluid_problem.Data.datum_dirichlet = data.u_F;
    biot_problem.Data.source_term = data.f_P;
    biot_problem.Data.source_mass = data.g_P;
    biot_problem.Data.datum_dirichlet = {data.u_P{1}; data.u_P{2}; data.p_P};
    pore_pressure_mean = options.pore_pressure_mean;
    if isequal(pore_pressure_mean, "exact")
        pore_pressure_mean = mesh_mean(data.p_P, biot_mesh);
    end
    [fluid_solution, biot_solution, state, global_system] = StokesBiotCoupling.solve( ...
        fluid_problem, biot_problem, fluid_mesh, biot_mesh, curve, data, ...
        Order = k, friction_bjs = options.friction_bjs, state = state, ...
        pore_pressure_mean = pore_pressure_mean);
end

%% errors at the last step
ff = fluid_solution.Unknowns.field;
bf = biot_solution.Unknowns.field;
err = @(field, mesh, exact) compute_l2_error(field, mesh, Datum(exact));
out.errors = struct( ...
    "sigma_F", err(ff(1), fluid_mesh, data.sigma_F), ...
    "u_F",     err(ff(2), fluid_mesh, data.u_F), ...
    "sigma_P", err(bf(1), biot_mesh, data.sigma_P), ...
    "u_P",     err(bf(2), biot_mesh, data.u_P), ...
    "q_P",     err(bf(4), biot_mesh, data.q_P), ...
    "p_P",     err(bf(5), biot_mesh, data.p_P));
out.h = h;
free = global_system.dof;
out.residue = global_system.residue / max(norm(global_system.rhs(free)), realmin);
out.solutions = struct("fluid", fluid_solution, "biot", biot_solution, ...
    "state", state, "global_system", global_system, ...
    "meshes", {{fluid_mesh, biot_mesh}}, "data", data);
end

function m = mesh_mean(f, mesh)
% the mean of the handle f over the mesh (its P0 projection)
proj = project(Datum(f), mesh, DubinerBasis(0));
m = sum(proj.dof(1, :) .* mesh.volume) * dubiner_basis(0, 0, 0) / sum(mesh.volume);
end
