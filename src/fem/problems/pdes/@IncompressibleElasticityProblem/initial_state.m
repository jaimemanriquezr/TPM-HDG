function state = initial_state(obj, mesh, k, options)
% Time-stepping state at t = obj.Data.time from the initial data.
%
%   state = initial_state(obj, mesh, k)                  (k defaults to obj.Order)
%   state = initial_state(obj, mesh, k, with_pressure = true)
%
% Projects Data.initial_velocity (a cell of handles {u_x; u_y}):
% u⁰ = Π u(0) onto the elements (L2proj2D) and û⁰ = P_M u(0) onto the
% faces (L2projSkeleton1d). When θ < 1 (the step reads a⁰) or
% with_pressure = true it also solves the acceleration problem
% (AccelerationMatrices, density > 0 required)
%   ρ a⁰ + ∇p⁰ = f(0) + ∇·σ_v⁰ [- ρ (β·∇)u⁰],   ∇·a⁰ = 0,
% with σ_v⁰ = 2μ ε_h(u⁰) the local lifting of (u⁰, û⁰) (lift_stress) and,
% when Data.convective_velocity is set, the convective term with β = u⁰
% (rho_c = density) in the momentum rows and, as the step's numerical
% flux, -ρ(û⁰ ⊗ β)n in the flux rows (not on traction faces). It reads
% Data.source_term = f(0), the traction datum Data.datum_neumann = q_N(0)
% and, on the Dirichlet faces, a⁰ = Data.datum_dirichlet_rate = ∂_t g(0):
% set them to their t = 0 values before the call. With
% Data.transfer_mapping the Dirichlet datum is transferred along the
% paths with the gradient of the local acceleration polynomial
% (TransferEssentialBoundaryConditions.field_gradient); that gradient is
% only O(h^k), and a⁰ then converges at about rate k instead of k + 1
% (observed for k >= 2; p⁰ and σ⁰ lose less). Componentwise
% Data.boundary_conditions keep their types; every Dirichlet component
% takes datum_dirichlet_rate (one handle pair for all tags) and every
% traction component its own data; they are not supported together with
% a transfer mapping. The interface_conditions option (an array of
% boundary-condition objects, default empty) is appended to these: it
% carries conditions on faces the default tags do not cover, e.g. the
% interface faces of a coupled problem (see
% FluidStructureCoupling.initial_state).
% Returns state.time = Data.time, state.step = 0, state.solid =
% struct.empty, state.diagnostics = struct.empty and state.fluid with
%   velocity           u⁰ (2*d2 x N)
%   velocity_skeleton  û⁰ (2*d1 x F)
%   acceleration       a⁰ = u_t(0) (2*d2 x N), or []
%   pressure           p⁰ (d2 x N), or []
%   stress             σ⁰ = σ_v⁰ - p⁰ I (PEERS dofs, d2peers x N), or []
% ([] when the acceleration problem is not solved).
arguments
    obj  (1,1) IncompressibleElasticityProblem
    mesh       SimplicialGrid
    k    (1,1) double = obj.Order
    options.with_pressure (1,1) logical = false
    options.interface_conditions = []
end
data = obj.Data;
assert(iscell(data.initial_velocity), "IncompressibleElasticityProblem:initialVelocity", ...
    "initial_state needs Data.initial_velocity as a cell of handles {u_x; u_y}.");
u0  = L2proj2D(mesh, data.initial_velocity, k);
u0h = L2projSkeleton1d(mesh, data.initial_velocity, k);
fluid = struct( ...
    "velocity",          u0,  ...
    "velocity_skeleton", u0h, ...
    "acceleration",      [],  ...
    "pressure",          [],  ...
    "stress",            []);
state = struct("time", data.time, "step", 0, "solid", struct.empty, ...
    "fluid", fluid, "diagnostics", struct.empty);
if obj.get_theta() == 1 && ~options.with_pressure
    return
end

assert(data.density > 0, "IncompressibleElasticityProblem:density", ...
    "The acceleration problem of initial_state requires density > 0.");
[a0, p0, sigma0] = acceleration_problem(data, mesh, k, u0, u0h, options.interface_conditions);
state.fluid.acceleration = a0;
state.fluid.pressure     = p0;
state.fluid.stress       = sigma0;
end

function [a0, p0, sigma0] = acceleration_problem(data, mesh, k, u0, u0h, interface_conditions)
N_ELEMENTS = size(mesh.elements, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
rho = data.density;
has_specs = isfield(data, "boundary_conditions") && ~isempty(data.boundary_conditions);
has_transfer = isfield(data, "transfer_mapping") && ~isempty(data.transfer_mapping);

% the momentum source moments (f(0) + ∇·σ_v⁰, v)_K
[sigma_v, ~, div_sigma] = lift_stress(mesh, k, data.viscosity, u0, u0h);
F = reshape([test_elements(mesh, data.source_term{1}, k);
             test_elements(mesh, data.source_term{2}, k)], 2*d2, N_ELEMENTS) + div_sigma;

% the convective term as a known source, β = u⁰: +(Cb u⁰ - BbLP û⁰) in
% the momentum rows and +BbLL û⁰ in the flux rows, rho_c = density, the
% flux skipped on the natural faces as the step does
flux_moments = [];
if isfield(data, "convective_velocity") && ~isempty(data.convective_velocity)
    natural_faces = mesh.neufaces(:);
    if has_specs
        for s = data.boundary_conditions(:).'
            if any(s.type == "neumann") && isKey(mesh.boundary_tags, s.tag)
                natural_faces = [natural_faces; reshape(get_face_by_tag(mesh, s.tag), [], 1)]; %#ok<AGROW>
            end
        end
    end
    beta = {u0(1:d2, :); u0(d2+1:end, :)};
    [Cb, BbPL, BbLL] = convective_matrix(mesh, k, beta, skip_faces = natural_faces);
    BbLP = rho*permute(BbPL, [2 1 3]);
    Cb = rho*Cb;
    BbLL = rho*BbLL;
    u0h_K = cat(1, reshape(u0h(1:d1, mesh.boundaries), 3*d1, 1, N_ELEMENTS), ...
                   reshape(u0h(d1+1:end, mesh.boundaries), 3*d1, 1, N_ELEMENTS));
    conv = [pagemtimes(Cb, reshape(u0(1:d2, :), d2, 1, N_ELEMENTS)) ...
              - pagemtimes(BbLP, u0h_K(1:3*d1, :, :));
            pagemtimes(Cb, reshape(u0(d2+1:end, :), d2, 1, N_ELEMENTS)) ...
              - pagemtimes(BbLP, u0h_K(3*d1+1:end, :, :))];
    F = F + reshape(conv, 2*d2, N_ELEMENTS);
    flux_moments = reshape([pagemtimes(BbLL, u0h_K(1:3*d1, :, :));
                            pagemtimes(BbLL, u0h_K(3*d1+1:end, :, :))], 6*d1, N_ELEMENTS);
end

closure = data.pressure_mean;
[lhs, flux, rhs] = AccelerationMatrices(k, mesh, rho, F, ...
    pressure_mean = closure, viscous_stress = sigma_v, flux_moments = flux_moments);
Pk = DubinerBasis(k);
global_spaces = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
if closure == "constrained"
    global_spaces = [global_spaces, ConstantBasis()];
end
problem = FEProblem(Order = k, Spaces = [ ...
    FESpace([Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
    FESpace(global_spaces)]);
problem.Mesh = mesh;
problem.Operator = AffineOperator(lhs + rhs) + flux;
problem.BoundaryConditions = [boundary_conditions(data, has_specs, has_transfer), ...
    interface_conditions];
solution = solve(problem, mesh, k);

a0 = solution.Unknowns.field(1).dof;
p0 = solution.Unknowns.field(2).dof;
% σ⁰ = σ_v⁰ - p⁰ I: the identity part lives in the σ_11, σ_22 Dubiner blocks
sigma0 = sigma_v;
sigma0(1:d2, :)        = sigma0(1:d2, :) - p0;
sigma0(3*d2 + (1:d2), :) = sigma0(3*d2 + (1:d2), :) - p0;
end

function bc = boundary_conditions(data, has_specs, has_transfer)
% the step's boundary conditions with the rate datum on the Dirichlet
% components and the traction data on the natural ones
rate = data.datum_dirichlet_rate;
if has_specs
    if has_transfer
        error("IncompressibleElasticityProblem:transferComponentwise", ...
            "initial_state does not support componentwise boundary_conditions " + ...
            "together with a transfer_mapping.");
    end
    specs = data.boundary_conditions;
    for i = 1:numel(specs)
        for c = 1:2
            if specs(i).type(c) == "dirichlet"
                specs(i).data{c} = rate{c};
            end
        end
    end
    bc = ComponentwiseBoundaryConditions(specs);
    return
end
gN = data.datum_neumann;
traction = boundary_condition("Neumann", x = "neumann", y = "neumann", ...
    x_data = gN{1}, y_data = gN{2});
if has_transfer
    bc = [ComponentwiseBoundaryConditions(traction), ...
          TransferEssentialBoundaryConditions.field_gradient(Datum(rate), "Dirichlet", ...
              data.transfer_mapping)];
else
    bc = ComponentwiseBoundaryConditions([ ...
        boundary_condition("Dirichlet", x_data = rate{1}, y_data = rate{2});
        traction]);
end
end
