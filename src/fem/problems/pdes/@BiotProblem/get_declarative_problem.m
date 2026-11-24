function fe = get_declarative_problem(obj, mesh, k)
%GET_DECLARATIVE_PROBLEM  The Biot step as a declarative FEProblem.
%
%   fe = get_declarative_problem(obj, mesh, k)
%
%   fe has Order = k (default obj.Order), the spaces of
%   BiotProblem.get_discrete_spaces (ids 1-5 condensed, id 6 the
%   three-component trace), Mesh = mesh, Data = obj.Data and
%       Operator = AffineOperator(lhs) + flux (+ rhs) + sources,
%   [lhs, flux, rhs] = BiotMatrices(...) with (cdt, history) from
%   biot_step_history(time_scheme, time_step, state.biot) and the sources
%   (f, v) ([-1, 2]) and (g, r) ([-1, 5]). An empty state.biot is the zero
%   history with cdt = time_step.
%
%See also: BiotProblem, BiotMatrices, biot_step_history
arguments
    obj  (1,1) BiotProblem
    mesh       SimplicialGrid
    k          double {mustBeScalarOrEmpty} = obj.Order
end
data = obj.Data;
P_k = DubinerBasis(k);
spaces = BiotProblem.get_discrete_spaces(k);

biot_state = struct.empty;
if ~isempty(data.state) && isfield(data.state, "biot")
    biot_state = data.state.biot;
end
if isempty(biot_state)
    % the zero history, as biot_step_history returns it for an empty state
    cdt = data.time_step;
    history = struct.empty;
else
    [cdt, history] = biot_step_history(data.time_scheme, data.time_step, biot_state);
end

[lhs, flux, rhs] = BiotMatrices(k, mesh, data.young_modulus, data.poisson_ratio, ...
    data.diffusion, biot_alpha = data.biot_alpha, storage = data.storage, ...
    effective_time_step = cdt, history = history, ...
    hooke_modulus = data.hooke_modulus, density = data.density, ...
    stabilization_darcy = data.stabilization_darcy, ...
    stabilization_solid = data.stabilization_solid);

fe = FEProblem(obj.Label, Order = k, Spaces = [ ...
    FESpace(spaces.volume, Condense = true), FESpace(spaces.hybrid)]);
fe.Mesh = mesh;
fe.Data = data;
fe.Operator = AffineOperator(lhs) + flux;
if ~isempty(rhs), fe.Operator = fe.Operator + rhs; end

w = TestFunction(P_k);
dx = Measure(mesh, 'dx');
f = data.source_term;
fe.Operator = fe.Operator + IndexedMatrix(IndexedArray(cat(1, ...
    double(Datum(f{1}) * w * dx), double(Datum(f{2}) * w * dx)), [-1, 2]));
fe.Operator = fe.Operator + IndexedMatrix(IndexedArray( ...
    double(Datum(data.source_mass) * w * dx), [-1, 5]));

if isfield(data, "boundary_conditions") && ~isempty(data.boundary_conditions)
    % componentwise specs act on eta_hat, and on p_hat when a spec sets a
    % third component (boundary_condition's p); otherwise p_hat is
    % natural (zero flux)
    fe.BoundaryConditions = ComponentwiseBoundaryConditions(data.boundary_conditions);
elseif ~isempty(data.transfer_mapping)
    fe.BoundaryConditions = TransferEssentialBoundaryConditions( ...
        Datum(data.datum_dirichlet), "Dirichlet", data.transfer_mapping, ...
        path_gradient = biot_path_gradient(k, data.young_modulus, ...
            data.poisson_ratio, data.biot_alpha, data.diffusion));
else
    fe.BoundaryConditions = StrongEssentialBoundaryConditions( ...
        Datum(data.datum_dirichlet), "Dirichlet", projection_order = k);
end
end
