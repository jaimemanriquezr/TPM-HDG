% test_input_handling.m
% Documented inputs either work or raise a named error, and no option is
% silently dropped. One check per input path:
%   - test_faces default subset is all faces;
%   - L2projSkeleton1d: the subset on the handle and cell paths;
%   - Datum.test with element_indices;
%   - IndexedArray - IndexedMatrix, and IndexedMatrix indexed assignment;
%   - Measure default type; WeightFunction base evaluate;
%   - StrongEssentialBoundaryConditions with two tags on a vector trace
%     equals the one-tag solve with the same datum;
%   - SimplicialGrid concatenation; RectangleMesh side tags;
%   - FluidStructureCoupling.solve with different fluid and solid time
%     steps;
%   - quadpoints_elements_2d / test_elements on an element subset;
%   - the 4-row (per-edge size) SimplicialGrid connectivity;
%   - submesh keeps a parent-boundary tag and retags only exposed faces;
%   - CircularCurve and the manufactured fluid measure at full precision;
%   - get_manufactured_solution on an unsupported interface;
%   - interface_meshes with a gap on a clockwise polygon;
%   - the one-step Taylor start, selected by the state, not the time.
clc; clear;

n_fail = 0;
mesh = UnitSquareMesh(.25);
n_faces = size(mesh.faces, 2);
g = @(x, y) sin(x).*y;

%% test_faces: the default subset is every face
k = 1;
n_fail = n_fail + check(isequal(test_faces(mesh, g, k), test_faces(mesh, g, k, subset = 1:n_faces)), ...
    "test_faces default subset");

%% L2projSkeleton1d: subsets
P = L2projSkeleton1d(mesh, g, k);
f = [2 5 7];
P_sub = L2projSkeleton1d(mesh, g, k, subset = f);
n_fail = n_fail + check(max(abs(P_sub - P(:, f)), [], "all") < 1e-14, "L2projSkeleton1d subset (handle)");
P_cell = L2projSkeleton1d(mesh, {g, g}, k, subset = f);
n_fail = n_fail + check(isequal(size(P_cell), [2*(k + 1), numel(f)]), "L2projSkeleton1d subset (cell)");

%% Datum.test with element_indices
I = Datum({@(x, y) x + 0*y}).test(mesh, DubinerBasis(1), element_indices = 1:3);
n_fail = n_fail + check(size(I, 3) == 3, "Datum.test element_indices");

%% IndexedArray - IndexedMatrix; IndexedMatrix indexed assignment
a = IndexedArray(magic(3), [1 2]);
B = IndexedMatrix(IndexedArray(eye(3), [1 2]));
n_fail = n_fail + check(isequal(double(a - B), magic(3) - eye(3)), "IndexedArray - IndexedMatrix");
n_fail = n_fail + check(throws(@() assign_entry(B), "IndexedMatrix:parenAssign"), ...
    "IndexedMatrix indexed assignment");

%% Measure default type; WeightFunction base evaluate
n_fail = n_fail + check(Measure(mesh).type == "dx", "Measure default type");
n_fail = n_fail + check(throws(@() WeightFunction().evaluate(zeros(2, 1), mesh), ...
    "WeightFunction:abstract"), "WeightFunction base evaluate");

%% StrongEssentialBoundaryConditions: two tags on a vector trace
u_bc = Datum({@(x, y) x.*y; @(x, y) x - y});
nu = 0.7;
src = {@(x, y) sin(pi*x).*cos(pi*y); @(x, y) exp(x).*y - exp(1)/2 + 1/2};
grid = StructuredSimplicialGrid(8, 8).add_stabilization(1);
grid_2 = grid.apply_tag("Inflow", @(x) abs(x(2, :) - 1) < 1e-12);
u_one = strong_stokes_velocity(grid, StrongEssentialBoundaryConditions(u_bc, "Dirichlet"), k, nu, src);
u_two = strong_stokes_velocity(grid_2, ...
    StrongEssentialBoundaryConditions([u_bc, u_bc], ["Dirichlet", "Inflow"]), k, nu, src);
d = max(abs(u_two - u_one), [], "all");
n_fail = n_fail + check(d < 1e-10, sprintf("strong BC, two tags vs one tag |du| = %.2e", d));

%% SimplicialGrid concatenation; RectangleMesh side tags
n_fail = n_fail + check(throws(@() [mesh, mesh], "SimplicialGrid:cat"), "SimplicialGrid concatenation");
n_fail = n_fail + check(~isempty(UnitSquareMesh(.25, top_side = "Inflow").get_face_by_tag("Inflow")), ...
    "RectangleMesh top_side = ""Inflow""");
n_fail = n_fail + check(throws(@() UnitSquareMesh(.25, top_side = "Bogus"), "RectangleMesh:unknownTag"), ...
    "RectangleMesh unknown side tag");

%% FluidStructureCoupling: fluid and solid time steps must agree
curve = FlatCurveHorizontal(0.5);
[fluid_mesh, solid_mesh] = fitted_split(curve, [.25 .25], 0);
[fp, sp, sources] = FluidStructureCoupling.load_problems(1, curve, 0.1, 1.0, ...
    "x", "0", ["1"; "0"], 1, 1.0, 1, 1.0, 0.25, 0.0);
sp.Data.time_step = 2*fp.Data.time_step;
n_fail = n_fail + check(throws(@() FluidStructureCoupling.solve(fp, sp, fluid_mesh, solid_mesh, ...
    curve, sources{:}, "steady"), "FluidStructureCoupling:timeStep"), "FSI time-step mismatch");

%% element subsets
[xy, ~, ~] = quadpoints_elements_2d(mesh, k);
[xy_sub, ~, ~] = quadpoints_elements_2d(mesh, k, [2 5]);
n_fail = n_fail + check(isequal(xy_sub{1}, xy{1}(:, [2 5])) && isequal(xy_sub{2}, xy{2}(:, [2 5])), ...
    "quadpoints_elements_2d subset");
I_all = test_elements(mesh, g, k);
I_sub = test_elements(mesh, g, k, subset = [2 5]);
n_fail = n_fail + check(max(abs(I_sub - I_all(:, :, [2 5])), [], "all") < 1e-15, "test_elements subset");

%% 4-row connectivity: per-edge sizes (0 = element_diameter)
square = SimplicialGrid([0 1 1 0; 0 0 1 1], [1:4; 2:4, 1; 1 1 1 1; 0.125 0 0 0], element_diameter = 0.5);
bottom = square.get_face_by_tag("Dirichlet");
on_bottom = all(abs(reshape(square.coordinates(2, square.faces(1:2, bottom)), 2, [])) < 1e-12, 1);
n_fail = n_fail + check(nnz(on_bottom) == 8, sprintf("4-row connectivity: %d bottom faces", nnz(on_bottom)));

%% submesh: a parent-boundary tag survives, the cut is tagged Interface
tagged = UnitSquareMesh(.25).apply_tag("Inflow", @(x) abs(x(1, :)) < 1e-12);
centroid_x = mean(reshape(tagged.coordinates(1, tagged.elements), 3, []));
left = tagged.submesh(find(centroid_x < 0.5));
iface = left.get_face_by_tag("Interface");
n_fail = n_fail + check(numel(left.get_face_by_tag("Inflow")) == numel(tagged.get_face_by_tag("Inflow")), ...
    "submesh keeps the Inflow tag");
n_fail = n_fail + check(~isempty(iface) && all(left.coordinates(1, left.faces(1:2, iface)) > 0.25 - 1e-12, "all"), ...
    "submesh tags only the exposed faces");

%% full-precision circle and fluid measure; unsupported interfaces
circle = CircularCurve([0.5; 0.5], 1/3);
n_fail = n_fail + check(max(abs(circle.get_2d_points(0) - [0.5 + 1/3; 0.5])) < 1e-15, "CircularCurve precision");
[~, ~, ~, pd] = FluidStructureCoupling.load_problems(1, circle, 0.1, 1.0, ...
    "x", "0", ["1"; "0"], 1, 1.0, 1, 1.0, 0.25, 0.0);
n_fail = n_fail + check(abs(pd.measure_F - (1 - pi/9)) < 1e-15, ...
    sprintf("manufactured measure_F error %.2e", abs(pd.measure_F - (1 - pi/9))));
n_fail = n_fail + check(throws(@() FluidStructureCoupling.load_problems(1, AirfoilCurve(), 0.1, 1.0, ...
    "x", "0", ["1"; "0"], 1, 1.0, 1, 1.0, 0.25, 0.0), "get_manufactured_solution:unsupportedInterface"), ...
    "manufactured solution on an airfoil");

%% gap on a clockwise polygon: the same meshes as counter-clockwise
P = [0 1 1 0; 0 0 1 1]*0.4 + 0.3;
[outer_cw, inner_cw] = interface_meshes(PolygonalCurve(fliplr(P)), "fitted", 0.1, gap = 0.01);
[outer_ccw, inner_ccw] = interface_meshes(PolygonalCurve(P), "fitted", 0.1, gap = 0.01);
n_fail = n_fail + check(abs(sum(outer_cw.volume) - sum(outer_ccw.volume)) < 1e-12 ...
    && abs(sum(inner_cw.volume) - sum(inner_ccw.volume)) < 1e-12, "gap on a clockwise polygon");

%% one-step Taylor start: selected by the state, independent of the time
solid = LinearElasticityProblem(1, density = 1, time_step = 0.1, ...
    initial_displacement = {@(x, y) x.*y; @(x, y) 0*x}, initial_velocity = {@(x, y) 0*x; @(x, y) y});
solid.Data.datum_dirichlet = solid.Data.initial_displacement;
state = initial_state(solid, mesh, 1);
solid.Data.state = state;
u_start = solve(solid, mesh, 1).Unknowns.field(2).dof;
solid.Data.time = 0.7;
u_late = solve(solid, mesh, 1).Unknowns.field(2).dof;
n_fail = n_fail + check(isequal(u_start, u_late), "one-step start independent of Data.time");
solid.Data.state.solid.displacement_previous = state.solid.displacement;
u_march = solve(solid, mesh, 1).Unknowns.field(2).dof;
n_fail = n_fail + check(max(abs(u_march - u_start), [], "all") > 1e-6, "one-step marching after the start");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_input_handling:failed", "%d checks failed", n_fail);
end

function u = strong_stokes_velocity(mesh, bc, k, nu, src)
% the declarative Stokes system of test_dsl_stokes with the given
% boundary conditions; returns the velocity dofs
Pk = DubinerBasis(k);
problem = FEProblem(Spaces = [ ...
    FESpace([Pk.to_tensor(), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
    FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])]);
problem.Mesh = mesh;
[lhs, flux, rhs] = StokesMatrices(k, mesh, nu, src);
problem.Operator = AffineOperator(lhs + rhs) + flux;
problem.BoundaryConditions = bc;
sol = solve(problem, mesh, k);
u = sol.Unknowns.field(2).dof;
end

function B = assign_entry(B)
B(1, 1) = 3;
end
