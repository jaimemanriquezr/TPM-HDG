%-------------------------------------------------------------------------
%Description: Zero history fields of the time-dependent solid's state.
%
%             Each volume history of state.solid (displacement u0,
%             velocity v0, acceleration a0, displacement_previous u_{-1})
%             is a dof array, or [] or the scalar 0 meaning "zero", and is
%             read independently of the others. For every scheme
%             (newmark; onestep start; onestep marching) and through
%             both the LinearElasticityProblem class (Af of get_matrices,
%             Data.state) and the LinearElasticityMatrices factory (the
%             history rhs, the state option):
%               (i)  u0 = [] (and u0 = 0) with the other histories arrays
%                    must equal the run with u0 = an explicit zeros array;
%               (ii) u0 an array with the other histories [] must equal
%                    the run with explicit zeros arrays (under onestep
%                    marching u_{-1} = 0: [] there selects the start).
%-------------------------------------------------------------------------
n_failed = 0;

rng(0);
k = 1;
mesh = StructuredSimplicialGrid(2, 2);
n_el = size(mesh.elements, 2);
d2 = nchoosek(k + 2, 2);
Z = zeros(2*d2, n_el);
U = rand(2*d2, n_el);
V = rand(2*d2, n_el);
A = rand(2*d2, n_el);
P = rand(2*d2, n_el);

common = {"density", 1.0, "time_step", 0.1};
schemes = { ...
    "newmark",           {"time_scheme", "newmark"}, ...
        struct("velocity", V, "acceleration", A), ...
        struct("velocity", [], "acceleration", []), ...
        struct("velocity", Z, "acceleration", Z);
    "onestep start",     {"time_scheme", "onestep"}, ...
        struct("velocity", V), ...
        struct("velocity", []), ...
        struct("velocity", Z);
    "onestep marching",  {"time_scheme", "onestep"}, ...
        struct("displacement_previous", P), ...
        struct("displacement_previous", 0.0), ...
        struct("displacement_previous", Z)};
paths = {"class", @class_rhs; "factory", @factory_rhs};

for i = 1:size(schemes, 1)
    [name, scheme, others, others_empty, others_zero] = schemes{i, :};
    fprintf("%s\n", name);
    for j = 1:size(paths, 1)
        [path, rhs] = paths{j, :};
        opts = [common, scheme];
        explicit = rhs(mesh, k, opts, with(others, "displacement", Z));
        n_failed = n_failed + check(name + ", " + path + ": (i) u0 = [], other histories arrays", ...
            @() isequal(rhs(mesh, k, opts, with(others, "displacement", [])), explicit));
        n_failed = n_failed + check(name + ", " + path + ": (i) u0 = 0, other histories arrays", ...
            @() isequal(rhs(mesh, k, opts, with(others, "displacement", 0.0)), explicit));
        n_failed = n_failed + check(name + ", " + path + ": (ii) u0 array, other histories []", ...
            @() isequal(rhs(mesh, k, opts, with(others_empty, "displacement", U)), ...
                        rhs(mesh, k, opts, with(others_zero,  "displacement", U))));
    end
end

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function s = with(s, name, value)
s.(name) = value;
end

function b = class_rhs(mesh, k, opts, solid)
pde = LinearElasticityProblem(k, opts{:});
pde.Data.state = struct("time", pde.Data.time, "step", 0, "solid", solid, ...
    "fluid", struct.empty, "diagnostics", struct.empty);
local_operators = pde.get_matrices(mesh, k);
b = local_operators{3};
end

function b = factory_rhs(mesh, k, opts, solid)
[~, ~, rhs] = LinearElasticityMatrices(k, mesh, 1e4, 0.1, opts{:}, state = solid);
b = double(rhs.Arrays{1});
end

function failed = check(name, condition)
% condition is evaluated here so that an error counts as a failure
try
    passed = condition();
    message = "";
catch err
    passed = false;
    message = " (" + err.message + ")";
end
failed = ~passed;
if failed
    fprintf(2, "  FAIL  %s%s\n", name, message);
else
    fprintf("  ok    %s\n", name);
end
end
