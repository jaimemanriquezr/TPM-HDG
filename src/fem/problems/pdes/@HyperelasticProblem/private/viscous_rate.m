function viscous = viscous_rate(pde, mesh, k, R_u, R_hat)
% The data of the Kelvin–Voigt term at the step being solved: with
% μ_v = Data.solid_viscosity > 0 and Data.density > 0 (Newmark), the
% discrete velocity gradient of the new step is affine in its discrete
% gradient l (PEERS dofs, d2peers x 1 x N),
%   l̇ = rate l + l_history,   rate = γ/(βΔt),
%   l_history = R_u w + R_û ŵ,
%   w = −rate u⁰ + (1 − γ/β) v⁰ + Δt (1 − γ/(2β)) a⁰,
% (w, ŵ) from the volume and skeleton histories of Data.state.solid (zero
% when absent). Returns struct(viscosity, rate, l_history, model), model
% = Data.viscous_model ("material" when absent; viscous_stress), or
% struct.empty without viscosity.
viscous = struct.empty;
data = pde.Data;
model = "material";
if isfield(data, "viscous_model")
    model = string(data.viscous_model);
end
if ~(isscalar(model) && ismember(model, ["material", "spatial"]))
    error("HyperelasticProblem:viscousModel", ...
        'Data.viscous_model must be "material" or "spatial".');
end
if ~isfield(data, "solid_viscosity") || data.solid_viscosity == 0 || data.density == 0
    return
end
beta = data.newmark_beta;
gamma = data.newmark_gamma;
dt = data.time_step;
rate = gamma/(beta*dt);
N_ELEMENTS = size(mesh.elements, 2);
N_FACES = size(mesh.faces, 2);
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
s = struct();
if isfield(data, "state") && ~isempty(data.state) && isfield(data.state, "solid")
    s = data.state.solid;
end
history = @(name, template) history_of(s, name, template);
volume = zeros(2*d2, N_ELEMENTS);
skeleton = zeros(2*d1, N_FACES);
w = -rate*history("displacement", volume) + (1 - gamma/beta)*history("velocity", volume) ...
    + dt*(1 - gamma/(2*beta))*history("acceleration", volume);
w_hat = -rate*history("displacement_skeleton", skeleton) ...
    + (1 - gamma/beta)*history("velocity_skeleton", skeleton) ...
    + dt*(1 - gamma/(2*beta))*history("acceleration_skeleton", skeleton);
l_history = pagemtimes(R_u, reshape(w, 2*d2, 1, N_ELEMENTS)) ...
    + pagemtimes(R_hat, element_traces(mesh, w_hat, d1));
viscous = struct("viscosity", data.solid_viscosity, "rate", rate, "l_history", l_history, ...
    "model", model);
end

function A = history_of(s, name, template)
% a state.solid array, zero (shaped like the template) when absent
A = template;
if isfield(s, name) && ~isempty(s.(name)) && ~isequal(s.(name), 0) && ~iscell(s.(name))
    A = s.(name);
end
end
