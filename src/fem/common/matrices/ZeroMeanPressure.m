function [H, RHS] = ZeroMeanPressure(obj, mesh, H, RHS)
% Enforces ∫_Ω p dΩ = 0 by appending a Lagrange multiplier
% row/column; no-op when data.pressure_mean = "pinned".
if isfield(obj.Data, "pressure_mean") && obj.Data.pressure_mean == "pinned"
    return
end
N_ELEMENTS = size(mesh.elements, 2);
H(end-N_ELEMENTS+1:end, end+1) = mesh.volume;
H(end+1, end-N_ELEMENTS:end-1) = mesh.volume;
RHS(end+1) = 0;
end
