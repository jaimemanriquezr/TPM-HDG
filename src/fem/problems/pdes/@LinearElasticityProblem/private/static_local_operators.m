function [A1, A2, Af] = static_local_operators(obj, mesh, k)
% the local operators without the inertial terms
obj.Data.density = 0;
local_operators = get_matrices(obj, mesh, k);
[A1, A2, Af] = deal(local_operators{:});
end
