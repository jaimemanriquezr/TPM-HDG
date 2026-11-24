function spaces = get_discrete_spaces(k)
%-------------------------------------------------------------------------
%Description: Returns FE space descriptors for the Darcy problem.
%             Volume:  vector Dubiner P_k^2 for the flux u;
%                      scalar Dubiner P_k  for pressure p.
%             Hybrid:  scalar Legendre P_k for pressure trace p̂.
%-------------------------------------------------------------------------
arguments (Input)
    k   double  % polynomial degree
end
arguments (Output)
    spaces  struct  % struct with 'volume' (volume spaces) and 'hybrid' (skeleton spaces)
end
if isempty(k); spaces = struct.empty; return; end
spaces = struct( ...
    "volume", [DubinerBasis(k).to_vector(2), DubinerBasis(k)], ...
    "hybrid", LegendreBasis(k));
end
