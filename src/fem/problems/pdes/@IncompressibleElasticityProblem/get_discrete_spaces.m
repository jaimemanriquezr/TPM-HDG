function spaces = get_discrete_spaces(k)
%-------------------------------------------------------------------------
%Description: Returns FE space descriptors for the incompressible elasticity
%             (stress-velocity) problem. Uses PEERS for stress.
%             Volume:  PEERS P_k             for stress σ;
%                      vector Dubiner P_k^2  for velocity u;
%                      scalar Dubiner P_k    for rotation γ.
%             Hybrid:  vector Legendre P_k^2 for velocity trace û;
%                      scalar P_0           for per-element pressure mean ρ.
%-------------------------------------------------------------------------
arguments (Input)
    k   double  % polynomial degree
end
arguments (Output)
    spaces  struct  % struct with 'volume' (volume spaces) and 'hybrid' (skeleton spaces)
end
if isempty(k); spaces = struct.empty; return; end
spaces = struct( ...
    "volume", [PEERSBasis(k), DubinerBasis(k).to_vector(2), DubinerBasis(k)], ...
    "hybrid", [LegendreBasis(k).to_vector(2), DubinerBasis(0)]);
end
