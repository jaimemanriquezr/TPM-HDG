function spaces = get_discrete_spaces(k)
%-------------------------------------------------------------------------
%Description: Returns FE space descriptors for the linear elasticity problem.
%             Uses the PEERS element for stress to ensure a locking-free discretization.
%             Volume:  PEERS P_k             for stress σ;
%                      vector Dubiner P_k^2  for displacement u;
%                      scalar Dubiner P_k    for rotation γ.
%             Hybrid:  vector Legendre P_k   for displacement trace û (2 components).
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
    "hybrid", LegendreBasis(k).to_vector(2));
end
