function gN = neumann_boundary_vector(mesh, varargin)
% NEUMANN_BOUNDARY_VECTOR Neumann datum tested on the boundary skeleton.
%   gN = neumann_boundary_vector(mesh, datum, basis) - Datum/AbstractFEBasis
%       form (DSL path).
%   gN = neumann_boundary_vector(mesh, q_fh, k, basis=@legendre_basis)
%       function-handle form (low-level path, k = degree).
if nargin >= 3 && isa(varargin{2}, 'double')
    gN = handle_neumann_boundary_vector(mesh, varargin{:});
else
    gN = datum_neumann_boundary_vector(mesh, varargin{:});
end
end

function gN = datum_neumann_boundary_vector(mesh, neumann_datum, basis)
arguments
    mesh SimplicialGrid
    neumann_datum Datum = @(x,y) 0 * (x .* y)
    basis AbstractFEBasis = LegendreBasis(0)
end
if isempty(mesh.neufaces)
    gN = [];
else
    gN = neumann_datum.test(mesh, basis);
    gN = gN(:, :, mesh.neufaces);
end
end

function gN = handle_neumann_boundary_vector(mesh, datum_neumann, k, options)
arguments
    mesh SimplicialGrid
    datum_neumann function_handle = @(x,y) 0 * (x .* y)
    k double = 0
    options.basis = @legendre_basis
end
if isempty(mesh.neufaces)
    gN = [];
else
    gN = test_faces(mesh, datum_neumann, k, subset=mesh.neufaces, basis=options.basis);
end
end
