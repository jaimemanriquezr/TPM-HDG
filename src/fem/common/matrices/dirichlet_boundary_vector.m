function gD = dirichlet_boundary_vector(mesh, varargin)
% DIRICHLET_BOUNDARY_VECTOR Dirichlet datum on the boundary skeleton.
%   gD = dirichlet_boundary_vector(mesh, datum, basis) - Datum/AbstractFEBasis
%       form (DSL path).
%   gD = dirichlet_boundary_vector(mesh, g_fh, k, basis=@legendre_basis)
%       function-handle form (low-level path, k = degree).
%   The two forms differ in quadrature: the Datum form projects with
%   Datum.project's default order (10), the handle form with order k
%   (k+1 Gauss points), so their dofs agree only up to quadrature error.
if nargin >= 3 && isa(varargin{2}, 'double')
    gD = handle_dirichlet_boundary_vector(mesh, varargin{:});
else
    gD = datum_dirichlet_boundary_vector(mesh, varargin{:});
end
end

function gD = datum_dirichlet_boundary_vector(mesh, dirichlet_datum, basis)
arguments
    mesh SimplicialGrid
    dirichlet_datum Datum = @(x, y) 0 * x .* y
    basis AbstractFEBasis = LegendreBasis(0)
end
if isempty(mesh.dirfaces)
    gD = [];
else
    proj_dirichlet_datum = dirichlet_datum.project(mesh, basis);
    gD = proj_dirichlet_datum.dof(:, mesh.dirfaces);
end
end

function gD = handle_dirichlet_boundary_vector(mesh, datum_dirichlet, k, options)
arguments
    mesh SimplicialGrid
    datum_dirichlet function_handle = @(x, y) 0 * x .* y
    k double = 0
    options.basis function_handle = @legendre_basis
end
if isempty(mesh.dirfaces)
    gD = [];
else
    gD = L2projSkeleton1d(mesh, datum_dirichlet, k, subset=mesh.dirfaces, basis=options.basis);
end
end
