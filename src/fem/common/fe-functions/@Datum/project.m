function [l2_proj, l2_proj_error, A_obj, b_obj] = project(obj, mesh, basis, options)
%PROJECT  L2 projection of a Datum onto a finite element basis.
%
%   l2_proj = project(obj, mesh, basis)
%   [l2_proj, l2_proj_error] = project(obj, mesh, basis, order=, faces=)
%   [~, ~, A_obj, b_obj] = project(obj, mesh, basis, ...)
%
%   Solves A_obj \ b_obj, the mass system of basis against obj, over the
%   elements of mesh when basis is two-dimensional or over its faces when
%   it is one-dimensional. l2_proj is an FEField on FESpace(basis, mesh).
%   basis defaults to LegendreBasis(0).
%
%   order is the quadrature order (default 10); faces restricts a
%   one-dimensional basis to a subset of faces (empty = all).
%   l2_proj_error, the L2 error of the projection, is computed only when
%   exactly two outputs are requested and is [] otherwise. A_obj and
%   b_obj are the bilinear and linear forms of the system.
arguments
    obj Datum
    mesh AbstractGrid
    basis AbstractFEBasis = LegendreBasis(0)
    options.order = 10;
    options.faces = [];     % facet subset for 1D bases (empty = all faces)
end

    u = TrialFunction(basis);
    v = TestFunction(basis);
    switch basis.GeometricInputDimension
        case 1
            dmu = Measure(mesh, 'ds');
        case 2
            dmu = Measure(mesh, 'dx');
    end
    dmu.order = options.order;
    dmu.faces = options.faces;
    if prod(basis.GeometricOutputDimension) > 1
        % vector/tensor values: contract the components
        b_obj = dot(obj, v) * dmu;
        A_obj = dot(u, v) * dmu;
    else
        b_obj = obj * v * dmu;
        A_obj = u * v * dmu;
    end
    l2_proj_dof = A_obj \ b_obj;
    if basis.GeometricInputDimension == 1
        % 'ds' gathers every face from both neighbours (boundary faces
        % repeated); the two slices coincide, keep one per face
        l2_proj_dof = l2_proj_dof(:, :, :, :, :, 1);
    end
space = FESpace(basis, mesh);
l2_proj = FEField(space, l2_proj_dof);

if nargout == 2
    l2_proj_error = l2_proj.compute_l2_error(mesh, obj);
else
    % callers requesting A_obj/b_obj (e.g. weak boundary conditions) do not
    % need the projection error; skip it (not defined for facet fields)
    l2_proj_error = [];
end
end
