function u_proj = project(obj, mesh, basis)
%PROJECT  L2 projection of an FE expression onto a basis on a mesh.
%
%   u_proj = project(obj, mesh, basis)
%
%   Solves the mass system of basis against obj (quadrature order 10, on
%   elements for a basis of geometric input dimension 2, on faces for 1)
%   and returns the FEField on FESpace(basis, mesh).
    arguments
        obj AbstractFEFunction
        mesh AbstractGrid
        basis AbstractFEBasis
    end
    switch basis.GeometricInputDimension
        case 2
            dmu = measure(mesh, 'element:interior');
        case 1
            dmu = measure(mesh, 'element:faces');
    end
    % the quadrature of Datum.project
    dmu.order = 10;
    u = get_trial_functions(basis);
    v = get_test_functions(basis);

    mass = dot(u, v) * dmu;
    if isa(basis, "TensorFEBasis")
        % vector/tensor values: contract the components
        rhs = dot(obj, v) * dmu;
    else
        rhs = obj * v * dmu;
    end
    dof = mass \ rhs;
    if basis.GeometricInputDimension == 1
        % 'ds' sees every face from both neighbours; keep one slice
        dof = dof(:, :, :, :, :, 1);
    end
    u_proj = FEField(FESpace(basis, mesh), dof);
end
