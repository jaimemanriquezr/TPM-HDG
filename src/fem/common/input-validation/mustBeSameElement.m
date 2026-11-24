function mustBeSameElement(basis, mesh)
%MUSTBESAMEELEMENT  Validate that a basis and a mesh have the same element type.
%   Only bases supported on elements are checked: face and scalar bases
%   (e.g. LegendreBasis, ConstantBasis) fit any mesh.
    if basis.GeometricInputDimension == mesh.geometric_dimension ...
            && ~isequal(mesh.ElementType, basis.ElementType)
        eid = 'ElementType:notEqual';
        msg = 'Element type of the mesh and AbstractFEBasis must be the same.';
        error(eid,msg)
    end
end