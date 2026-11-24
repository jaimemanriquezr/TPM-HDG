function mustBeCompatible(basis2, basis1)
%MUSTBECOMPATIBLE  Validate that two bases have equal geometric dimensions and element type.
%   The element type is compared only for bases supported on elements;
%   face and scalar bases fit any element.
    eid = 'GeometricDimension:notEqual';
    msg = 'Geometric dimensions of AbstractFEBasis objects are not compatible.';
    for geo_dim = ["GeometricInputDimension", "GeometricOutputDimension"]
        assert(isequal(basis1.(geo_dim), basis2.(geo_dim)), eid, msg);
    end

    if basis1.Support == "element" && ~isequal(basis1.ElementType, basis2.ElementType)
        eid = 'ElementType:notEqual';
        msg = 'Element type of AbstractFEBasis objects are not compatible.';
        error(eid,msg)
    end
end