classdef (InferiorClasses = {?TrialFunction, ?TestFunction, ?FEField, ...
        ?FEFunction, ?LinearCombination, ?Product, ?Datum, ?Measure, ...
        ?MeshFunction, ?ExteriorNormal, ?Stabilization}) ...
        InterfaceNormal
    % InterfaceNormal(mesh): the outward interface normal of the side
    % owning `mesh`, evaluated at the patch quadrature points. Only
    % meaningful inside dot(...) in a transmission form; see
    % ExtendedFunction.
    properties
        mesh SimplicialGrid
    end
    methods
        function obj = InterfaceNormal(mesh)
            obj.mesh = mesh;
        end

        function c = mtimes(a, b)
            % promote raw trials on contact (phat2*n2_)
            if isa(a, "InterfaceNormal")
                c = mtimes(ExtendedFunction.coerce(b), a);
            else
                c = mtimes(ExtendedFunction.coerce(a), b);
            end
        end

        function c = dot(a, b)
            if isa(a, "InterfaceNormal")
                c = dot(ExtendedFunction.coerce(b), a);
            else
                c = dot(ExtendedFunction.coerce(a), b);
            end
        end
    end
end
