classdef TestFunction < FEField
    %TESTFUNCTION  Placeholder field marking the test function of a form.
    %
    %   v = TestFunction(basis, add_to_signature)
    %
    %   An FEField on basis (default DubinerBasis(0)) without dof. Its
    %   signature is "FEField:Test" + add_to_signature + ":" + the basis
    %   GeometricInputDimension + ":" + the basis Order.
    %
    %See also: TrialFunction, FEField
    methods
        function obj = TestFunction(basis, add_to_signature)
            arguments
                basis = DubinerBasis(0)
                add_to_signature string = "";
            end
            obj = obj@FEField(basis);
            obj.signature = "FEField:Test" + add_to_signature ...
                + ":" + basis.GeometricInputDimension ...
                + ":" + basis.Order;
        end

    end
end
