classdef TrialFunction < FEField
    %TRIALFUNCTION  Placeholder field marking the trial function of a form.
    %
    %   u = TrialFunction(basis, add_to_signature)
    %
    %   An FEField on basis (default DubinerBasis(0)) without dof. Its
    %   signature is "FEField:Trial" + add_to_signature + ":" + the basis
    %   GeometricInputDimension + ":" + the basis Order.
    %
    %See also: TestFunction, FEField
    methods
        function obj = TrialFunction(basis, add_to_signature)
            arguments
                basis = DubinerBasis(0)
                add_to_signature string = "";
            end
            obj = obj@FEField(basis);
            obj.signature = "FEField:Trial" + add_to_signature ...
                            + ":" + basis.GeometricInputDimension ...
                            + ":" + basis.Order;
        end

    end
end