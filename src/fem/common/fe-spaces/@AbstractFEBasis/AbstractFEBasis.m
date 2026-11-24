classdef (Abstract) AbstractFEBasis < matlab.mixin.Heterogeneous
    % Abstract base of the finite element bases; arrays of different
    % concrete bases can be mixed (heterogeneous arrays).
    %
    % Subclasses supply Order, Dimension and the geometric constants
    % ElementType, GeometricInputDimension (0, 1 or 2: global, facet or
    % element support), GeometricOutputDimension, IsDiscontinuous and
    % IsPiolaTransformed, together with evaluate.
    %
    % Sealed methods (callable on arrays of bases):
    %   [x_local, w] = get_quadrature(bases, dim)   reference quadrature
    %   get_trial_functions(bases)   one TrialFunction per basis, with
    %                                signature suffix "-i" for the i-th
    %   get_test_functions(bases)    likewise, TestFunction
    %   get_field_functions(bases)   one FEField per basis, returned as
    %                                an array in the first output
    %   get_field_functions(bases, id1, id2, ...)   the same for the
    %                                sub-array bases(id_j), one output per id
    %
    % See also: ScalarFEBasis, TensorFEBasis, PEERSBasis
    properties (Abstract)
        Order double
    end

    properties (Hidden)
        id double = [];
        IsCondensed logical = false;
    end

    properties (Dependent)
        % Support kind, inferred from geometry (2D mesh):
        %   "scalar"  gid = 0: one global dof, no geometric support
        %   "face"    gid = 1: facet trace
        %   "element" gid = 2: element-wise
        % The single authority for the dimension -> kind rule.
        Support string
    end

    methods
        varargout = evaluate_derivatives(obj, n, points, mesh, options)

        function s = get.Support(obj)
            switch obj.GeometricInputDimension
                case 0
                    s = "scalar";
                case 1
                    s = "face";
                case 2
                    s = "element";
            end
        end
    end

    methods (Sealed)
        [x_local, w] = get_quadrature(bases, dim)

        function varargout = get_field_functions(obj, id)
            arguments
                obj AbstractFEBasis
            end

            arguments (Repeating)
                id double
            end

            if ~isempty(id)
                varargout = cell(length(id), 1);
                for ix = 1:length(id)
                    varargout{ix} = get_field_functions(obj(id{ix}));
                end
            else
                varargout = cell(length(obj), 1);
                for ix = 1:length(obj)
                    varargout{ix} = FEField(obj(ix));
                end
                varargout{1} = [varargout{:}];
            end
        end

        function varargout = get_trial_functions(obj)
            n_out = length(obj);
            varargout = cell(1, n_out);
            for i = 1:n_out
                varargout{i} = TrialFunction(obj(i), "-" + i);
            end
        end

        function varargout = get_test_functions(obj)
            n_out = length(obj);
            varargout = cell(1, n_out);
            for i = 1:n_out
                varargout{i} = TestFunction(obj(i), "-" + i);
            end
        end
    end
end
