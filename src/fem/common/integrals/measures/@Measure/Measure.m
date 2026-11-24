classdef Measure
    %MEASURE  Integration measure of the form language.
    %
    %   dmu = Measure(domain, type, weight, order) integrates over the
    %   AbstractGrid domain. type is a string (default 'dx'); the
    %   quadrature evaluator recognizes
    %     'dx'   the elements, weighted by element area
    %     'dS'   the boundary of each element, per element side
    %     'ds'   the mesh faces, each seen from both neighbouring elements
    %            (the property faces restricts it to a subset of faces;
    %            empty means all faces)
    %     'de'   the segments of a one-dimensional mesh, weighted by length
    %   weight (default 1.0) multiplies the measure and order (default [],
    %   chosen from the integrand) sets the quadrature degree. All
    %   arguments are optional.
    %
    %   dmu = Measure(preset = "dx" | "ds" | "dS" | "de") sets the type
    %   alone and leaves the domain empty; the resulting measure is bound
    %   to a mesh later, in the form that uses it.
    %
    %   Use  form = integrand * dmu  to build an integral. A test
    %   expression without data (v, c*v, w*v) integrates as a linear form.
    %
    %See also: Datum, FEField
    properties
        domain AbstractGrid = SimplicialGrid.empty;
        type string
        weight double = 1.0;
        order double = [];
        faces double = [];      % facet subset for 'ds' (empty = all faces)
    end

    methods
        function obj = Measure(domain, type, weight, order, options)
            arguments
                domain AbstractGrid = SimplicialGrid.empty
                type string = 'dx';
                weight double = 1.0;
                order double = []
                options.preset = string.empty;
            end

            % a preset (dx, ds, dS, de) sets the measure type alone; without
            % one the positional arguments set type, domain, weight and order
            switch options.preset
                case 'dx'
                    obj.type = 'dx';
                case 'ds'
                    obj.type = 'ds';
                case 'dS'
                    obj.type = 'dS';
                case 'de'
                    obj.type = 'de';
                otherwise
                    obj.type = type;
                    obj.domain = domain;
                    obj.weight = weight;
                    obj.order = order;
            end
        end
        
        function I = mtimes(obj, integrand)
            if isa(integrand, "AbstractFEFunction") && ~isa(integrand, "LinearForm") ...
                    && contains(integrand.signature, "Test") ...
                    && ~contains(integrand.signature, "Trial")
                % a test expression without data (v, c*v, w*v for a
                % weight w) is the linear form of itself
                integrand = LinearForm(@(t, varargin) t.evaluate(varargin{:}), integrand);
            end
            if isempty(obj.domain)
                I = IntegralForm(integrand, obj);
            else
                if isa(integrand, "LinearForm")
                    I = integrate(integrand, obj.domain, ...
                                    measure_type=obj.type, ...
                                    weight=obj.weight, ...
                                    order=obj.order, ...
                                    faces=obj.faces);
                else
                    I = integrate(integrand, Datum(1.0), obj.domain, ...
                                    measure_type=obj.type, ...
                                    weight=obj.weight, ...
                                    order=obj.order, ...
                                    faces=obj.faces);
                end
            end
        end
    end
end