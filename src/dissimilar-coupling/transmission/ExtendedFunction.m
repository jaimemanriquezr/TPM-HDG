classdef (InferiorClasses = {?TrialFunction, ?TestFunction, ?FEField, ...
        ?FEFunction, ?LinearCombination, ?Product, ?Datum, ?Measure, ...
        ?MeshFunction, ?ExteriorNormal, ?Stabilization}) ...
        ExtendedFunction
    % A function extended beyond its mesh by the TPM extension
    % operator E of the element polynomial past the computational
    % boundary; Extrapolate and PathIntegral both apply E. A sum of
    % terms, each a primitive trial function of one subproblem
    % evaluated in one of four ways over the induced patch partition:
    %   "onface"       bare local trial at the patch face points
    %   "extrapolate"  Extrapolate(u): at the interface points x + m(x)
    %   "path"         PathIntegral(u): int_0^1 (W u)(x + t m) . m dt
    %   "trace"        skeleton-supported trial, a direct trace block
    % Terms route by their trial's own mesh (each side's functions come
    % from its subproblem's space); dot(., InterfaceNormal(mesh))
    % contracts vector values with that side's outward patch normal;
    % Stabilization() factors multiply by the patch tau. A form is
    % completed by  expr * mu * Measure(mesh_i, "dI")  where the
    % measure's mesh names the test side.
    properties
        % normal: "" none | "dot" full contraction (vector . n) |
        % "mtimes" matrix action (sigma*n row-wise; scalars promote p*n)
        % op: "" | "trace" (tr applied to the evaluated value, before
        % any matrix weight, as in the c2*trace(sg)*eye(2) reconstruction)
        terms struct = struct("coeff", {}, "kind", {}, "trial", {}, ...
            "weight", {}, "stab", {}, "normal", {}, "op", {})
        test = []          % TestFunction of the row (set by expr*mu)
        measure = []       % Measure(mesh_i, "dI") (set by form*dI)
    end

    methods
        function obj = ExtendedFunction(x)
            arguments
                x = []
            end
            if ~isempty(x)
                obj.terms = ExtendedFunction.to_terms(x);
            end
        end

        function c = plus(a, b)
            c = ExtendedFunction.coerce(a);
            b = ExtendedFunction.coerce(b);
            c.terms = [c.terms, b.terms];
        end

        function c = minus(a, b)
            c = plus(a, uminus(ExtendedFunction.coerce(b)));
        end

        function obj = uminus(obj)
            obj = ExtendedFunction.coerce(obj);
            for i = 1:numel(obj.terms)
                obj.terms(i).coeff = -obj.terms(i).coeff;
            end
        end

        function c = mtimes(a, b)
            if isa(b, "ExtendedFunction") && ~isa(a, "ExtendedFunction")
                c = mtimes(b, a);
                return
            end
            c = ExtendedFunction.coerce(a);
            if isa(b, "double") && isscalar(b)
                for i = 1:numel(c.terms)
                    c.terms(i).coeff = b * c.terms(i).coeff;
                end
            elseif isa(b, "TestFunction")
                assert(isempty(c.test), "ExtendedFunction:test", ...
                    "the form already carries a test function");
                c.test = b;
            elseif isa(b, "Measure")
                assert(b.type == "dI", "ExtendedFunction:measure", ...
                    "transmission forms take a Measure(mesh, ""dI"")");
                assert(~isempty(b.domain), "ExtendedFunction:measure", ...
                    "the dI measure must carry the test side's mesh");
                c.measure = b;
            elseif isa(b, "Stabilization")
                for i = 1:numel(c.terms)
                    c.terms(i).stab = true;
                end
            elseif isa(b, "Datum")
                error("ExtendedFunction:mtimes", ...
                    "add Datum sources into the sum, do not multiply");
            elseif isa(b, "InterfaceNormal")
                % matrix action with the side's outward normal:
                % tensor*n contracts the second value index (sigma n);
                % scalar*n promotes to the vector p n
                c = mark_normal(c, b, "mtimes");
            else
                error("ExtendedFunction:mtimes", ...
                    "unsupported factor of class %s", class(b));
            end
        end

        function c = times(a, b)
            c = mtimes(a, b);
        end

        function c = dot(a, b)
            if isa(b, "InterfaceNormal")
                c = mark_normal(ExtendedFunction.coerce(a), b, "dot");
            elseif isa(a, "InterfaceNormal")
                c = mark_normal(ExtendedFunction.coerce(b), a, "dot");
            else
                error("ExtendedFunction:dot", ...
                    "dot on a transmission expression expects an InterfaceNormal");
            end
        end
    end

    methods (Access = private)
        function c = mark_normal(c, nrm, mode)
            for i = 1:numel(c.terms)
                assert(isequal(c.terms(i).trial.mesh, nrm.mesh), ...
                    "ExtendedFunction:normal", ...
                    "InterfaceNormal mesh does not match the term's side");
                assert(c.terms(i).normal == "", "ExtendedFunction:normal", ...
                    "term already carries a normal contraction");
                c.terms(i).normal = mode;
            end
        end
    end

    methods (Static)
        function [cs, st, W] = parse_weight(c)
            % decompose a coefficient node into (scalar, stab flag,
            % matrix action): double | Stabilization | ConstantWeight |
            % Products thereof
            cs = 1.0; st = false; W = [];
            switch class(c)
                case "double"
                    cs = c;
                case "Stabilization"
                    st = true;
                    % scalar coefficients fold into element_weights
                    % (WeightFunction.mtimes double case)
                    assert(isscalar(c.element_weights), ...
                        "ExtendedFunction:weight", ...
                        "non-scalar Stabilization element_weights");
                    cs = c.element_weights;
                case "ConstantWeight"
                    if isscalar(c.weight)
                        cs = c.weight;
                    else
                        W = c.weight;
                    end
                case "Product"
                    for i = 1:numel(c.factors)
                        [ci, si, Wi] = ExtendedFunction.parse_weight(c.factors(i));
                        cs = cs * ci; st = st || si;
                        if isempty(W), W = Wi;
                        elseif ~isempty(Wi), W = W * Wi; end
                    end
                otherwise
                    error("ExtendedFunction:weight", ...
                        "unsupported coefficient of class %s", class(c));
            end
        end

        function obj = coerce(x)
            % always repackage into the base class so sums of
            % Extrapolate/PathIntegral nodes concatenate homogeneously
            if isa(x, "ExtendedFunction")
                obj = ExtendedFunction();
                obj.terms = x.terms;
                obj.test = x.test;
                obj.measure = x.measure;
            else
                obj = ExtendedFunction(x);
            end
        end

        function terms = to_terms(x, coeff, stab, weight, op)
            % flatten a general node into transmission terms; the
            % primitive must remain visible (TrialFunction, or a
            % LinearCombination / Product of them with ConstantWeight or
            % Stabilization factors)
            arguments
                x
                coeff (1,1) double = 1.0
                stab (1,1) logical = false
                weight double = []   % pending left matrix action
                op (1,1) string = ""
            end
            one = @(kind, trial, c, st, W, op) struct("coeff", c, "kind", kind, ...
                "trial", trial, "weight", W, "stab", st, "normal", "", "op", op);
            switch class(x)
                case "TrialFunction"
                    assert(~isempty(x.mesh), "ExtendedFunction:mesh", ...
                        "transmission trials must carry their mesh: " + ...
                        "build them via get_trial_functions(problem)");
                    switch [x.basis.Support]
                        case "face"
                            kind = "trace";
                            assert(isempty(weight), "ExtendedFunction:weight", ...
                                "matrix weights on trace trials are not supported");
                        case "scalar"
                            % coupling-owned scalar (the alpha unknown):
                            % a direct global column of value 1
                            kind = "scalar";
                        otherwise
                            kind = "onface";
                    end
                    terms = one(kind, x, coeff, stab, weight, op);
                case "TracedFunction"
                    assert(op == "", "ExtendedFunction:trace", ...
                        "nested trace operators are not supported");
                    terms = ExtendedFunction.to_terms(x.u, coeff, stab, ...
                        weight, "trace");
                case "FEField"
                    % known trace data (the kinematic history h0): a
                    % skeleton FEField contributes to the RHS
                    assert([x.basis.Support] == "face", ...
                        "ExtendedFunction:data", ...
                        "FEField data terms must live on the skeleton");
                    terms = one("trace_data", x, coeff, stab, weight, op);
                case "Datum"
                    % interface source data: evaluated at the test
                    % side's interface points, moved to the RHS
                    terms = one("data", x, coeff, stab, weight, op);
                case "LinearCombination"
                    terms = struct("coeff", {}, "kind", {}, "trial", {}, ...
                        "weight", {}, "stab", {}, "normal", {}, "op", {});
                    for i = 1:numel(x.summands)
                        % coefficients may themselves be weight nodes
                        % (tau_*(a - b) scales them by the WeightFunction)
                        [cs, st, W] = ExtendedFunction.parse_weight(x.coefficients(i));
                        Wc = weight;
                        if isempty(Wc), Wc = W;
                        elseif ~isempty(W), Wc = Wc * W; end
                        terms = [terms, ExtendedFunction.to_terms( ...
                            x.summands(i), coeff * cs, stab || st, Wc, op)]; %#ok<AGROW>
                    end
                case "Product"
                    fe_part = [];
                    for i = 1:numel(x.factors)   % no arrayfun (size overloads)
                        f_i = x.factors(i);
                        if isa(f_i, "Stabilization") || isa(f_i, "ConstantWeight")
                            % parse_weight also recovers scalars folded
                            % into a weight's element_weights (nu*tau_)
                            [ci, si, Wi] = ExtendedFunction.parse_weight(f_i);
                            coeff = coeff * ci;
                            stab = stab || si;
                            if isempty(weight), weight = Wi;
                            elseif ~isempty(Wi), weight = weight * Wi; end
                        else
                            assert(isempty(fe_part), "ExtendedFunction:product", ...
                                "a product may carry only one FE operand");
                            fe_part = f_i;
                        end
                    end
                    terms = ExtendedFunction.to_terms(fe_part, coeff, stab, weight, op);
                otherwise
                    error("ExtendedFunction:coerce", ...
                        "cannot use a %s in a transmission form " + ...
                        "(the primitive trial must stay visible)", class(x));
            end
        end
    end
end
