classdef (Abstract) AbstractFEFunction < matlab.mixin.Heterogeneous
    % Base class of the symbolic FE expressions (bases, fields, data,
    % products, sums, forms). A subclass implements
    % evaluate(obj, points, mesh, options); the operators (+, -, *, \, /,
    % .^, transpose, dot, trace, deviator) build new AbstractFEFunction
    % objects that evaluate their operands on demand. The signature string
    % carries role tags ("Trial", "Test"), the geometric dimension and an
    % id; read_signature parses it.

    properties
        signature string = 'AbstractFEFunction'
    end

    methods (Abstract)
        varargout = evaluate(obj, varargin)
        % evaluate(obj, points, mesh, options)
    end

    methods
        function obj = uplus(obj)
            obj = +1 * obj;
        end
        function m_obj = uminus(obj)
            m_obj = -1 * obj;
        end

        obj_x = mtimes(obj, x)
        obj_x = mldivide(obj, x)
        obj_x = mrdivide(obj, x)
        obj_x = plus(obj, x)
        obj_x = minus(obj, x)
        obj_x = power(obj, n);
        
        obj_t = transpose(obj);
        obj_t = ctranspose(obj);
        
        obj_dot_x = dot(obj, x)
        function obj_norm = norm(obj)
            obj_norm = dot(obj, obj);
        end
        
        u_proj = project(obj, mesh, basis);

        tr = trace(obj)
        dev = deviator(obj)

        p_bdy = evaluate_boundary(obj, x, mesh, options)

        I = integrate(u, v, mesh, options)

        [s, dim, order, id] = read_signature(obj)
        
        function varargout = value_size(obj, varargin)
            % the value shape of the function (2x1 vector, 2x2 tensor,
            % ...), probed by a dummy evaluation; builtin size() is the
            % object-array size
            [~, dim] = read_signature(obj);
            obj_0 = obj.evaluate(zeros(dim, 1));
            obj_0 = zeros(size(obj_0, [3 4]));
            [varargout{1:nargout}] = size(obj_0, varargin{:});
        end
    end

    methods (Static, Sealed, Access = protected)
       function cobj = convertObject(DomClass,objToConvert)
           switch DomClass
            case "WeightFunction"
                cobj = WeightFunction(objToConvert);
            case "MeshFunction"
                cobj = MeshFunction(objToConvert);
           end
       end
    end
    
    methods (Static)
        function [trial, test] = is_trial_or_test(u, v)
            u_v = [u, v];
            s = [u_v.signature];

            has_trial = contains(s, "Trial");
            if any(has_trial)
                trial = u_v(has_trial);
            else
                trial = FEFunction.empty;
            end

            has_test = contains(s, "Test");
            if any(has_test)
                test = u_v(has_test);
            else
                test = FEFunction.empty;
            end
        end
        
        function signature = match_signature(u, v)
            if isequal(u.signature, v.signature)
                signature = u.signature;
            elseif contains(u.signature, ["Trial", "Test"]);
                signature = u.signature;
            elseif contains(v.signature, ["Trial", "Test"])
                signature = v.signature;
            else
                signature = 'FEFunction';
            end
        end
        
        function u_op_v = apply(op, u, v)
            if ismember(func2str(op), ["times", "mtimes", "dot"])
                [trial, test] = AbstractFEFunction.is_trial_or_test(u, v);
                if ~isempty(trial) && ~isempty(test)
                    % operands as written; the form reads the roles from
                    % the signatures
                    u_op_v = BilinearForm(u, v, op);
                    return
                elseif ~isempty(test)
                    % data x test, either order: a linear form, op applied
                    % as written; the test enters as the form's argument
                    % (so a sum of tests can be split per summand)
                    if contains(u.signature, "Test")
                        eval_form = @(t, varargin) AbstractFEFunction.pagefun(op, ...
                            t.evaluate(varargin{:}), v.evaluate(varargin{:}));
                    else
                        eval_form = @(t, varargin) AbstractFEFunction.pagefun(op, ...
                            u.evaluate(varargin{:}), t.evaluate(varargin{:}));
                    end
                    u_op_v = LinearForm(eval_form, test);
                    return
                elseif isscalar(trial) && isequal(func2str(op), "mtimes")
                    % a scalar dot(a, b) of data times the trial, either
                    % order: the structured product trial*a.'*b, so that on
                    % the element boundary each factor is evaluated by its
                    % own evaluate_boundary; the role is the trial's
                    if contains(u.signature, "Trial")
                        s = v;
                    else
                        s = u;
                    end
                    if isa(s, "Product") && numel(s.factors) == 2 ...
                            && isequal(func2str(s.product), "AbstractFEFunction.dot_") ...
                            && ~contains(s.signature, ["Trial", "Test"])
                        u_op_v = Product([trial, s.factors(1).', s.factors(2)], ...
                            @AbstractFEFunction.mtimes_);
                        u_op_v.signature = trial.signature;
                        return
                    end
                end
            end
            eval_fun = @(varargin) AbstractFEFunction.pagefun(op, u.evaluate(varargin{:}), v.evaluate(varargin{:}));
            u_op_v = FEFunction(eval_fun, AbstractFEFunction.match_signature(u, v));
        end
        
        function c = pagefun(op, varargin)
            switch func2str(op)
                case 'mtimes'
                    c = AbstractFEFunction.mtimes_(varargin{:});
                case 'dot'
                    c = AbstractFEFunction.dot_(varargin{:});
                case 'ctranspose'
                    c = AbstractFEFunction.ctranspose_(varargin{:});
                case 'transpose'
                    c = AbstractFEFunction.transpose_(varargin{:});
                case 'mldivide'
                    c = AbstractFEFunction.mldivide_(varargin{:});
                case 'mrdivide'
                    c = AbstractFEFunction.mrdivide_(varargin{:});
                otherwise
                    c = op(varargin{:});
            end
        end

        function c = mtimes_(a, b)
            a_ = permute(a, [3 4 1 2 5 6 7]);
            b_ = permute(b, [3 4 1 2 5 6 7]);
            c = permute(pagemtimes(a_, b_), [3 4 1 2 5 6 7]);
        end

        function c = dot_(a, b)
            a_ = permute(a, [3 4 1 2 5 6 7]);
            b_ = permute(b, [3 4 1 2 5 6 7]);
            pagedot = @(a, b) sum(a.*b, 1:2);
            c = permute(pagedot(a_, b_), [3 4 1 2 5 6 7]);
        end

        function c = ctranspose_(a)
            a_ = permute(a, [3 4 1 2 5 6 7]);
            c = permute(pagectranspose(a_), [3 4 1 2 5 6 7]);
        end

        function c = transpose_(a)
            a_ = permute(a, [3 4 1 2 5 6 7]);
            c = permute(pagetranspose(a_), [3 4 1 2 5 6 7]);
        end

        function c = mldivide_(a, b)
            a_ = permute(a, [3 4 1 2 5 6 7]);
            b_ = permute(b, [3 4 1 2 5 6 7]);
            c = permute(pagemldivide(a_, b_), [3 4 1 2 5 6 7]);
        end

        function c = mrdivide_(a, b)
            a_ = permute(a, [3 4 1 2 5 6 7]);
            b_ = permute(b, [3 4 1 2 5 6 7]);
            c = permute(pagemrdivide(a_, b_), [3 4 1 2 5 6 7]);
        end
    end
end