classdef IntegralForm
    % A form integrated over a measure: the pair (integrand, measure),
    % normally created by  form * dmu  with a LinearForm or BilinearForm
    % and a Measure.
    %
    %   I = IntegralForm(integrand, measure)
    %
    % The Hidden property sign (+1 or -1) records the sign of the term.
    % id is the id of the integrand. Arithmetic builds arrays of terms:
    % A + B is [A, B], A - B is [A, -B], and -A flips the sign of every
    % element. assemble(I, mesh) integrates the terms over mesh.
    %
    % See also: Measure, LinearForm, BilinearForm
    properties
        integrand LinearForm
        measure Measure
    end

    properties (Hidden)
        sign = 1.0;
    end

    properties (Dependent)
        id
    end

    methods
        function obj = IntegralForm(integrand, measure)
            obj.integrand = integrand;
            obj.measure = measure;
        end

        function id = get.id(obj)
            id = obj.integrand.id;
        end

        function p_obj = uplus(obj)
            p_obj = obj;
        end

        function m_obj = uminus(obj)
            m_obj = obj;
            for i = 1:numel(m_obj)
                m_obj(i).sign = -m_obj(i).sign;
            end
        end

        function A_plus_B = plus(A, B)
            A_plus_B = [A, B];
        end

        function A_minus_B = minus(A, B)
            A_minus_B = [A, -B];
        end

    end
end
