classdef AbstractCurve
    % Root of the curve taxonomy. The TPM machinery (mesh carving,
    % induced partitions, transmission quadrature) consumes curves
    % through the query surface below; both lineages implement it:
    % ParametrizedCurve analytically, LevelSetCurve through the induced
    % arc-length parametrization built by discretize.
    %
    %   t_range                     parameter interval [t0 t1]
    %   is_periodic()               closed curve: the parameter wraps
    %   get_2d_points(t)            (1, :, :) params -> (2, :, :) points
    %   map_points_to_interface(X)  points -> parameter of closest curve point
    %   tangent(t)                  gamma'(t); function-handle property on
    %                               ParametrizedCurve, method on LevelSetCurve
    %   unit_normal(t)              unit normal at parameter values;
    %                               outward on closed curves
    %
    % The shared machinery built on the query surface lives here:
    % map_mesh_to_interface, get_partitions, get_sliver_jacobian, distance.
    methods
        segments = map_mesh_to_interface(obj, mesh, options)
        [partition_1, partition_2] = get_partitions(obj, mesh_1, mesh_2, options)
        J = get_sliver_jacobian(obj, mesh, partition, options)

        function [d, t, P] = distance(obj, X)
            % Distance from the points X (2 x n) to the curve, with the
            % parameter t (1 x n) and position P (2 x n) of the closest
            % curve point.
            X = reshape(X, 2, []);
            t = reshape(map_points_to_interface(obj, X), 1, []);
            P = reshape(get_2d_points(obj, t), 2, []);
            d = vecnorm(X - P);
        end

        % query-surface stubs: subclasses implement (kept as erroring
        % stubs, not Abstract members, because tangent is a property on
        % one lineage and a method on the other)
        function t = map_points_to_interface(obj, points) %#ok<STOUT,INUSD>
            error("AbstractCurve:notImplemented", ...
                "%s must implement map_points_to_interface", class(obj));
        end

        function p = get_2d_points(obj, t) %#ok<STOUT,INUSD>
            error("AbstractCurve:notImplemented", ...
                "%s must implement get_2d_points", class(obj));
        end

        function n = unit_normal(obj, t) %#ok<STOUT,INUSD>
            error("AbstractCurve:notImplemented", ...
                "%s must implement unit_normal", class(obj));
        end

        function tf = is_periodic(obj) %#ok<STOUT>
            error("AbstractCurve:notImplemented", ...
                "%s must implement is_periodic", class(obj));
        end
    end
end
