classdef (Abstract) BoundaryConditions < matlab.mixin.Heterogeneous
    %BOUNDARYCONDITIONS  Abstract list of boundary data attached to tagged faces.
    %
    %   bc = BoundaryConditions(values, tags)
    %
    %   values is a Datum array and tags a string array of the same length;
    %   values(i) is imposed on the boundary faces carrying tags(i). Both
    %   default to empty. bc = add(bc, value, tag) appends one entry.
    %
    %   [x_hat, g_hat, H_hat, dof_list, cache] = apply(bc, problem) evaluates
    %   every condition of a (possibly heterogeneous) array on problem and
    %   returns, over the global dof numbering, the strongly imposed values
    %   x_hat, the weak-condition load g_hat and matrix H_hat, and dof_list
    %   with fields essential.weak, essential.strong, natural and free.
    %   Options: S (condensed system, reused by transfer rows) and cache
    %   (transfer blocks from a previous apply on the same problem).
    %   Unsupported subclasses raise BoundaryConditions:unknownClass.
    %
    %See also: StrongEssentialBoundaryConditions, WeakEssentialBoundaryConditions,
    %   ComponentwiseBoundaryConditions, TransferEssentialBoundaryConditions
    properties
        values Datum
        tags string
    end

    methods
        function obj = BoundaryConditions(values, tags)
            arguments
                values = Datum.empty;
                tags = string.empty;
            end
            obj.values = values;
            obj.tags = tags;
        end
        
        function obj = add(obj, value, tag)
            arguments
                obj BoundaryConditions
                value Datum
                tag string
            end
            obj.values = [obj.values, value];
            obj.tags = [obj.tags, tag];
        end
    end

    methods (Sealed)
        % Sealed so that a heterogeneous array of conditions (e.g.
        % [ComponentwiseBoundaryConditions, TransferEssentialBoundaryConditions])
        % can be applied in one call; the loop below dispatches per element.
        function [x_hat, g_hat, H_hat, dof_list, cache] = apply(obj, problem, options)
            % CACHE.blocks{i} holds the transfer blocks of condition i
            % (TransferEssentialBoundaryConditions.evaluate); pass it back
            % as options.cache to reuse them on the same problem.
            arguments
                obj
                problem
                options.S struct = struct.empty      % condensed system of problem (condense_global_system), reused by the transfer rows
                options.cache struct = struct.empty  % transfer blocks per condition, from a previous apply on the same problem
            end
            cache = options.cache;
            if isempty(cache)
                cache = struct("blocks", {repmat({struct.empty}, 1, numel(obj))});
            end
            % system size spans all global supports; essential data can
            % only attach to the skeleton-supported (trace) entries
            L = global_dof_layout(problem.GlobalSpaces, problem.Mesh);
            n_dof = L.n_dof;
            H_hat = sparse(n_dof, n_dof);
            [x_hat, g_hat] = deal(zeros(n_dof, 1));
            
            dof_list = struct('essential', struct('weak', [], 'strong', []), ...
                              'natural', []);
            ess_dofs_strong_list = cell(length(obj), 1);
            ess_dofs_weak_list = cell(length(obj), 1);
            for i = 1:length(obj)
                bc = obj(i);
                switch class(bc)
                    case "StrongEssentialBoundaryConditions"
                        [ess_dofs_strong, ess_val_strong] = bc.evaluate(problem);
                        ess_dofs_strong = cell2mat(cellfun(@(d) d(:), ess_dofs_strong(:), "UniformOutput", false));
                        x_hat(ess_dofs_strong) = vertcat(ess_val_strong{:});
                        ess_dofs_strong_list{i} = ess_dofs_strong;
                    case {"WeakEssentialBoundaryConditions", ...
                          "TransferEssentialBoundaryConditions"}
                        if isa(bc, "TransferEssentialBoundaryConditions")
                            [ess_dofs_weak, ess_val_weak_H, ess_val_weak_g, cache.blocks{i}] = ...
                                bc.evaluate(problem, S = options.S, blocks = cache.blocks{i});
                        else
                            [ess_dofs_weak, ess_val_weak_H, ess_val_weak_g] = bc.evaluate(problem);
                        end
                        ess_dofs_weak_list{i} = ess_dofs_weak;
                        H_hat = H_hat + ess_val_weak_H;
                        g_hat(ess_dofs_weak) = ess_val_weak_g;
                    case "ComponentwiseBoundaryConditions"
                        [sd, sv, nd, nv] = bc.evaluate(problem);
                        x_hat(sd) = sv;
                        % plain assignment: within one object the last spec
                        % wins on a shared face (as the low-level path), and
                        % across objects the natural dofs are disjoint from
                        % the weak-essential ones written above
                        g_hat(nd) = nv;
                        ess_dofs_strong_list{i} = sd;
                        dof_list.natural = [dof_list.natural; nd];
                    otherwise
                        error("BoundaryConditions:unknownClass", "Unsupported boundary condition class %s.", class(bc));
                end
            end
            dof_list.essential.weak = [ess_dofs_weak_list{:}];
            dof_list.essential.strong = [ess_dofs_strong_list{:}];
            dof_list.free = setdiff(1:n_dof, dof_list.essential.strong);
        end
    end

end