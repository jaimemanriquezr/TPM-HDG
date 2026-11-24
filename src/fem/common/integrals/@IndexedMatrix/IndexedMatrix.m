classdef IndexedMatrix < matlab.mixin.indexing.RedefinesParen & ...
                         matlab.mixin.indexing.RedefinesBrace
    % A block matrix: a list of IndexedArray blocks (Arrays), each tagged
    % with its [trial test] id. double(A) places every block at the
    % cumulative offsets of its ids and returns one paged array; paren
    % indexing acts on that array, brace indexing selects blocks. Sums
    % add the blocks of equal id (see match) and keep the others apart.
    properties
        Arrays (:, 1) cell
    end

    properties (Dependent)
        id
    end
    
    methods 
        function obj = IndexedMatrix(arrays)
            arguments (Repeating)
                arrays
            end
            if length(arrays) > 1
                obj.Arrays = arrays(:);
            elseif isscalar(arrays)
                arrays = arrays{:};
                if isa(arrays, 'cell')
                    obj.Arrays = arrays;
                elseif isa(arrays, 'IndexedArray')
                    obj.Arrays = {arrays};
                else
                    obj.Arrays = {IndexedArray(arrays)};
                end
            else
                obj.Arrays = {IndexedArray()};
            end
        end
        
        function id = get.id(obj)
            id = cellfun(@(A) A.id, obj.Arrays, "UniformOutput", false);
        end
        
        function obj = assemble(obj, varargin); end
    end
    
    methods (Access=protected)
        % https://se.mathworks.com/help/matlab/ref/matlab.mixin.indexing.redefinesparen-class.html
        function varargout = parenReference(obj, indexOp)
            double_obj = double(obj);
            double_obj = double_obj.(indexOp(1));
            if isscalar(indexOp)
                varargout{1} = double_obj;
                return;
            end
            [varargout{1:nargout}] = double_obj.(indexOp(2:end));
        end

        function obj = parenAssign(obj, indexOp, varargin) %#ok<INUSD>
            error("IndexedMatrix:parenAssign", "IndexedMatrix does not support indexed assignment.");
        end

        function n = parenListLength(obj, indexOp, ctx)
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj, indexOp(3:end), ctx);
        end

        function obj = parenDelete(obj, indexOp)
            obj.Arrays.(indexOp) = [];
        end

        function varargout = braceReference(obj, indexOp)
            obj.Arrays = obj.Arrays(indexOp(1).Indices{:});
            if isscalar(indexOp)
                varargout{1} = obj;
                return;
            end
            [varargout{1:nargout}] = obj.(indexOp(2:end));
        end

        function obj = braceAssign(obj, indexOp, varargin)
            if isempty(obj)
                obj = varargin{1};
            end
            if isscalar(indexOp)
                assert(nargin==3);
                rhs = varargin{1};
                obj.Arrays.(indexOp) = rhs.Arrays;
                return;
            end
            [obj.(indexOp(2:end))] = varargin{:};
        end

        function n = braceListLength(obj, indexOp, ctx)
            if numel(indexOp) <= 2
                n = 1;
                return;
            end
            containedObj = obj.(indexOp(1:2));
            n = listLength(containedObj, indexOp(3:end), ctx);
        end
    end
    
    methods (Access=public)
        function out = cat(dim, varargin)
            numCatMats = nargin - 1;
            newMats = {};
            for ix = 1:numCatMats
                if isa(varargin{ix}, 'IndexedMatrix')
                    newMats = [newMats varargin{ix}.Arrays];
                else
                    newMats = [newMats IndexedArray(varargin{ix})];
                end
            end
            out = IndexedMatrix(newMats{:});
        end
        
        function varargout = size(obj, varargin)
            [varargout{1:nargout}] = size(double(obj), varargin{:});
        end

        function TF = isempty(obj)
            % do not materialize double(obj): the stored blocks may have
            % inhomogeneous trailing (element/face) dimensions, which only
            % compatible subsets extracted via braces can be assembled from
            TF = isempty(obj.Arrays) || all(cellfun(@isempty, obj.Arrays));
        end
    end
        methods (Static, Access=public)
            function obj = empty()
                obj = IndexedMatrix({});
            end
        end
end