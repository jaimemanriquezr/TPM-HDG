function tau = S_identity(obj, c)
%S_IDENTITY  Constant stabilization parameter on every element face.
%
%   tau = S_identity(mesh, c)
%
%   Returns c (default 1) times an array of ones of the size of
%   mesh.boundaries, i.e. 3-by-n_elements.
arguments
    obj AbstractGrid
    c (1,1) double = 1.
end

tau = c * (obj.boundaries).^0;
end
