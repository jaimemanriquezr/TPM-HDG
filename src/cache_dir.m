function path = cache_dir()
%CACHE_DIR Directory for generated artefacts, created on first use.
%   PATH = CACHE_DIR() returns the absolute path of the toolbox's cache/
%   directory, creating it if it does not exist. When the toolbox folder
%   is not writable (an installed .mltbx or MPM package), it returns
%   tpmhdg-cache under tempdir instead.
%
%   For generated files such as the symbolic manufactured-solution
%   data, which is derived on every call (FluidStructureCoupling writes
%   a private scratch file there and deletes it after loading). The
%   contents are reproducible; the directory is gitignored.
%
%   Example:
%       filename = cache_dir() + "/stokes_solutions.txt";
%
%   See also: FluidStructureCoupling.get_manufactured_solution
path = fullfile(fileparts(fileparts(mfilename("fullpath"))), "cache");
if ~is_writable_folder(path)
    path = fullfile(tempdir, "tpmhdg-cache");
    if ~isfolder(path)
        mkdir(path);
    end
end
path = string(path);
end

function ok = is_writable_folder(path)
% true if path exists, or can be created, and is writable by the user
ok = isfolder(path);
if ~ok
    [ok, ~, ~] = mkdir(path);   % message outputs keep a failure silent
end
if ok
    [~, attributes] = fileattrib(path);
    ok = attributes.UserWrite;
end
end
