function dmu = measure(mesh, type, options)
%MEASURE  Measure on a mesh selected by name.
%
%   dmu = measure(mesh, type, weight=1)
%
%   type is "element:interior" (default, volume_measure), "element:boundary"
%   or "element:boundaries" (boundary_measure) or "element:face" /
%   "element:faces" (the 'ds' Measure with weight 1). The weight option
%   applies to the first two.
    arguments
        mesh AbstractGrid
        type string = "element:interior";
        options.weight = 1.0;
    end
    switch type
        case "element:interior"
            dmu = volume_measure(mesh, options.weight);
        case {"element:boundaries", "element:boundary"}
            dmu = boundary_measure(mesh, options.weight);
        case {"element:face", "element:faces"}
            dmu = Measure(mesh, "ds", 1.0);
    end
end