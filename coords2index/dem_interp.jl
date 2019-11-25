using Plots
using PyCall
using Images
using Statistics
scipy_interp = pyimport("scipy.interpolate");

include("load_s1slc_ann.jl");
# TODO: virker ikke med x_step > 1
function dem_interp(x_value, y_value, z_value, x_grid_view, y_grid_view; x_step=1, y_step=1)
    # create grid to evaluate interpolation to
    x_tmp = [convert(Int, i) for i in collect(x_grid_view.start:x_step:x_grid_view.stop)]
    y_tmp = [convert(Int, i) for i in collect(y_grid_view.start:y_step:y_grid_view.stop)]

    y_grid = Array{Int}(undef, length(y_tmp), 0)
    x_grid = Array{Int}(undef, length(x_tmp), 0)

    for i in 1:length(x_tmp)
        y_grid = hcat(y_grid, y_tmp)
    end

    for i in 1:length(y_tmp)
        x_grid = hcat(x_grid, x_tmp)
    end
    x_grid = x_grid';

    # gather x and y values in one variable
    xy_value = hcat(x_value, y_value)

    # change input float value to work with scipy
    z_value = [convert(Float64, i) for i in z_value];

    # interpolate z values to input grid view
    return scipy_interp.griddata(xy_value, z_value, (x_grid, y_grid), method="linear");
end
