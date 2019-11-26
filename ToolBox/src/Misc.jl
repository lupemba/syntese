module Misc

using PyCall
scipy_interp = pyimport("scipy.interpolate");

# TODO: virker ikke med x_step > 1
# TODO: Initialise whole y.x-grid instead of hcat
function interp2grid(x_point, y_point, value, x_range, y_range)
    # create grid to evaluate interpolation to
    x,y = grid(collect(x_range),collect(y_range))
    x = reshape(x, :)
    y = reshape(y, :)
    # interpolate z values to input grid view
    return scipy_interp.griddata(hcat(x_point, y_point), value, (x, y), method="linear");
end


function grid(index_1,index_2)
    n_1 = length(index_1)
    n_2 = length(index_2)

    index1_grid =  Array{typeof(index_1[1])}(undef, n_1, n_2)
    index2_grid =  Array{typeof(index_2[1])}(undef, n_1, n_2)

    for i in 1:n_2
        index1_grid[:,i] .= index_1
    end

    for j in 1:n_1
        index2_grid[j,:] .= index_2
    end

    return index1_grid, index2_grid
end


function flatten(index_1,index_2,data_2d)
    index_1_array, index_2_array = grid(index_1,index_2)
    index_1_array = reshape(index_1_array, :)
    index_2_array = reshape(index_2_array, :)
    data_array = reshape(data_2d, :)
    return index_1_array, index_2_array, data_array
end

"""
    print2maps_co(lat,lon,name="Corner",color="#FF0000")

    Print lat,lon in format there can be copy pasted to maps.co.
    maps.co is an easy way to plot points on a interactive map.
"""
function print2maps_co(lat,lon,name="Corner",color="#FF0000")
    for i=1:length(lat)
        println(lat[i],",",lon[i],",",name,",",color)
    end
end


end
