module Misc

using Base.Cartesian
using PyCall
scipy_interp = pyimport("scipy.interpolate");
ndimage = pyimport("scipy.ndimage");



# credit https://github.com/aamini/FastConv.jl
# direct version (do not check if threshold is satisfied)
@generated function fastconv(E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote

        retsize = [size(E)...] + [size(k)...] .- 1
        retsize = tuple(retsize...)
        ret = zeros(T, retsize)

        convn!(ret,E,k)
        return ret

    end
end

# credit https://github.com/aamini/FastConv.jl
# in place helper operation to speedup memory allocations
@generated function convn!(out::Array{T}, E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote
        @inbounds begin
            @nloops $N x E begin
                @nloops $N i k begin
                    (@nref $N out d->(x_d + i_d - 1)) += (@nref $N E x) * (@nref $N k i)
                end
            end
        end
        return out
    end
end


function interp_grid(index1,index2,data_2d,index1_out, index2_out)
    F = scipy_interp.interp2d(collect(index2),collect(index1),data_2d)
    return F(collect(index2_out),collect(index1_out))
end

"""
    resample(view_in, data, index1_out, index2_out,order=1)

    Interpolate 2d complex data to points given by index1_out, index2_out.

    # Arguments
    - `view::Array{UnitRange,1}`: the view of the data
    - `data::Array{Complex,2}: The complex data
    - `index1_out::Array{Float,N}`: Array with the first index of output points
    - `index2_out::Array{Float,N}`: Array with the second index of output points
    - `order::Int`: order of the polynomial interpolation.
"""
function resample(view_in, data, index1_out, index2_out,order=1)
    index = [index1_out.-view_in[1].start, index2_out.-view_in[2].start]
    data_real = ndimage.map_coordinates(real.(data), index, order=order, mode="nearest")
    data_imag = ndimage.map_coordinates(imag.(data), index, order=order, mode="nearest")
    return data_real .+ data_imag.*im
end

"""
    interp(x_point, y_point, value, x, y,method="linear")

    Julia wrapper for scipy.interpolate.griddata
        interpolates from x_point, y_point, to x,y.

"""
function interp(x_point, y_point, value, x, y,method="linear")
    # interpolate z values to input grid view
    return scipy_interp.griddata(hcat(x_point, y_point), value, (x, y), method);
end

"""
    grid(index_1,index_2)

    Creates two 2d grids from two 1d arrays.
"""
function grid(index_1,index_2)
    n_1 = length(index_1)
    n_2 = length(index_2)

    index1_grid =  Array{eltype(index_1)}(undef, n_1, n_2)
    index2_grid =  Array{eltype(index_2)}(undef, n_1, n_2)

    for i in 1:n_2
        index1_grid[:,i] .= index_1
    end

    for j in 1:n_1
        index2_grid[j,:] .= index_2
    end

    return index1_grid, index2_grid
end

"""
    flatten(index_1,index_2;data_2d)

    Gives two arrays with all combinations of index1 and index2.
    If data_2d is given then it is also returned as an array
"""
function flatten(index_1,index_2)
    index_1_array, index_2_array = grid(index_1,index_2)
    index_1_array = reshape(index_1_array, :)
    index_2_array = reshape(index_2_array, :)
    return index_1_array, index_2_array
end

function flatten(index_1,index_2,data_2d)
    index_1_array, index_2_array = flatten(index_1,index_2)
    data_array = reshape(data_2d, :)
    return index_1_array, index_2_array, data_array
end


"""
    reproject_dem(input_dem_path, output_dem_path; from_projection="4326+5773",
                  to_projection="4979", overwrite=false)

Reproject Digital Elevation Model (DEM) from one EPSG projection to another EPSG projection.

# Arguments
- `input_dem_path`: String; input DEM
- `output_dem_path`: String; output DEM
- `from_projection`: String; from EPSG projection, default: 4326+5773 (using EGM96)
- `to_projection`: String; to EPSG projection, default: 4979 (WGS84)
- `overwrite`: Bool; set to `true` overwrite original file.

# Examples:
```jldoctest
julia> f = open("dem_path.txt")
julia> input_dem_path = readlines(f)[1]
julia> tmp = split(input_dem_path, ".")
julia> output_dem_path = tmp[1] * "_reprojected." * tmp[2]
julia> from_projection = "4326+5773"
julia> to_projection = "4979"
julia> reproject_dem(input_dem_path, output_dem_path, from_projection, to_projection)
Process(`gdalwarp input_dem_path/srtm_38_01.tif output_dem_path/srtm_38_01_reprojected.tif -s_srs EPSG:4326+5773
-t_srs EPSG:4979`, ProcessExited(0))
```

# Returns
- None, new file created in output_dem_path
"""
function reproject_dem(input_dem_path, output_dem_path; from_projection="4326+5773", to_projection="4979", overwrite=false)
    if overwrite
        call = `gdalwarp -overwrite $input_dem_path $output_dem_path -s_srs EPSG:$from_projection -t_srs EPSG:$to_projection`
    else
        call = `gdalwarp $input_dem_path $output_dem_path -s_srs EPSG:$from_projection -t_srs EPSG:$to_projection`
    end
    run(call)
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
