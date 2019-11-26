module SlcUtil
import Statistics
import Images

using PyCall
scipy_interp = pyimport("scipy.interpolate");


export SlcRaw, show_img, original_view, footprint

## Structs

"""
    SlcRaw(view, meta, data)

    # Arguments
    - `view::Array{UnitRange,1}`: view defines the subset of the .tiff in data
    - `meta::Dict`: Annotations as a Dict
    - `data::Array{Complex,2}`: Complex pixel values
"""
struct SlcRaw
    view::Array{UnitRange,1}
    meta::Dict
    data::Array{Complex,2}
end

## oveload
function Base.:show(io::IO, img::SlcRaw)
    represent = "SlcRaw::"
    represent *= " " *  img.meta["mission_id"]
    represent *= "-" * img.meta["mode"]
    represent *= "-" * img.meta["polarisation"]
    represent *= "-" * "Swath" * string(img.meta["swath"])
    represent *= " " * string(img.meta["t_0"])
    represent *= " " * "View:[" * string(img.view[1]) * "," * string(img.view[2]) * "]"
    print(io,represent)
end

function Base.:size(img::SlcRaw)
    return size(img.data)
end

function Base.:size(img::SlcRaw,dim)
    return size(img.data,dim)
end

"""
    img[I1,I2]

    Slice ::SlcRaw like an array.

    !!!!!!!
    The keyword end and "open" slice does not work eg.
    [5:end] and [5:].

    Use instead:
    img[5:size(img,1),5:size(img,2)]
    !!!!!!!!
"""
function Base.:getindex(img::SlcRaw,I1::UnitRange{Int}, I2::UnitRange{Int})
    # slice data
    data = img.data[I1,I2]

    # update view
    delta_1 = I1.stop-I1.start
    delta_2 = I2.stop-I2.start
    view_1_start =  (img.view[1].start + I1.start-1)
    view_2_start =  (img.view[2].start + I2.start-1)
    view = [view_1_start:(view_1_start+delta_1),view_2_start:(view_2_start+delta_2)]

    # return sliced SclRaw
    return SlcRaw(view, img.meta, data)
end


## utilitize

"""
    show_img(img,max_quantile=0.98)

    Scales images and show it using Image.Gray()
"""
function show_img(img,max_quantile=0.98)
    return Images.Gray.(abs.(img)./Statistics.quantile(reshape(abs.(img), :), max_quantile))
end
function show_img(img::SlcRaw,max_quantile=0.98)
    return show_img(img.data,max_quantile)
end


"""
    original_view(s1_ann)

    Get the view of the entire tiff file associated with the annotations
    "s1_ann"
"""
function original_view(s1_ann)
    return [1:(s1_ann["lines_per_burst"]*s1_ann["burst_count"]),1:s1_ann["samples_per_burst"]]
end
function original_view(img::SlcRaw)
    return original_view(img.meta)
end


"""
    _raw_coords(s1_ann,line,sample)

    Gives approximated coordinates for line,sample
"""
function _raw_coords(s1_ann,line,sample)
    scipy_interp = PyCall.pyimport("scipy.interpolate");
    line_goepoints = s1_ann["geolocation"]["line"]
    sample_goepoints = s1_ann["geolocation"]["sample"]
    latitude_goepoints = s1_ann["geolocation"]["latitude"]
    longitude_goepoints = s1_ann["geolocation"]["longitude"]


    point_latitude = scipy_interp.griddata(hcat(line_goepoints,sample_goepoints),
                                            latitude_goepoints, (line,sample), method="linear")
    point_logitude = scipy_interp.griddata(hcat(line_goepoints,sample_goepoints),
                                            longitude_goepoints, (line,sample), method="linear")

    return point_latitude,point_logitude
end

"""
    footprint(s1_ann,view)

    Gives APPROXIMATED coordinates for the footprint based on
    the geolocations in the annotaions.

    !!!
    Theses coordinates can easily be 100 pixels off
    !!!

"""
function footprint(s1_ann,view)
    # Get corners
    line = [view[1].start,view[1].stop,view[1].stop,view[1].start]
    sample = [view[2].start,view[2].start,view[2].stop,view[2].stop]
        return _raw_coords(s1_ann,line,sample)
end
function footprint(img::SlcRaw)
    return footprint(img.meta,img.view)
end


"""
    mosaic_view(s1_ann,view)

    Gives the view in mosaic geometry (same as sali).
    Mosaic geometry is when the lines are counted from the first
    burst and with equal spacing. It correospond to collapsing all overlapping
    lines. Sample range is the same
"""
function mosaic_view(s1_ann,view)
    # get info
    start_line = view[1].start
    end_line = view[1].stop
    lines_per_burst = s1_ann["lines_per_burst"]
    mosaic_lines = s1_ann["burst_meta"]["first_line_mosaic"]

    # Find start line
    burst = Int(floor((start_line-1)/lines_per_burst))+1  # minus one because of 0 vs 1 index
    line_burst = start_line - (burst-1)*lines_per_burst
    mosaic_start = mosaic_lines[burst]+(line_burst-1)
    # Check the burst after has a line there is before the first burst
    if end_line > (burst*lines_per_burst)
        mosaic_start = minimum([mosaic_start,mosaic_lines[burst+1]])
    end

    # find end line
    burst = Int(floor((end_line-1)/lines_per_burst))+1  # minus one because of 0 vs 1 index
    line_burst = end_line - (burst-1)*lines_per_burst
    mosaic_stop = mosaic_lines[burst]+(line_burst-1)
    # Check the burst before has a line there is after the last burst
    if start_line < ((burst-1)*lines_per_burst)
        mosaic_stop = maximum([mosaic_stop,mosaic_lines[burst-1]+lines_per_burst-1])
    end

    return [mosaic_start:mosaic_stop,view[2]]
end
function mosaic_view(img::SlcRaw)
    mosaic_view(img.meta,img.view)
end





end
