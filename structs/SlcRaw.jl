import Statistics
import Images
import PyCall


## Structs
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


function Base.:getindex(img::SlcRaw,I1::UnitRange{Int}, I2::UnitRange{Int})
    data = img.data[I1,I2]

    delta_1 = I1.stop-I1.start
    delta_2 = I2.stop-I2.start

    view_1_start =  (img.view[1].start + I1.start-1)
    view_2_start =  (img.view[2].start + I2.start-1)

    view = [view_1_start:(view_1_start+delta_1),view_2_start:(view_2_start+delta_2)]
    return SlcRaw(view, img.meta, data)
end


## utilitize
function show_img(img,max_quantile=0.98)
    return Images.Gray.(abs.(img)./Statistics.quantile(reshape(abs.(img), :), max_quantile))
end


function show_img(img::SlcRaw,max_quantile=0.98)
    return show_img(img.data,max_quantile)
end



function raw_coords(s1_ann,line,sample)
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


function get_footprint(s1_ann,view)
    line = [view[1].start,view[1].stop,view[1].stop,view[1].start]
    sample = [view[2].start,view[2].start,view[2].stop,view[2].stop]
        return raw_coords(s1_ann,line,sample)
end


function get_footprint(img::SlcRaw)
    return get_footprint(img.meta,img.view)
end




function get_mosaic_view(s1_ann,view)
    start_line = view[1].start
    end_line = view[1].stop
    lines_per_burst = s1_ann["lines_per_burst"]
    mosaic_lines = s1_ann["burst_meta"]["first_line_mosaic"]

    burst = Int(floor((start_line-1)/lines_per_burst))+1  # minus one because of 0 vs 1 index
    line_burst = start_line - (burst-1)*lines_per_burst
    mosaic_start = mosaic_lines[burst]+(line_burst-1)

    # Check the burst after has a line there is before the first burst
    if end_line > (burst*lines_per_burst)
        mosaic_start = minimum([mosaic_start,mosaic_lines[burst+1]])
    end


    burst = Int(floor((end_line-1)/lines_per_burst))+1  # minus one because of 0 vs 1 index
    line_burst = end_line - (burst-1)*lines_per_burst
    mosaic_stop = mosaic_lines[burst]+(line_burst-1)

    # Check the burst before has a line there is after the last burst
    if start_line < ((burst-1)*lines_per_burst)
        mosaic_stop = maximum([mosaic_stop,mosaic_lines[burst-1]+lines_per_burst-1])
    end

    return [mosaic_start:mosaic_stop,view[2]]
end



function get_mosaic_view(img::SlcRaw)
    get_mosaic_view(img.meta,img.view)
end




function print2maps_co(lat,lon,name="Corner",color="#FF0000")
    for i=1:length(lat)
        println(lat[i],",",lon[i],",",name,",",color)
    end
end

# Load funtctions
function load_s1slc_data(path,view)
    rasterio = PyCall.pyimport("rasterio")
    dataset = rasterio.open(path)

    # subtract one because array[a:b] in python retuns the (a+1)'th to b'th element
    window = ((view[1].start-1,view[1].stop),(view[2].start-1,view[2].stop))
    return dataset.read(1, window=window)
end
