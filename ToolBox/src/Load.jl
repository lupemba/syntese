module Load
import PyCall
import EzXML
import XMLDict
import Dates
rasterio = PyCall.pyimport("rasterio")


export s1slc_data, s1slc_ann, pod


"""
    s1slc_data(path,view)

    Load Sentinel 1 SLC tiff using the view. Same as loading the the whole tiff
    as img and then slicing it as img[view]. It just saves memory.
"""
function s1slc_data(path,view)
    dataset = rasterio.open(path)

    # subtract one because array[a:b] in python retuns the (a+1)'th to b'th element
    window = ((view[1].start-1,view[1].stop),(view[2].start-1,view[2].stop))
    return dataset.read(1, window=window)
end




"""
    load_s1slc_ann(path)

    Load annotation of s1 SLC products and

    # Arguments
    - `path::String`: path to annotation

    # Output
    - `s1_meta::Dict`: Dict with relevant meta info
"""
function s1slc_ann(path)
    ## open xml files
    doc = EzXML.readxml(path)
    meta_dict = XMLDict.xml_dict(doc)

    # Create empty dict for results
    s1_meta  = Dict{String,Any}()
    burst_meta  = Dict{String,Any}()
    geolocation_dict  = Dict{String,Any}()

    # get start and stop time
    t_start = meta_dict["product"]["adsHeader"]["startTime"]
    t_stop = meta_dict["product"]["adsHeader"]["stopTime" ]
    s1_meta["t_0"] = _str2date(t_start)
    s1_meta["t_start"] = _str_date2float(t_start,s1_meta["t_0"])
    s1_meta["t_stop"] = _str_date2float(t_stop,s1_meta["t_0"]);

    # get general info
    s1_meta["mission_id"] = meta_dict["product"]["adsHeader"]["missionId"]
    s1_meta["product_type"] = meta_dict["product"]["adsHeader"]["productType"]
    s1_meta["polarisation"] = meta_dict["product"]["adsHeader"]["polarisation"]
    s1_meta["mode"] = meta_dict["product"]["adsHeader"]["mode"]
    s1_meta["swath"] = parse(Int, string(meta_dict["product"]["adsHeader"]["swath"][end]))
    s1_meta["image_number"] = meta_dict["product"]["adsHeader"]["imageNumber"]
    s1_meta["absolute_orbit_number"] = meta_dict["product"]["adsHeader"]["absoluteOrbitNumber"]
    s1_meta["mission_data_id"] = meta_dict["product"]["adsHeader"]["missionDataTakeId"]
    s1_meta["pass"] = meta_dict["product"]["generalAnnotation"]["productInformation"]["pass"]


    # get infor for "imageInformation"
    img_info = meta_dict["product"]["imageAnnotation"]["imageInformation"];
    s1_meta["range_pixel_spacing"] = parse(Float64,img_info["rangePixelSpacing"])
    s1_meta["azimuth_frequency"] = parse(Float64,img_info["azimuthFrequency"])
    s1_meta["slant_range_time"] = parse(Float64,img_info["slantRangeTime"])
    s1_meta["incidence_angle_mid"] = parse(Float64,img_info["incidenceAngleMidSwath"])

    # Get infor about burst
    swath_timing = meta_dict["product"]["swathTiming"]
    s1_meta["lines_per_burst"] =  parse(Int,swath_timing["linesPerBurst"])
    s1_meta["samples_per_burst"] =  parse(Int,swath_timing["samplesPerBurst"])
    s1_meta["burst_count"] = parse(Int,swath_timing["burstList"][:count])

    # Get info for each burst
    burst = swath_timing["burstList"]["burst"]
    burst_meta["burst_times"] = [_str_date2float(elem["azimuthTime"],s1_meta["t_0"]) for elem in burst]
    burst_meta["fist_valid_pixel"] = [parse.(Int,split(elem["firstValidSample"][""])) for elem in burst]
    burst_meta["last_valid_pixel"] = [parse.(Int,split(elem["lastValidSample"][""])) for elem in burst];

    # create a array with info about what line the first line in each burst corrosponds to in a mosaic
    first_line_mosaic = 1 .+(burst_meta["burst_times"] .- s1_meta["t_start"]) .*s1_meta["azimuth_frequency"]
    #println("test",first_line_mosaic) = test[1.0, 1343.0, 2684.0, 4027.0, 5368.0, 6710.0, 8053.0, 9393.0, 10735.0, 12078.0]
    burst_meta["first_line_mosaic"] = round.(Int,first_line_mosaic)


    # Get GeoLocations
    geolocation_list = meta_dict["product"]["geolocationGrid"]["geolocationGridPointList"]["geolocationGridPoint"]
    # Add one to line and sample because the annotation is zero index based.
    geolocation_dict["line"] = [parse(Int,elem["line"]) for elem in geolocation_list] .+1
    geolocation_dict["sample"] = [parse(Int,elem["pixel"]) for elem in geolocation_list] .+1
    geolocation_dict["latitude"] = [parse(Float64,elem["latitude"]) for elem in geolocation_list]
    geolocation_dict["longitude"] = [parse(Float64,elem["longitude"]) for elem in geolocation_list]

    # Collect info
    s1_meta["burst_meta"] = burst_meta
    s1_meta["geolocation"] = geolocation_dict


    # Sentinel 1 is right looking
    s1_meta["right_looking"] = true

    return s1_meta
end


"""
    load_pod(path,t_0)

    Load preciese orbit files

    # Arguments
    - `path::String`: path to precise orbits file
    - `t_0::DateTime`: Reference time

    # Output
    - `osv:: Array{float}(Nx6)`: six arrays with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `t_sv::Array{float}(N)`: time of each orbit state relative to t_0 in seconds.
"""
function pod(path,t_0)

    # Load data as dict
    doc = EzXML.readxml(path)
    pod_dict = XMLDict.xml_dict(doc)

    # Acces orbit state vectors
    osv_dict = pod_dict["Earth_Explorer_File"]["Data_Block"]["List_of_OSVs"]["OSV"];

    # get vectors
    tags = ["X","Y","Z","VX","VY","VZ"]
    osv = [[parse(Float64,elem[tag][""]) for elem in osv_dict] for tag in tags];
    osv = hcat(osv[1],osv[2],osv[3],osv[4],osv[5],osv[6])
    # get times
    t_sv = [_str_date2float(elem["UTC"][5:end],t_0) for elem in osv_dict]

    return osv,t_sv
end









"""
    _str2date(time)

Convert a string of format "2017-03-15T05:39:50.703105" to a DateTime with an
accuarcy of minutes
"""
function _str2date(time)
    y = parse(Int,time[1:4])
    m = parse(Int,time[6:7])
    d = parse(Int,time[9:10])
    h = parse(Int,time[12:13])
    mi = parse(Int,time[15:16])

    return Dates.DateTime(y, m, d, h, mi)
end

"""
    _str_date2float(time,t_0)

    # Arguments
    - `time::String`: date of format of format "2017-03-15T05:39:50.703105"
    - `t_0::DateTime`: Reference time

    # Output
    - `dt::Float64`: time relative to t_0 in seconds.
"""
function _str_date2float(time,t_0)
    # find time diff with date time down to min
    t_i = _str2date(time)
    dt_i = Dates.Second(t_i-t_0).value

    # convert seconds seperately to get float pression
    s = parse(Float64,time[18:end])

    return dt_i + s
end

#TODO Add check for row column bounds
function dem(dem_path, latlon_window; nan_fill= NaN, padding=[0,0])

    dem_annotations = rasterio.open(dem_path)
    # find the corresponding slc corner in row, col of the dem and add padding
    (max_row, min_col) = dem_annotations.index(latlon_window[2][1], latlon_window[1][1])
    (min_row, max_col) = dem_annotations.index(latlon_window[2][2], latlon_window[1][2]);

    # make intervals with padding for .read's window function
    row_interval = (min_row - padding[1], max_row + padding[1])
    col_interval = (min_col - padding[2], max_col + padding[2])

    # load subset of dem
    dem_data = dem_annotations.read(1, window=(row_interval, col_interval))

    dem_data[dem_data .== (-32768)] .= nan_fill;

    #dem_view = [row_interval[1]:row_interval[2], col_interval[1]:col_interval[2]]
    transform = dem_annotations.get_transform()

    rows = collect(1:dem_annotations.height);
    columns = collect(1:dem_annotations.width);
    lon = transform[1] .+ rows .* transform[2];
    lat  = transform[4] .+ columns .* transform[6];

    #dem = dem[index1,index2];
    lat = lat[row_interval[1]:row_interval[2]-1]
    lon = lon[col_interval[1]:col_interval[2]-1];

    # return subset of dem and dem_view
    return lat, lon, dem_data
end





end
