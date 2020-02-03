module Load
import PyCall
import EzXML
import XMLDict
import Dates
rasterio = PyCall.pyimport("rasterio")


export slc_data, slc_ann, precise_orbit


"""
    slc_data(path, view, satellite="s1")

    Load Sentinel 1 SLC tiff using the view. Same as loading the the whole tiff
    as img and then slicing it as img[view]. It just saves memory.
"""
function slc_data(path, view, satellite="s1")
    if satellite != "s1"
        println("Warning this function is only for sentinel 1 SCL images. Updates will come later")
    end
    dataset = rasterio.open(path)

    # subtract one because array[a:b] in python retuns the (a+1)'th to b'th element
    window = ((view[1].start-1,view[1].stop),(view[2].start-1,view[2].stop))
    return dataset.read(1, window=window)
end

"""
    slc_meta(path, satellite="s1")

    Load annotation of s1 SLC products and

    # Arguments
    - `path::String`: path to annotation

    # Output
    - `s1_meta::Dict`: Dict with relevant meta info
"""
function slc_meta(path, satellite="s1")
    if satellite != "s1"
        println("Warning this function is only for sentinel 1 SCL images. Updates will come later")
    end
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
    s1_meta["t_start"] = _str_date2float(t_start, s1_meta["t_0"])
    s1_meta["t_stop"] = _str_date2float(t_stop, s1_meta["t_0"]);

    # get general info
    s1_meta["mission_id"] = meta_dict["product"]["adsHeader"]["missionId"]
    s1_meta["product_type"] = meta_dict["product"]["adsHeader"]["productType"]
    s1_meta["polarisation"] = meta_dict["product"]["adsHeader"]["polarisation"]
    s1_meta["mode"] = meta_dict["product"]["adsHeader"]["mode"]
    s1_meta["swath"] = parse(Int, string(meta_dict["product"]["adsHeader"]["swath"][end]))
    s1_meta["radar_frequency"] = parse(Float64,
            meta_dict["product"]["generalAnnotation"]["productInformation"]["radarFrequency"])
    s1_meta["range_sampling_rate"] = parse(Float64,
            meta_dict["product"]["generalAnnotation"]["productInformation"]["rangeSamplingRate"])
    s1_meta["image_number"] = meta_dict["product"]["adsHeader"]["imageNumber"]
    s1_meta["absolute_orbit_number"] = meta_dict["product"]["adsHeader"]["absoluteOrbitNumber"]
    s1_meta["mission_data_id"] = meta_dict["product"]["adsHeader"]["missionDataTakeId"]
    s1_meta["pass"] = meta_dict["product"]["generalAnnotation"]["productInformation"]["pass"]
    s1_meta["azimuth_steering_rate"] = parse(Float64, meta_dict["product"]["generalAnnotation"]["productInformation"]["azimuthSteeringRate"])


    # get info for "imageInformation"
    img_info = meta_dict["product"]["imageAnnotation"]["imageInformation"];
    s1_meta["range_pixel_spacing"] = parse(Float64,img_info["rangePixelSpacing"])
    s1_meta["azimuth_frequency"] = parse(Float64,img_info["azimuthFrequency"])
    s1_meta["slant_range_time"] = parse(Float64,img_info["slantRangeTime"])
    s1_meta["incidence_angle_mid"] = parse(Float64,img_info["incidenceAngleMidSwath"])
    s1_meta["azimuth_pixel_spacing"] = parse(Float64,img_info["azimuthPixelSpacing"])
    s1_meta["number_of_samples"] = parse(Int, img_info["numberOfSamples"])

    # Get info about burst
    swath_timing = meta_dict["product"]["swathTiming"]
    s1_meta["lines_per_burst"] =  parse(Int,swath_timing["linesPerBurst"])
    s1_meta["samples_per_burst"] =  parse(Int,swath_timing["samplesPerBurst"])
    s1_meta["burst_count"] = parse(Int,swath_timing["burstList"][:count])

    s1_meta["azimuth_time_interval"] = parse(Float64, meta_dict["product"]["imageAnnotation"]["imageInformation"]["azimuthTimeInterval"])

    # Get info for each burst
    burst = swath_timing["burstList"]["burst"]
    burst_meta["burst_times"] = [_str_date2float(elem["azimuthTime"],s1_meta["t_0"]) for elem in burst]
    burst_meta["fist_valid_pixel"] = [parse.(Int,split(elem["firstValidSample"][""])) for elem in burst]
    burst_meta["last_valid_pixel"] = [parse.(Int,split(elem["lastValidSample"][""])) for elem in burst];
    burst_mid_times = burst_meta["burst_times"] .+ s1_meta["lines_per_burst"]/(2*s1_meta["azimuth_frequency"])

    # A ordered dictionary of all the Dopple Centroid polynomial estimates
    burst_meta["dc_estimate_list"] = meta_dict["product"]["dopplerCentroid"]["dcEstimateList"]["dcEstimate"]

    # select the polynomials and t0's closest to mid burst time
    data_dc_polynomial = Array{Float64, 2}(undef, s1_meta["burst_count"], 3)
    data_dc_t0 = Array{Float64, 1}(undef, s1_meta["burst_count"])
    dc_time_diff = Array{Float64, 1}(undef, length(burst_meta["dc_estimate_list"]))
    for i in 1:s1_meta["burst_count"]
        for j in 1:length(burst_meta["dc_estimate_list"])
            dc_time = _str_date2float(burst_meta["dc_estimate_list"][j]["azimuthTime"], s1_meta["t_0"])
            dc_time_diff[j] = abs(dc_time - burst_mid_times[i])
        end
        best_dc_index = argmin(dc_time_diff)
        data_dc_polynomial[i, :] = [parse(Float64, param) for param in split(meta_dict["product"]["dopplerCentroid"]["dcEstimateList"]["dcEstimate"][best_dc_index]["dataDcPolynomial"][""])]
        data_dc_t0[i] = parse(Float64, meta_dict["product"]["dopplerCentroid"]["dcEstimateList"]["dcEstimate"][best_dc_index]["t0"])
    end
    burst_meta["data_dc_polynomial"] = data_dc_polynomial
    burst_meta["data_dc_t0"] = data_dc_t0

    # A ordered dictionary of all the azimuth fm rate polynomials
    burst_meta["azimuth_fm_rate_list"] = meta_dict["product"]["generalAnnotation"]["azimuthFmRateList"]["azimuthFmRate"]

    # select the polynomials and t0's closest to mid burst time
    azimuth_fm_rate_polynomial = Array{Float64, 2}(undef, s1_meta["burst_count"], 3)
    azimuth_fm_rate_t0 = Array{Float64, 1}(undef, s1_meta["burst_count"])
    fm_time_diff = Array{Float64, 1}(undef, length(burst_meta["azimuth_fm_rate_list"]))
    for i in 1:s1_meta["burst_count"]
        for j in 1:length(burst_meta["azimuth_fm_rate_list"])
            fm_time = Load._str_date2float(burst_meta["azimuth_fm_rate_list"][j]["azimuthTime"], s1_meta["t_0"])
            fm_time_diff[i] = abs(fm_time - burst_mid_times[i])
        end
        best_dc_index = argmin(fm_time_diff)

        azimuth_fm_rate_polynomial[i, :] = [parse(Float64, param) for param in split(meta_dict["product"]["generalAnnotation"]["azimuthFmRateList"]["azimuthFmRate"][best_dc_index]["azimuthFmRatePolynomial"][""])]
        azimuth_fm_rate_t0[i] = parse(Float64, meta_dict["product"]["generalAnnotation"]["azimuthFmRateList"]["azimuthFmRate"][best_dc_index]["t0"])
    end

    burst_meta["azimuth_fm_rate_polynomial"] = azimuth_fm_rate_polynomial
    burst_meta["azimuth_fm_rate_t0"] = azimuth_fm_rate_t0



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
    precise_orbit(path,t_0)

    Load preciese orbit files

    # Arguments
    - `path::String`: path to precise orbits file
    - `t_0::DateTime`: Reference time

    # Output
    - `state_vectors:: Array{float}(Nx6)`: six arrays with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `time_state_vectors::Array{float}(N)`: time of each orbit state relative to t_0 in seconds.
"""
function precise_orbit(path,t_0)

    # Load data as dict
    doc = EzXML.readxml(path)
    pod_dict = XMLDict.xml_dict(doc)

    # Acces orbit state vectors
    state_vectors_dict = pod_dict["Earth_Explorer_File"]["Data_Block"]["List_of_OSVs"]["OSV"];

    # get vectors
    tags = ["X","Y","Z","VX","VY","VZ"]
    state_vectors = [[parse(Float64,elem[tag][""]) for elem in state_vectors_dict] for tag in tags];
    state_vectors = hcat(state_vectors...)
    # get times
    time_state_vectors = [_str_date2float(elem["UTC"][5:end],t_0) for elem in state_vectors_dict]

    return state_vectors,time_state_vectors
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
function _str_date2float(time, t_0)
    # find time diff with date time down to min
    t_i = _str2date(time)
    dt_i = Dates.Second(t_i-t_0).value

    # convert seconds seperately to get float pression
    s = parse(Float64,time[18:end])

    return dt_i + s
end

"""
    dem(path, lat_lon_window; nan_fill= NaN, padding=[0,0],nan_value =-32768)

    Load preciese orbit files

    # Arguments
    - `path::String`: path to precise orbits file
    - `lat_lon_window::((Float64, Float64)(Float64, Float64))`: Window that specify what subset to load ((min_lat,max_lat)(min_lon,max_lon))
    - `nan_fill::Float64`: Value to replace NaN values
    - `padding::[Int, Int]`: [lat_padding,lon_padding] - loads the extra padding (in pixels) to the image
    - `nan_value::Float64`: Value to be replaced with nan_fill

    # Output
    - `lat:: Array{Float64}(k)`: Lattitudes.
    - `lon::Array{Float64}(l)`: Longitudes.
    - `dem_data::Array{Float64}(kxl)`: Heights.
"""
function dem(path, lat_lon_window; nan_fill= NaN, padding=[0,0],nan_value =-32768)

    dem_annotations = rasterio.open(path)
    # find the corresponding slc corner in row, col of the dem and add padding
    (max_row, min_col) = dem_annotations.index(lat_lon_window[2][1], lat_lon_window[1][1])
    (min_row, max_col) = dem_annotations.index(lat_lon_window[2][2], lat_lon_window[1][2]);


    # make intervals with padding for .read's window function
    min_row = min_row - padding[1]
    max_row = max_row + padding[1]
    min_col = min_col - padding[2]
    max_col = max_col + padding[2]


    if min_row < 0
        min_row = 0
        println("Warning min row are not in the picture.")
    end

    if max_row > dem_annotations.height
        max_row = dem_annotations.height
        println("Warning max row are not in the picture.")
    end

    if min_col < 0
        min_col = 0
        println("Warning min column are not in the picture.")
    end

    if max_col > dem_annotations.width
        max_col = dem_annotations.width
        println("Warning max column are not in the picture.")
    end


    row_interval = (min_row , max_row)
    col_interval = (min_col , max_col)

    # load subset of dem
    dem_data = dem_annotations.read(1, window=(row_interval, col_interval))

    dem_data[dem_data .== (nan_value)] .= nan_fill;

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
