import EzXML
import XMLDict
include("convert_time.jl")

"""
    load_s1slc_ann(path)

    Load annotation of s1 SLC products and

    # Arguments
    - `path::String`: path to annotation

    # Output
    - `s1_meta::Dict`: Dict with relevant meta info
"""
function load_s1slc_ann(path)
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
    s1_meta["t_0"] = str2date(t_start)
    s1_meta["t_start"] = str_date2float(t_start,s1_meta["t_0"])
    s1_meta["t_stop"] = str_date2float(t_stop,s1_meta["t_0"]);
    
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

    # Get infor about burst
    swath_timing = meta_dict["product"]["swathTiming"]
    s1_meta["lines_per_burst"] =  parse(Int,swath_timing["linesPerBurst"])
    s1_meta["samples_per_burst"] =  parse(Int,swath_timing["samplesPerBurst"])
    s1_meta["burst_count"] = parse(Int,swath_timing["burstList"][:count])

    # Get info for each burst
    burst = swath_timing["burstList"]["burst"]
    burst_meta["burst_times"] = [str_date2float(elem["azimuthTime"],s1_meta["t_0"]) for elem in burst]
    burst_meta["fist_valid_pixel"] = [parse.(Int,split(elem["firstValidSample"][""])) for elem in burst]
    burst_meta["last_valid_pixel"] = [parse.(Int,split(elem["lastValidSample"][""])) for elem in burst];

    # create a array with info about what line the first line in each burst corrosponds to in a mosaic
    first_line_mosaic = 1 .+(burst_meta["burst_times"] .- s1_meta["t_start"]) .*s1_meta["azimuth_frequency"]
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

    return s1_meta
end
