module SlcUtil
import Statistics
import Images
import Colors
include("Misc.jl")
import .Misc

using PyCall
scipy_interp = pyimport("scipy.interpolate");


export SlcRaw, show_img, original_view, footprint, phase_ramp





"""
    complex_coherence(master, slave, flat, kernel, view)

    Calculates the complex coherence by multilooking

    # Arguments
    - master::Array{Complex{Float64},2}: The master data mosaic
    - slave::Array{Complex{Float64},2}: The coregistrated slave data mosaic
    - flat::Array{Complex{Float64},2}: The flat earth and topografi phase
    - kernel::Array{Number,2}: Kernel to multilooking
    - view::Array{range}(2): - (line_range,sample_range), of the master data

    # Output
    - complex_coherence::Array{Complex{Float64},2}: Complex coherence
    - master_intensity::Array{Float64,2}: Multilooked master intensity
    - slave_intensity::Array{Complex{Float64},2}: Multilooked slave intensity
    - lines::Range : line range of results
    - sample::Range : sample range of results
"""
function complex_coherence(master, slave, flat, kernel, view)
    # Define relevant image signals
    signal_1 = master .* conj.(slave) .* flat
    signal_2 = abs2.(master)
    signal_3 = abs2.(slave)
    kernel_1 = convert.(eltype(signal_1),kernel)
    kernel_2 = convert.(eltype(signal_2),kernel)
    kernel_3 = convert.(eltype(signal_3),kernel)

    # Compute real and imaginary parts seperately
    interferogram =  Misc.fastconv(signal_1, kernel_1)
    master_intensity = Misc.fastconv(signal_2, kernel_2)
    slave_intensity = Misc.fastconv(signal_3, kernel_3)

    # Compute the complex coherence
    complex_coherence = interferogram ./ (sqrt.(master_intensity .* slave_intensity));

    # Cut away padded areas
    no_padd = (size(kernel, 1):size(master, 1)-size(kernel, 1)+1, size(kernel, 2):size(master, 2)-size(kernel, 2)+1)
    complex_coherence = complex_coherence[no_padd...]
    master_intensity = master_intensity[no_padd...]./ length(kernel_2)
    slave_intensity = slave_intensity[no_padd...]./ length(kernel_3)

    #TODO change lines, samples to view
    # Pixel positions in line, sample
    lines = ((1-(size(kernel)[1]-1)/2) :1: (size(master,1) +(size(kernel)[1]-1)/2)) .+ (view[1].start-1)
    lines = lines[no_padd[1]]
    samples = ((1-(size(kernel)[2]-1)/2) :1: (size(master,2) +(size(kernel)[2]-1)/2)) .+ (view[2].start-1)
    samples = samples[no_padd[2]]

    return complex_coherence, master_intensity, slave_intensity, lines, samples
end


"""
    mosaic(data,view_data,meta)

    Remove the burst structure of the data. Equvilant to ESA snap
    TOPS deburst

    # Arguments
    - data::Array{Complex{Float64},2}: The data
    - view::Array{range}(2): - (line_range,sample_range), of the data
    - `meta::Dict`:  Meta information from the Sentinel-1 SLC image

    # Output
    - mosaic_data::Array{Complex{Float64},2}: mosaiced data
    - mosaic_view::Array{Float64,2}:(line_range,sample_range), view of the mosaiced
                                    data
"""
function mosaic(data,view_data,meta)
# Find the burst in the image
start_burst = ceil(Int,(view_data[1].start)/meta["lines_per_burst"])
end_burst = ceil(Int,(view_data[1].stop-1)/meta["lines_per_burst"])

#find the overlaps
overlap = meta["lines_per_burst"].+meta["burst_meta"]["first_line_mosaic"][1:end-1] .- meta["burst_meta"]["first_line_mosaic"][2:end]
d_ovelap_m = floor.(Int,overlap .*0.5)
d_ovelap_p = overlap .- d_ovelap_m

# Initialise start line
start_line = 1

# Check if the first burst have enough line to be included
lines_in_first_burst = 1 + start_burst* meta["lines_per_burst"]-(view_data[1].start)
view_start_line = meta["burst_meta"]["first_line_mosaic"][start_burst] +(meta["lines_per_burst"]- lines_in_first_burst)
if lines_in_first_burst <= d_ovelap_m[start_burst]
    # If not jump to next burst
    start_burst += 1
    # update valuse
    start_line = lines_in_first_burst +1
    view_start_line = meta["burst_meta"]["first_line_mosaic"][start_burst]
end

# Check if the end burst have enough line to be included
lines_in_end_burst = view_data[1].stop - (end_burst-1)* meta["lines_per_burst"]

if end_burst == 1
    println("To Be made")
end

if lines_in_end_burst <= d_ovelap_p[end_burst-1]
    # if not skip the end burst
    end_burst -= 1
    lines_in_end_burst =  meta["lines_per_burst"]
end
# make end line
view_end_line = meta["burst_meta"]["first_line_mosaic"][end_burst]+lines_in_end_burst -1

# Check for single burst
if end_burst == start_burst
    println("To Be made")
end

# Make the mosaic view
mosaic_view = (view_start_line:view_end_line,view_data[2])

# initialize array for results
mosaic_data = Array{eltype(data)}(undef,length.(mosaic_view)...)

# Find initilize mosaic with the first burst exept a bit of the end part
end_line = meta["lines_per_burst"]*start_burst-d_ovelap_m[start_burst]-view_data[1].start
mosaic_data[1:(end_line-start_line)+1,:] .= data[start_line:end_line,:];

# update
line_temp = (end_line-start_line)+2

# loop throuth the middle burst and cut out the section needed
for k in (start_burst+1):(end_burst-1)
    start_line = 1+d_ovelap_p[k-1]+(k-1)*meta["lines_per_burst"]-view_data[1].start
    end_line = meta["lines_per_burst"]*k-d_ovelap_m[k]-view_data[1].start
    # add the result to mosaic
    mosaic_data[line_temp:line_temp + (end_line-start_line),:] .=  data[start_line:end_line,:]
    line_temp  = line_temp +(end_line-start_line)+1
end

# handle the last burst
start_line = 1+d_ovelap_p[end_burst-1]+(end_burst-1)*meta["lines_per_burst"]-view_data[1].start
end_line = (end_burst-1)*meta["lines_per_burst"]+lines_in_end_burst-(view_data[1].start-1)
mosaic_data[line_temp:line_temp+(end_line-start_line),:] .=  data[start_line:end_line,:];

return mosaic_data,mosaic_view
end




"""
    get_burst_corners(burst_number::Int, meta::Dict)

Computes the corner indices of a burst.

# Arguments
- `burst_number::Int`: The number of the burst of interest.
- `meta::Dict`:  Meta information from the Sentinel-1 SLC image

# Output

# Examples:
```jldoctest
julia> meta = Load.slc_meta(path_meta_1);
julia> burst_number = 1
julia> first_index_row, last_index_row, first_index_col, last_index_col = get_burst_corners(burst_number, meta)
julia> print(first_index_row, ", ",last_index_row, ", ",first_index_col, ", ",last_index_col)
32, 1495, 946, 24512
```
"""
function get_burst_corners(burst_number, meta)
    first_index_row = findfirst(x->x!=-1, meta["burst_meta"]["fist_valid_pixel"][burst_number])
    last_index_row = findlast(x->x!=-1, meta["burst_meta"]["last_valid_pixel"][burst_number])
    first_index_col = meta["burst_meta"]["fist_valid_pixel"][burst_number][first_index_row]
    last_index_col = meta["burst_meta"]["last_valid_pixel"][burst_number][last_index_row]
    return first_index_row, last_index_row, first_index_col, last_index_col
end



"""
    phase_ramp(linesArray{Int,1}, samples::Array{Int,1}, burst_number::Int, meta::Dict, precise_orbit:Dict)

Computes the phase ramp (phi) for the given burst number for input lines and samples.

# Arguments
- `lines::Array{Int,1}`: The (Nx1)-lines of interest. Can also be a view i.e. UnitRange{Int64}
- `samples::Array{Int,1}`: The (Mx1)-samples of interest. Can also be a view i.e. UnitRange{Int64}
- `burst_number::Int`: The number of the burst of interest.
- `meta::Dict`:  Meta information from the Sentinel-1 SLC image
- `v_mid:Dict`: Satellite speed at mid burst time

# Output
- ramp::Array{Float64,1}: Array with phase ramp for chosen lines and samples in list format

# Examples:
```jldoctest
julia> meta = Load.slc_meta(path_meta_1);
julia> precise_orbit = Load.precise_orbit(path_pod_1, meta["t_0"]);
julia> burst_number = 1
julia> ramp = phase_ramp(750:900, 1037:2000, burst_number, meta, precise_orbit);
1525×23567 Array{Float64,1}:
 6.17566e8  6.17565e8  6.17564e8  …  5.89003e8  5.89001e8  5.89e8
```

#Notes
Equation reference: Miranda, 2017: "Definition of the TOPS SLC deramping function for products generated by the S-1 IPF"
"""
function phase_ramp(lines, samples, burst_number, meta, v_mid)
    # fish out constants and parameters
    c = 299792458
    k_psi = meta["azimuth_steering_rate"] * pi/180
    dc_coef = meta["burst_meta"]["data_dc_polynomial"][burst_number, :]
    dc_tau0 = meta["burst_meta"]["data_dc_t0"][burst_number]
    fm_coef = meta["burst_meta"]["azimuth_fm_rate_polynomial"][burst_number, :]
    fm_tau0 = meta["burst_meta"]["azimuth_fm_rate_t0"][burst_number]
    f_c = meta["radar_frequency"]
    lines_per_burst = meta["lines_per_burst"]
    number_of_samples = meta["number_of_samples"]
    Delta_t_s = meta["azimuth_time_interval"]
    Delta_tau_s = 1/meta["range_sampling_rate"]
    tau_0 = meta["slant_range_time"]
    v_s = v_mid

    # Temporary functions, allows different x=tau inputs
    k_a(x, fm_param, x0) = fm_param[1] .+ fm_param[2].*(x .- x0) .+ fm_param[3].*(x .- x0).^2 # Doppler FM rate, Eqn. 11
    f_etac(x, dc_param, x0) = dc_param[1] .+ dc_param[2].*(x .- x0) .+ dc_param[3].*(x .- x0).^2; # Doppler centroid freq Eqn. 13

    tau = tau_0 .+ (samples .- 1) .* Delta_tau_s # Slant range time of ith sample, Eqn. 12

    # Doppler rate equations
    k_s = 2 * v_s/c * f_c * k_psi; # Doppler rate from antenna scanning, Eqn. 4
    alpha = 1 .- k_s ./ k_a(tau, fm_coef, fm_tau0); # conversion factor, Eqn. 3
    k_t = k_s ./ alpha; # Doppler Centroid Rate, Eqn. 2

    # Doppler azimuth time equations
    eta_c = - f_etac(tau, dc_coef, dc_tau0) ./ k_a(tau, fm_coef, fm_tau0); # Beam centre crossing time, Eqn. 7
    tau_mid = tau_0 + number_of_samples/2 * Delta_tau_s

    eta_ref = eta_c .- (- f_etac(tau_mid, dc_coef, dc_tau0)/k_a(tau_mid, fm_coef, fm_tau0)); # Reference time, Eqn. 6
    line_in_burst = lines .- lines_per_burst*(burst_number-1)
    eta = -lines_per_burst/2*Delta_t_s .+ (line_in_burst .- 1/2 ) .* Delta_t_s

    # Compute the phase ramp added the modulation term
    ramp = pi * k_t .* (eta .- eta_ref).^2 .+ 2 * pi .* f_etac(tau, dc_coef, dc_tau0) .* (eta .- eta_ref); # Eqn. 14
    return ramp
end



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
    img1 = abs.(img)./Statistics.quantile(reshape(abs.(img), :), max_quantile)
    img1[img1 .> 1] .= 1
    return Colors.Gray.(img1)
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


"""
    plot_phase(img)

    Plot a color image of the phase of img
    - `img::Array{Complex}`:
"""
function plot_phase(img)
    phase = (angle.(img) .+pi)./(2*pi)

    return Colors.RGB{Float32}.(1 .-phase.^2,4 .*(phase .-  phase.^2),phase.^2)
end


"""
    _phase_colorbar(n=17)

    Plot a color bar that fits  plot_phase(img)
"""
function _phase_colorbar(n=17)
    phase = 0:0.003:1
    display(Colors.RGB{Float32}.(1 .-phase.^2,4 .*(phase .-  phase.^2),phase.^2))
    println(" -π"*" "^n,"-π/2"*" "^n,"0"*" "^n,"π/2"*" "^n,"π")
end



"""
    plot_water(img,water,max_quantile=0.98)

    Plot image with water in blue.
"""
function plot_water(img,water,max_quantile=0.98)
    gray = abs.(img)./Statistics.quantile(reshape(abs.(img), :), max_quantile)
    gray[water].= 0
    gray[gray.>1] .= 0
    return Colors.RGB{Float32}.(gray,gray,gray.+water)
end




"""
    calibrate_slave_data(data, view,lut, calibration_dict , kind = "sigma")

    !!! This function is for coregistered slave data!!!
    !!! For not resampled data use calibrate_data()!!!
    Calibrate the data to backscatter.

    # Arguments
    - data::Array{Complex{Float64},2}: The data complex data or the amplitude.
    - view::Array{range}(2): - (line_range,sample_range), of the data
    - `calibration_dict`: calibration data as returned by Load.slc_calibration()
    - kind::String : Must match a key in calibration_dict

- `lut::Dict`:  Look up table
    # Output
    - calibrated_data::Array{Complex{Float64},2}: mosaiced data
"""
function calibrate_slave_data(data, view,lut, calibration_dict , kind = "sigma")

    master_line, master_sample = Misc.flatten(view...)
    # interpolate slave line and sample
    slave_line = Misc.interp_grid(lut["master_line"] ,lut["master_sample"],
    reshape(lut["slave_line"],(length(lut["master_line"]),length(lut["master_sample"])))
    ,view[1], view[2])
    slave_sample = Misc.interp_grid(lut["master_line"] ,lut["master_sample"],
        reshape(lut["slave_sample"],(length(lut["master_line"]),length(lut["master_sample"])))
        ,view[1], view[2]);
    slave_line = reshape(slave_line,:)
    slave_sample= reshape(slave_sample,:)


    return calibrate_data(data, slave_line, slave_sample, calibration_dict, kind)
end


"""
    calibrate_data(data, lines, samples, calibration_dict, kind = "sigma")

    !!! Do not use for resampled data!!!
    !!! see calibrate_slave_data() instead!!!
    Calibrate the data to backscatter.

    # Arguments
    - data::Array{Complex{Float64},2}: The data complex data or the amplitude.
    - lines::range: - range or Array of lines in the data
    - samples::range: - range or Array of samples in the data
    - `calibration_dict`: calibration data as returned by Load.slc_calibration()
    - kind::String : Must match a key in calibration_dict

    # Output
    - calibrated_data::Array{Complex{Float64},2}: mosaiced data
"""
function calibrate_data(data, lines, samples, calibration_dict, kind = "sigma")
    cal_line, cal_sample = Misc.flatten(calibration_dict["line"],calibration_dict["sample"])
    cal_value = Misc.interp(cal_line,cal_sample, reshape(calibration_dict[kind],:), lines , samples)
    cal_value = reshape(cal_value,size(data));
    return data./cal_value
end


"""
    temporal_filter(images, k)

    Apply a simple temporal sepckle filter. Described in
    "Filtering of Multichannel SAR Images" by Shaun Quegan and Jiong Jiong Yu (2001)
    Equation 14

    # Arguments
    - images::Array{Array{Float64,2},1}: Array of coregistered intensity images
    - k::Int: Size of kernel

    # Output
    - filtered_images::Array{Array{Float64,2},1}: Array of filtered images

"""
function temporal_filter(images, k)
    correction = zeros(size(images[1]));
    kernel = ones(k,k)
    pad = round(Int,(k-1)/2)
    n = length(images)
    sigma = similar(images)


    for i in 1:n

        sigma[i] = Misc.fastconv(images[i],kernel)[1+pad:end-pad,1+pad:end-pad]./k^2;
        correction .+= images[i]./ sigma[i]
    end

    correction = correction./n

    return [elem.*correction for elem in sigma]
end


end
