# LUT unit test
include("../ToolBox/ToolBox.jl")
using .ToolBox
using .Geometry
using .Load
using .Misc
using .SlcUtil
import EzXML
import XMLDict
using Dates
using Plots
using JSON

function seconds_since_midnight(date_time)
    #try
        if typeof(date_time) == DateTime
            date_time_string = Dates.format(date_time, "yyyy-mm-ddTHH:MM:SS")
        elseif typeof(date_time) == String
            date_time_string = date_time
        end
        # Compute time from DateTime input.
        split_time_string = split(split(date_time_string, "T")[end], ":")
        hour_digits, minute_digits = parse.(Int, split_time_string[1:2])
        second_decimals_string = split(split_time_string[3], ".")
        if length(second_decimals_string) == 1
            sec_digits  = parse.(Int, second_decimals_string[1])
            time_now = Time(hour_digits, minute_digits, sec_digits)
        elseif length(second_decimals_string) == 2
            sec_digits, milisec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3]])
            time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits)
        elseif length(second_decimals_string) == 3
            sec_digits, milisec_digits, microsec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3], second_decimals_string[2][4:end]])
            time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits, microsec_digits)
        elseif length(second_decimals_string) == 4
            sec_digits, milisec_digits, microsec_digits, nanosec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3], second_decimals_string[2][4:6], second_decimals_stringsecond_decimals_string[2][6:end]])
            time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits, microsec_digits, nanosec_digits)
        else
            println("Only down to nanaseconds supported")
        end

        # Convert midnight of date to Time object
        time_midnight = Time(DateTime(split(date_time_string, "T")[1]))

        # get integer number of nanoseconds between input time and midnight
        time_since_midnight = time_now - time_midnight
        time_since_midnight = time_since_midnight / Nanosecond(1) * 10^(-9)
        return time_since_midnight
    #catch
        #println("Expected type: DateTime or String with format yyyy-mm-ddTHH:MM:SS")
        #typeassert(date_time, String)
        #throw(TypeError(seconds_since_midnight, "Input should be String or DateTime with format yyyy-mm-ddTHH:MM:SS", String, date_time))
    end
end

## LOADING
doc = EzXML.readxml(m_meta_path)
m_meta_dict = XMLDict.xml_dict(doc)
m_startTime = meta_dict["product"]["adsHeader"]["startTime"]

json_path = "/Users/eyu/Google Drive/DTU/10_semester/Persistent_scaterer/phase bug investigation/forEigil20200407/test.json"
gamma_lut_dictionary = Dict()
open(json_path, "r") do f
    global gamma_lut_dictionary
    dicttxt = read(f)  # file information to string
    #print(dicttxt)
    gamma_lut_dictionary = JSON.parse(String(dicttxt))  # parse and transform data
end

dem_path = "/Users/eyu/local_data/data/srtm_38_01/srtm_38_01.tif"
master_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170408T053951_20170408T054019_005065_008DBC_AEEF.SAFE"
slave_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170420T053952_20170420T054019_005240_0092C6_3820.SAFE"

m_data_path, m_meta_path, m_calibration_path = Load.slc_paths(master_safe_path, "VV", 3);
s_data_path, s_meta_path, s_calibration_path = Load.slc_paths(slave_safe_path, "VV", 3);

m_meta = Load.slc_meta(m_meta_path);
s_meta = Load.slc_meta(s_meta_path);
m_pod = Load.precise_orbit(Load.pod_path(m_meta["t_0"], m_meta["mission_id"],
                        "/Users/eyu/local_data/data/phase_bug/POD"), m_meta["t_0"])
s_pod = Load.precise_orbit(Load.pod_path(s_meta["t_0"], s_meta["mission_id"],
                        "/Users/eyu/local_data/data/phase_bug/POD"), s_meta["t_0"])

azimuth_start_index = convert(Int, round((20394.149330 - seconds_since_midnight(m_startTime)) * m_meta["azimuth_frequency"])) + 1
azimuth_end_index = convert(Int, round((20402.797055 - seconds_since_midnight(m_startTime)) * m_meta["azimuth_frequency"])) + 1
range_start_index = convert(Int, round((902747.0461 * 2 / c - m_meta["slant_range_time"]) * m_meta["range_sampling_rate"])) + 1
range_end_index = convert(Int, round((961752.5220 * 2 / c - m_meta["slant_range_time"]) * m_meta["range_sampling_rate"])) + 1


# padding=2
# index = [1100,1100]
# view =[(index[1]-padding):(index[1]+padding), (index[2]-padding):(index[2]+padding)]
view = [2000:6000, 1000:2000]
gamma_view = [1343:5550, 1:1500]
# view = [1:3*1524, 1:1500]
mosaic_view = gamma_view # SlcUtil.get_mosaic_view(m_meta, view)

footprint = SlcUtil.footprint(m_meta, view)
latlon_window = ((minimum(footprint[1]),maximum(footprint[1])),(minimum(footprint[2]),maximum(footprint[2])))
#latlon_window = ((56.22, 56.52), (8.45455059726242, 8.95))
dem = Load.dem(dem_path, latlon_window; nan_fill= 0, padding=[90,90]);

# input
master_view = gamma_view # mosaic_view #[2000-1524:6000-1524, 4801:7801]
meta = (m_meta, s_meta)
precise_orbit = (m_pod, s_pod)
dem = dem
stride = (1,1)

## initilize lut
lut  = Dict{String,Any}()

line = collect(master_view[1])
sample = collect(master_view[2])

# Get master line and sample
lut["master_line"] = line
lut["master_sample"] = sample

master_line,master_sample = Misc.flatten(line,sample)


# Get heights
lat_dem, lon_dem, heights = Misc.flatten(dem...)
line_sample_dem = to_line_sample(hcat(lat_dem,lon_dem),heights,precise_orbit[1]...,meta[1])
lut["heights"] = Misc.interp(line_sample_dem[:,1], line_sample_dem[:,2], heights,
                            master_line, master_sample)
@assert sum(isnan.(lut["heights"])) == 0

# Get latitude and longitude
lat_lon = to_lat_lon(hcat(master_line, master_sample), lut["heights"], precise_orbit[1]...,meta[1])
lut["latitude"] = lat_lon[:,1]
lut["longitude"] = lat_lon[:,2]
@assert sum(isnan.(lut["latitude"])) == 0
@assert sum(isnan.(lut["longitude"])) == 0

## TEST LAT LON
lat_grid = reshape(lat_lon[:, 1], (length(line), length(sample)))
lon_grid = reshape(lat_lon[:, 2], (length(line), length(sample)))

lon_idx = gamma_lut_dictionary[1]["values"]
lat_idx = gamma_lut_dictionary[2]["values"]

lat_values = Array{Float64}(undef, 100, 1)
lon_values = Array{Float64}(undef, 100, 1)

for i in range(1, 100)
    idx1 = lon_idx[i]
    idx2 = lat_idx[i]
    lat_values[i] = lat_grid[idx2, idx1]
end

for i in range(1, 100)
    idx1 = lon_idx[i]
    idx2 = lat_idx[i]
    lon_values[i] = lon_grid[idx2, idx1]
end

maximum(abs.(lat_values - gamma_lut_dictionary[3]["values"]))
maximum(abs.(lon_values - gamma_lut_dictionary[4]["values"]))

scatter(lat_values, (lat_values - gamma_lut_dictionary[3]["values"]), xlabel = "Latitude [deg]", ylabel = "Error [deg]")
scatter(lon_values, (lon_values - gamma_lut_dictionary[4]["values"]), xlabel = "Longitude [deg]", ylabel = "Error [deg]")

##
# Get slave line and sample
# line_sample = to_line_sample(hcat(lut["latitude"],lut["longitude"]),lut["heights"],precise_orbit[2]...,meta[2]);
# function to_line_sample(
# input:
lat_lon = hcat(lut["latitude"],lut["longitude"])
height = lut["heights"]
state_vectors = precise_orbit[2][1]
time_state_vectors = precise_orbit[2][2]
meta = meta[2]
c = 299792458

t_0 = meta["t_0"]
t_start = meta["t_start"]
t_stop = meta["t_stop"]
compare = [meta["range_sampling_rate"], meta["azimuth_frequency"], meta["slant_range_time"]]

inv_range_pixel_spacing = 2 * 6.4345241e+07 / c  # (2*meta["range_sampling_rate"])/c
azimuth_frequency =  486.4863103 # meta["azimuth_frequency"]
r_near =  902747.0461  # meta["slant_range_time"]  *c/2
deg2rad = pi/180

state_vectors_poly, state_vectors_mean, state_vectors_std = Geometry.satellite_trajectory(state_vectors, time_state_vectors, t_start, t_stop)

line_sample = Array{Float64}(undef,size(lat_lon)[1],2)

## to_line_sample FOR LOOP
# for i in [1,2,3] # 1:size(lat_lon)[1]
for i in 1:size(lat_lon)[1]
    point_xyz = Geometry.ellipsoid2xyz(lat_lon[i,1]* deg2rad, lat_lon[i,2]* deg2rad ,height[i])
    time,range = Geometry.zero_doppler_bisect(point_xyz, t_start, t_stop,
                                    state_vectors_poly, state_vectors_mean, state_vectors_std)
    line_sample[i,1] = 1 + (time - t_start) * azimuth_frequency
    line_sample[i,2] = 1 + (range - r_near) * inv_range_pixel_spacing
end

line_grid = reshape(line_sample[:, 1], (length(line), length(sample)))
sample_grid = reshape(line_sample[:, 2], (length(line), length(sample)))

## Lav en convert to time since midnight funktion:
doc = EzXML.readxml(s_meta_path)
meta_dict = XMLDict.xml_dict(doc)
startTime = meta_dict["product"]["adsHeader"]["startTime"]

##

json_path = "/Users/eyu/Google Drive/DTU/10_semester/Persistent_scaterer/phase bug investigation/forEigil20200407/test.json"
gamma_lut_dictionary = Dict()
open(json_path, "r") do f
    global gamma_lut_dictionary
    dicttxt = read(f)  # file information to string
    #print(dicttxt)
    gamma_lut_dictionary = JSON.parse(String(dicttxt))  # parse and transform data
end


sample_indices = gamma_lut_dictionary[1]["values"]  # sample indices
line_indices = gamma_lut_dictionary[2]["values"]  # sample indices

# line_grid and sample_grid are a meshgrid, all columns or rows are the same

azimuth_times = seconds_since_midnight(startTime) .+ (line_grid[line_indices, 1]) / s_meta["azimuth_frequency"]
slant_range = (s_meta["slant_range_time"] .+ (sample_grid[1, sample_indices]) / s_meta["range_sampling_rate"]) .* c / 2

range_times_1 = (s_meta["slant_range_time"] .+ (sample_grid[1, 1]) / s_meta["range_sampling_rate"]) .* c / 2
range_times_100 = (s_meta["slant_range_time"] .+ (sample_grid[1, 100]) / s_meta["range_sampling_rate"]) .* c / 2
abs(range_times_1 - gamma_range_times[1])
abs(range_times_100 - gamma_range_times[100])

gamma_azimuth_times = gamma_lut_dictionary[5]["values"]
gamma_slant_range = gamma_lut_dictionary[6]["values"]

max_error_azimuth = maximum(abs.(azimuth_times - gamma_azimuth_times))
max_error_slant_range = maximum(abs.(range_times - gamma_range_times))
# sortperm(abs.(slant_range - gamma_slant_range))

# points in subset
scatter(azimuth_times, slant_range, xlabel = "Azimuth time [s]", ylabel = "Slant range, [m]",  markerstrokewidth=0)
scatter!(gamma_azimuth_times, gamma_slant_range, markersize=2.5, markerstrokewidth=0)

sorted_idices = ((line_grid[sort(line_indices), 1])) #sortperm(abs.(azimuth_times - gamma_azimuth_times))
sorted_times = seconds_since_midnight(startTime) .+ sorted_idices / s_meta["azimuth_frequency"]
sorted_gamma_times = gamma_azimuth_times[sortperm(line_indices)]

scatter(sort(line_indices), abs.(sorted_times - sorted_gamma_times), xlabel = "Line", ylabel = "Error, [s]")

sorted_idices = ((line_grid[sort(line_indices), 2])) #sortperm(abs.(azimuth_times - gamma_azimuth_times))
sorted_times = (s_meta["slant_range_time"] .+ (sorted_idices) / s_meta["range_sampling_rate"]) .* c / 2
sorted_gamma_times = gamma_azimuth_times[sortperm(sample_indices)]

scatter(sort(line_indices), abs.(sorted_times - sorted_gamma_times), xlabel = "Azimuth time [s]", ylabel = "Error, [s]")


# error wrt to time and range:
scatter(abs.(azimuth_times - gamma_azimuth_times), xlabel = "Azimuth time [s]", ylabel = "Error, [s]")
scatter(slant_range, abs.(slant_range - gamma_slant_range), xlabel = "Slant range, [m]", ylabel = "Error, [m]")

# error wrt to time and range (indices on x):
scatter(line_grid[line_indices, 1], abs.(azimuth_times - gamma_azimuth_times), xlabel = "Line", ylabel = "Error, [s]")
scatter(line_grid[line_indices, 2], abs.(slant_range - gamma_slant_range), xlabel = "Line", ylabel = "Error, [m]")

s_meta[]

## zero_doppler_bisect WHILE LOOP
point_xyz = Geometry.ellipsoid2xyz(lat_lon[i, 1] * deg2rad, lat_lon[i, 2] * deg2rad, height[i])
# time_, range_ = Geometry.zero_doppler_bisect(point_xyz, t_start, t_stop,
                                # state_vectors_poly, state_vectors_mean, state_vectors_std)
# function zero_doppler_bisect(point_xyz, azimuth_start_time, azimuth_stop_time,
                             # state_vectors_poly, state_vectors_mean, state_vectors_std, debug = 0)
azimuth_start_time = t_start
azimuth_stop_time = t_stop



# Initial bisection parameters, start iteration with entire data interval
search_interval_duration = search_interval_end - search_interval_start
small_time_interval = 1e-6
max_iter = 100
iter = 0
line_of_sight = 0
trial_time = 0
state_vectors = [0.,0.,0.,0.,0.,0.]

# search interval
data_duration = azimuth_stop_time - azimuth_start_time
search_interval_start = (azimuth_start_time - data_duration)
search_interval_end = (azimuth_stop_time + data_duration)

while (search_interval_duration > small_time_interval) && (iter < max_iter)
    global search_interval_end; global search_interval_start;
    search_interval_duration = (search_interval_end - search_interval_start) / 2
    global trial_time = search_interval_start + search_interval_duration

    global state_vectors = Geometry.polyval_state_vectors(state_vectors_poly, trial_time, state_vectors_mean, state_vectors_std)
    trial_sat_position = state_vectors[1:3]
    global line_of_sight = point_xyz .- trial_sat_position
    println(line_of_sight)
    # Compute the squint angle
    trial_sat_velocity = state_vectors[4:6]
    sin_squint_angle = line_of_sight' * trial_sat_velocity

    if (sin_squint_angle < 0)
        search_interval_end = trial_time
    else
        search_interval_start = trial_time
    end
    global iter += 1;
end

if iter >= max_iter
    print("Error, max iteration reached")
end

output = [trial_time,sqrt(line_of_sight' * line_of_sight)]


## BACK

line_sample[i,1] = 1 + (time - t_start) * azimuth_frequency
line_sample[i,2] = 1 + (range - r_near) * inv_range_pixel_spacing
# end

line_sample[1, 1]
line_sample[1, 2]

line_grid = reshape(line_sample[:, 1], (length(line), length(sample)))
sample_grid = reshape(line_sample[:, 2], (length(line), length(sample)))

rebel_line_1 = line_grid[1,1]
rebel_sample_1 = sample_grid[1,1]

gamma_line_1 = 2.21722950540623
gamma_sample_1 = 17.7534334735424
