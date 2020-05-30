include("../ToolBox/ToolBox.jl")
using .ToolBox
using .Geometry
using .Load
using .Misc
using .UnitTest
using Plots
## LOADING
path_to_gamma_json = "/Users/eyu/Google Drive/DTU/10_semester/Persistent_scaterer/phase bug investigation/forEigil20200407/srdem_subset.json"
gamma_dem = UnitTest.gamma_dem_json(json_path=path_to_gamma_json)

master_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170408T053951_20170408T054019_005065_008DBC_AEEF.SAFE"
m_data_path, m_meta_path, m_calibration_path = Load.slc_paths(master_safe_path, "VV", 3);
m_meta = Load.slc_meta(m_meta_path);
m_pod = Load.precise_orbit(Load.pod_path(m_meta["t_0"], m_meta["mission_id"],
                        "/Users/eyu/local_data/data/phase_bug/POD"), m_meta["t_0"])
m_start_time, m_stop_time = UnitTest.meta_start_datetime(m_meta_path)

# coded like this for now,
# TODO: look into loading gamma .par and get these params out
c = 299792458
gamma_meta = Dict()
gamma_meta["t_start"] = m_meta["t_start"] + 20394.149330 - UnitTest.seconds_since_midnight(m_start_time)
gamma_meta["t_stop"] = m_meta["t_stop"] + 20402.797055 - UnitTest.seconds_since_midnight(m_stop_time)
gamma_meta["right_looking"] = true
gamma_meta["incidence_angle_mid"] = 43.5926
gamma_meta["range_sampling_rate"] = 6.4345241e+07
gamma_meta["azimuth_frequency"] = 486.4863103
gamma_meta["slant_range_time"] = 902747.0461 * 2 / c

# input
master_view = [1:size(gamma_dem)[1], 1:size(gamma_dem)[2]]
stride = (1,1)

## FIRST PART OF LUT, INTERPOLATE HEIGHTS
lut  = Dict{String,Any}()

line = collect(master_view[1])
sample = collect(master_view[2])

# Get master line and sample
master_line, master_sample = Misc.flatten(line, sample)

line_sample = hcat(master_line, master_sample)
heights = vec(gamma_dem)
state_vectors = m_pod[1]
time_state_vectors = m_pod[2]

lat_lon = Geometry.to_lat_lon(line_sample, heights, state_vectors, time_state_vectors, gamma_meta; c = 299792458)

json_path = "/Users/eyu/Google Drive/DTU/10_semester/Persistent_scaterer/phase bug investigation/forEigil20200407/test.json"
gamma_lut_dictionary = Dict()
open(json_path, "r") do f
    global gamma_lut_dictionary
    dicttxt = read(f)  # file information to string
    #print(dicttxt)
    gamma_lut_dictionary = JSON.parse(String(dicttxt))  # parse and transform data
end

gamma_lats = gamma_lut_dictionary[3]["values"]
gamma_lons = gamma_lut_dictionary[4]["values"]

lat_grid = reshape(lat_lon[:, 1], (length(line), length(sample)))
lon_grid = reshape(lat_lon[:, 2], (length(line), length(sample)))

indices_1 = gamma_lut_dictionary[1]["values"]
indices_2 = gamma_lut_dictionary[2]["values"]

latitude_difference = Array{Float64}(undef, 100, 1)
longitude_difference = Array{Float64}(undef, 100, 1)

for i in range(1, stop=100)
    latitude_difference[i] = gamma_lats[i] - lat_grid[indices_2[i], indices_1[i]]
    longitude_difference[i] = gamma_lons[i] - lon_grid[indices_2[i], indices_1[i]]
end

scatter(gamma_lons, longitude_difference)
scatter(gamma_lats, latitude_difference)

maximum(longitude_difference) * cosd(56) * 111120
maximum(latitude_difference) * 111120


scatter(gamma_lons, longitude_difference .* cosd(56) .* 111120)
scatter(gamma_lats, latitude_difference .* 111120)

using JLD

save("longitude_difference.jld", "diff_in_deg", longitude_difference)
save("latitude_difference.jld", "diff_in_deg", latitude_difference)

using StatsPlots

histogram(randn(10))


histogram(latitude_difference)
# histogram(latitude_difference)
