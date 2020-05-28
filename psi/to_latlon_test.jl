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

# padding=2
# index = [1100,1100]
# view =[(index[1]-padding):(index[1]+padding), (index[2]-padding):(index[2]+padding)]
view = [2000:6000, 1000:2000]
gamma_view = [1343:5550, 1:1500]  # I think the sample range is arbitrary

# view = [1:3*1524, 1:1500]
mosaic_view = gamma_view # SlcUtil.get_mosaic_view(m_meta, gamma_view)

footprint = SlcUtil.footprint(m_meta, view)
latlon_window = ((minimum(footprint[1]),maximum(footprint[1])),(minimum(footprint[2]),maximum(footprint[2])))
#latlon_window = ((56.22, 56.52), (8.45455059726242, 8.95))
dem = Load.dem(dem_path, latlon_window; nan_fill= 0, padding=[90,90]);




# input
master_view = gamma_view #[2000-1524:6000-1524, 4801:7801]
meta = (m_meta, s_meta)
precise_orbit = (m_pod, s_pod)
dem = dem
stride = (1,1)


# initilize lut
lut  = Dict{String,Any}()

line = collect(master_view[1])
sample = collect(master_view[2])

# Get master line and sample
lut["master_line"] = line
lut["master_sample"] = sample

master_line,master_sample = Misc.flatten(line,sample)


# Get heights
lat_dem,lon_dem, heights = Misc.flatten(dem...)
line_sample_dem = to_line_sample(hcat(lat_dem,lon_dem),heights,precise_orbit[1]...,meta[1])
lut["heights"] = Misc.interp(line_sample_dem[:,1], line_sample_dem[:,2], heights,
                            master_line, master_sample)
@assert sum(isnan.(lut["heights"])) == 0

# Get latitude and longitude
# lat_lon = to_lat_lon(hcat(master_line, master_sample), lut["heights"], precise_orbit[1]...,meta[1])

## FUNCTION to_latlon STARTS HERE
# function to_lat_lon(line_sample, height, state_vectors, time_state_vectors, meta; c = 299792458)
#lat_lon = to_lat_lon(hcat(master_line, master_sample), lut["heights"], precise_orbit[1]...,meta[1])

 #input:
line_sample = hcat(master_line, master_sample)
height = lut["heights"]
state_vectors = precise_orbit[1][1]
time_state_vectors = precise_orbit[1][2]
meta = meta[1]
c = 299792458

# Hardcode Gamma params:

# put old params in list for comparison:
compare = [m_meta["t_start"], m_meta["t_stop"], m_meta["right_looking"],
           m_meta["incidence_angle_mid"], m_meta["range_sampling_rate"],
           m_meta["azimuth_frequency"], m_meta["slant_range_time"]]

# Notes:
# - meta["incidence_angle_mid"] = 43.5932 differs from 43.5926 (Gamma)-> conclusion: not so relevant, difference is very small
t_gamma_start = 20394.149330
t_gamma_end = 20402.797055
t_start = m_meta["t_start"]
t_stop = m_meta["t_stop"]
sign_angle  = m_meta["right_looking"] ? 1 : -1
theta_0 = sign_angle * 43.5926 * pi / 180 # sign_angle*abs(m_meta["incidence_angle_mid"]*pi/180)
range_pixel_spacing = c / (2 * 6.4345241e+07)   # c/(2*m_meta["range_sampling_rate"])
inv_azimuth_frequency =  1 / 486.4863103  # 1/m_meta["azimuth_frequency"]
r_near =  902747.0461  # m_meta["slant_range_time"]  *c/2 # m

state_vectors_poly, state_vectors_mean, state_vectors_std = Geometry.satellite_trajectory(state_vectors, time_state_vectors, t_start, t_stop)
rad2deg = 180/pi

lat_lon = Array{Float64}(undef,size(line_sample)[1],2)


for i in 1:size(line_sample)[1] #[1, 100, 1000, 2000, 4000]
    time =  t_start + (line_sample[i, 1] - 1) * inv_azimuth_frequency#  (line_sample[i, 1] - 1) * inv_azimuth_frequency
    range_x = r_near + (line_sample[i, 2] - 1) * range_pixel_spacing#  (line_sample[i, 2] - 1) * range_pixel_spacing

    state_vectors_0 = Geometry.polyval_state_vectors(state_vectors_poly, time, state_vectors_mean, state_vectors_std)
    x_sat = state_vectors_0[1:3]
    v_sat = state_vectors_0[4:6]

    x_0 = Geometry.ellipsoid_intersect(x_sat, Geometry.approx_line_of_sight(x_sat, v_sat, theta_0))
    x = Geometry.solve_radar(range_x, height[i], x_0, x_sat, v_sat)
    point =  Geometry.xyx2ellipsiod(x...)
    lat_lon[i,1] = point[1] * rad2deg
    lat_lon[i,2] = point[2] * rad2deg

end


json_path = "/Users/eyu/Google Drive/DTU/10_semester/Persistent_scaterer/phase bug investigation/forEigil20200407/test.json"
gamma_lut_dictionary = Dict()
open(json_path, "r") do f
    global gamma_lut_dictionary
    dicttxt = read(f)  # file information to string
    #print(dicttxt)
    gamma_lut_dictionary = JSON.parse(String(dicttxt))  # parse and transform data
end

lat_grid = reshape(lat_lon[:, 1], (length(line), length(sample)))
lon_grid = reshape(lat_lon[:, 2], (length(line), length(sample)))



# test lat_lon against gamma:
lat_gamma_1 = 56.5034
lon_gamma_1 = 8.6401
lat_gamma_100 = 56.4913 # 56.4918
lon_gamma_100 = 8.6363 # 8.6307
lat_gamma_1000 = 56.3809 # 56.3862
lon_gamma_1000 = 8.6014 # 8.5452
lat_gamma_2000 = 56.2582 # 56.2688
lon_gamma_2000 = 8.5629 # 8.4509
lat_gamma_4000 = 56.0128 # 56.0337
lon_gamma_4000 = 8.4873 # 8.2661


lat_rebel_1 = lat_grid[1,1]
lon_rebel_1 = lon_grid[1,1]
lat_rebel_100 = lat_grid[100,1]
lon_rebel_100 = lon_grid[100,1]
lat_rebel_1000 = lat_grid[1000,1]
lon_rebel_1000 = lon_grid[1000,1]
lat_rebel_2000 = lat_grid[2000,1]
lon_rebel_2000 = lon_grid[2000,1]
lat_rebel_4000 = lat_grid[4000,1]
lon_rebel_4000 = lon_grid[4000,1]

plot([1, 100, 1000, 2000, 4000], [[lat_rebel_1, lat_rebel_100, lat_rebel_1000, lat_rebel_2000, lat_rebel_4000], [lat_gamma_1, lat_gamma_100, lat_gamma_1000, lat_gamma_2000, lat_gamma_4000]], title="Latitude", label = ["Rebel" "Gamma"], lw = 3)
plot([1, 100, 1000, 2000, 4000], [[lon_rebel_1, lon_rebel_100, lon_rebel_1000, lon_rebel_2000, lon_rebel_4000], [lon_gamma_1, lon_gamma_100, lon_gamma_1000, lon_gamma_2000, lon_gamma_4000]], title="Longitude", label = ["Rebel" "Gamma"], lw = 3)


doc = EzXML.readxml(m_meta_path)
meta_dict = XMLDict.xml_dict(doc)

startTime = meta_dict["product"]["adsHeader"]["startTime"]
stopTime = meta_dict["product"]["adsHeader"]["stopTime" ]
# t0 =
# s1_meta["t_start"] = Load._str_date2float(t_start, s1_meta["t_0"])


split_time_string = split(split(startTime, "T")[end], ":")
hour_digits, minute_digits = parse.(Int, split_time_string[1:2])
second_decimals_string = split(split_time_string[3], ".")
sec_digits, milisec_digits, microsec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3], second_decimals_string[2][4:end]])

t_gamma_start = 20394.149330
t_gamma_end = 20402.797055
dt_gamma = 2.0555563e-03

time_midnight = Time(DateTime(split(startTime, "T")[1]))

time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits, microsec_digits)
time_since_midnight = time_now - time_midnight
in_secs = time_since_midnight / Nanosecond(1) * 10^(-9)

l_mosaic_gamma_start = (t_gamma_start - in_secs) / dt_gamma + 1
l_mosaic_gamma_end = (t_gamma_end - in_secs) / dt_gamma + 1

gamma_view = [1343:5550, 1:1500]  # I think the sample range is arbitrary


# sample number from slant range:
range_gamma_end = 961752.5220
range_gamma_start = 902747.0461 # similar to:
range_start = m_meta["slant_range_time"] * c / 2  # equal to gamma
final_range_index = (961752.5220 - m_meta["slant_range_time"] * c / 2) / 2.329562 + 1 # should be 25330 - check!
