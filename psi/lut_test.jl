# LUT unit test
include("../ToolBox/ToolBox.jl")
using .ToolBox
using .Geometry
using .Load
using .Misc
using .SlcUtil

"""
    look_up_table(master_view,meta,precise_orbit,dem)

    Create a look up table for a given view in the master image

    # Arguments
    - `master_view::Array{range}(2)`: - (line_range,sample_range), The ranges can have step sizes other
                                    than 1. The look up table is computed for all pixels in the view
    - `meta::Array{Dict}(2)`: - (master_meta, slave_meta)
    - `precise_orbit`: - (master_precise_orbit, slave_precise_orbit) as returned by Load.precise_orbit()
    - `time_state_vectors::Array{float}(L)`: time of each orbit state relative to t_0 in seconds.
    - `dem`: (lat, lon, dem_data) as returned by Load.dem()
    - `stride::Tuple(2)`: The stride in line and sample.
    # Output
    - `lut::Dict`: Dict with, [master_line,master_sample,heights,latitude,longitude,slave_line,slave_sample]
"""
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
lat_lon = to_lat_lon(hcat(master_line, master_sample), lut["heights"], precise_orbit[1]...,meta[1])
lut["latitude"] = lat_lon[:,1]
lut["longitude"] = lat_lon[:,2]
@assert sum(isnan.(lut["latitude"])) == 0
@assert sum(isnan.(lut["longitude"])) == 0

# Get slave line and sample
line_sample = to_line_sample(hcat(lut["latitude"],lut["longitude"]),lut["heights"],precise_orbit[2]...,meta[2]);
lut["slave_line"] = line_sample[:,1]
lut["slave_sample"] = line_sample[:,2]
@assert sum(isnan.(lut["slave_line"])) == 0
@assert sum(isnan.(lut["slave_sample"])) == 0

ind = 1
println(lut["master_line"][ind], ", ",lut["master_sample"][ind], ", ",lut["heights"][ind], ", ",
      lut["latitude"][ind], ", ",lut["longitude"][ind], ", ",lut["slave_line"][ind], ", ",
      lut["slave_sample"][ind])

m_data = Load.slc_data(s_data_path, master_view)

# m_data_flat = Misc.flatten(m_data...)

_, _, m_imag = Misc.flatten(line, sample, imag(m_data))
_, _, m_real = Misc.flatten(line, sample, real(m_data))

_ = 0

(lon, lat) = 56.025409041065025, 7.995479786469652
(lon, lat) = 56.498503489908686, 8.568195522233822
