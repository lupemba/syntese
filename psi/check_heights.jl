# check heights
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

## LOADING GAMMA DEM
json_path = "/Users/eyu/Google Drive/DTU/10_semester/Persistent_scaterer/phase bug investigation/forEigil20200407/srdem_subset.json"
srdem_subset = Dict()
open(json_path, "r") do f
    global srdem_subset
    dicttxt = read(f)  # file information to string
    #print(dicttxt)
    srdem_subset = JSON.parse(String(dicttxt))  # parse and transform data
end

gamma_dem = Array{Float64}(undef, size(srdem_subset)..., size(srdem_subset[1])...)
for i in range(1, size(srdem_subset)...)
    gamma_dem[i, :] = transpose(srdem_subset[i])
end


dem_path = "/Users/eyu/local_data/data/srtm_38_01/srtm_38_01_reprojected.tif"
master_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170408T053951_20170408T054019_005065_008DBC_AEEF.SAFE"
slave_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170420T053952_20170420T054019_005240_0092C6_3820.SAFE"

m_data_path, m_meta_path, m_calibration_path = Load.slc_paths(master_safe_path, "VV", 3);
s_data_path, s_meta_path, s_calibration_path = Load.slc_paths(slave_safe_path, "VV", 3);

m_meta = Load.slc_meta(m_meta_path);
m_pod = Load.precise_orbit(Load.pod_path(m_meta["t_0"], m_meta["mission_id"],
                        "/Users/eyu/local_data/data/phase_bug/POD"), m_meta["t_0"])

doc = EzXML.readxml(s_meta_path)
s_meta_dict = XMLDict.xml_dict(doc)
s_start_time = s_meta_dict["product"]["adsHeader"]["startTime"]


l_mosaic_gamma_start = (20394.149330 - seconds_since_midnight(m_start_time)) / dt_gamma + 1
l_mosaic_gamma_end = (20402.797055 - seconds_since_midnight(m_start_time)) / dt_gamma + 1
# From Gamma LUT
azimuth_start_index = convert(Int, round((20394.763347 - seconds_since_midnight(s_start_time)) * s_meta["azimuth_frequency"])) + 1
azimuth_end_index = convert(Int, round((20403.411072 - seconds_since_midnight(s_start_time)) * s_meta["azimuth_frequency"])) + 1
range_start_index = convert(Int, round((902747.0461 * 2 / c - s_meta["slant_range_time"]) * s_meta["range_sampling_rate"])) + 1
range_end_index = convert(Int, round((961752.5220 * 2 / c - s_meta["slant_range_time"]) * s_meta["range_sampling_rate"])) + 1

gamma_view = [azimuth_start_index:azimuth_end_index, range_start_index:1500]
footprint = SlcUtil.footprint(m_meta, [1000:6000, 1:2000])
latlon_window = ((minimum(footprint[1]),maximum(footprint[1])),(minimum(footprint[2]),maximum(footprint[2])))
dem = Load.dem(dem_path, latlon_window; nan_fill = 39, padding=[90,90]);

# input
master_view = gamma_view # mosaic_view #[2000-1524:6000-1524, 4801:7801]
stride = (1,1)

## FIRST PART OF LUT, INTERPOLATE HEIGHTS
lut  = Dict{String,Any}()

line = collect(master_view[1])
sample = collect(master_view[2])

# Get master line and sample
lut["master_line"] = line
lut["master_sample"] = sample

master_line,master_sample = Misc.flatten(line,sample)

# Get heights
lat_dem, lon_dem, heights = Misc.flatten(dem...)
line_sample_dem = to_line_sample(hcat(lat_dem, lon_dem), heights, m_pod..., m_meta)
lut["heights"] = Misc.interp(line_sample_dem[:,1], line_sample_dem[:,2], heights,
                            master_line, master_sample)
@assert sum(isnan.(lut["heights"])) == 0

## COMPARE TO GAMMA

height_grid = reshape(lut["heights"], (length(line), length(sample)))

SlcUtil.show_img(gamma_dem)
SlcUtil.show_img(height_grid)

maximum(gamma_dem)
maximum(height_grid)
minimum(gamma_dem)
minimum(height_grid)

1/length(gamma_dem) * sum(gamma_dem)
1/length(height_grid) *sum(height_grid)

SlcUtil.show_img(abs.(height_grid .- gamma_dem))

## Conclusions:
# Perhaps Gamma uses 30-meter STRM?
# Not accurate cut-out area so use Gamma STRM for testing.
