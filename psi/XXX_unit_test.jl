include("../ToolBox/ToolBox.jl")
using .ToolBox
using .Geometry
using .Load
using .Misc
#using Colors
dem_path = "/Users/eyu/local_data/data/srtm_38_01/srtm_38_01_reprojected.tif"
master_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170408T053951_20170408T054019_005065_008DBC_AEEF.SAFE"
slave_safe_path = "/Users/eyu/local_data/data/phase_bug/BB/S1B_IW_SLC__1SDV_20170420T053952_20170420T054019_005240_0092C6_3820.SAFE"

m_data_path, m_meta_path, m_calibration_path = Load.slc_paths(master_safe_path, "VV", 3);
s_data_path, s_meta_path, s_calibration_path = Load.slc_paths(slave_safe_path, "VV", 3);

m_meta = Load.slc_meta(m_meta_path);
s_meta = Load.slc_meta(s_meta_path);
meta = (m_meta, s_meta);
m_pod = Load.precise_orbit(Load.pod_path(m_meta["t_0"], m_meta["mission_id"],
                        "/Users/eyu/local_data/data/phase_bug/POD"), m_meta["t_0"])
s_pod = Load.precise_orbit(Load.pod_path(s_meta["t_0"], s_meta["mission_id"],
                        "/Users/eyu/local_data/data/phase_bug/POD"), s_meta["t_0"])

pod = (m_pod, s_pod)

padding=2
index = [1100,1100]
view =[(index[1]-padding):(index[1]+padding), (index[2]-padding):(index[2]+padding)]
mosaic_view = SlcUtil.mosaic_view(m_meta, view)

footprint = SlcUtil.footprint(m_meta, view)
latlon_window = ((minimum(footprint[1]),maximum(footprint[1])),(minimum(footprint[2]),maximum(footprint[2])))
dem = Load.dem(dem_path, latlon_window; nan_fill= 0, padding=[90,90]);

lut = look_up_table(mosaic_view, meta, pod, dem)
















SlcUtil.show_img(Load.slc_data(m_data_path, view))

footprint = SlcUtil.footprint(m_meta, view)
latlon_window = ((minimum(footprint[1]),maximum(footprint[1])),(minimum(footprint[2]),maximum(footprint[2])))
lat_dem,lon_dem, dem = Load.dem(dem_path, latlon_window; nan_fill= 0, padding=[90,90]);
