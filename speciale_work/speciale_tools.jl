
# Fits with View [8000:10000,1200:6000]; of vejle
crs_transform =[[0.00035  0  9.43]; 
                [0  -0.00025  55.8305]; 
                [ 0.0 0  1]] 

lat_array  = ((1:800).-1) .* crs_transform[2,2] .+ crs_transform[2,3]
lon_array = ((1:850).-1) .* crs_transform[1,1] .+ crs_transform[1,3];



function geo_refferece_data(data_list,lat,lon, dem, precise_orbit, meta, data_view, dims)
    lat_dem, lon_dem, heights = Misc.flatten(dem...)
    heights = Misc.interp(lat_dem, lon_dem, heights, lat,lon);
    geo_lisa = to_line_sample(hcat(lat,lon),heights,precise_orbit...,meta);
    geo_ref = [Misc.resample(data_view, data, geo_lisa[:,1], geo_lisa[:,2]) for data in data_list];
    geo_ref = [reshape(elem,dims) for elem in geo_ref]
end