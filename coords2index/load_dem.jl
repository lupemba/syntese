include("/Users/eyu/Google Drive/DTU/9_semester/synthesis/code/syntese/coords2index/load_s1slc_ann.jl");

function load_dem(dem_annotations, lat_corners_slc, lon_corners_slc; approx_ellipsoid_value=false, nan_fill=false, lat_padding=0, lon_padding=0)

    lat_interval_slc = (minimum(lat_corners_slc), maximum(lat_corners_slc))
    lon_interval_slc = (minimum(lon_corners_slc), maximum(lon_corners_slc))

    # find the corresponding slc corner in row, col of the dem and add padding
    (max_row, min_col) = dem_annotations.index(lon_interval_slc[1], lat_interval_slc[1])
    (min_row, max_col) = dem_annotations.index(lon_interval_slc[2], lat_interval_slc[2]);

    # make intervals with padding for .read's window function
    row_interval = (min_row - lat_padding, max_row + lat_padding)
    col_interval = (min_col - lon_padding, max_col + lon_padding)

    # load subset of dem
    dem = dem_annotations.read(1, window=(row_interval, col_interval))

    if nan_fill != false
        dem[dem .== (-32768)] .= nan_fill;
    end

    if approx_ellipsoid_value != false
        dem = dem .+ approx_ellipsoid_value
    end

    dem_view = [row_interval[1]:row_interval[2], col_interval[1]:col_interval[2]]

    # return subset of dem and dem_view
    return dem, dem_view
end
