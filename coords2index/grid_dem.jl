# INPUT
# skal tage dem_view som input, output kun lat lon

function grid_dem(dem, dem_annotations, dem_view)
    #dem = dem_annotations.read(1)
    transform = dem_annotations.get_transform()

    rows = collect(1:dem_annotations.height);
    columns = collect(1:dem_annotations.width);
    lon = transform[1] .+ rows .* transform[2];
    lat  = transform[4] .+ columns .* transform[6];

    #dem = dem[index1,index2];
    lat = lat[dem_view[1][1]:dem_view[1][end]-1]
    lon = lon[dem_view[2][1]:dem_view[2][end]-1];

    lat_grid = Array{Float64}(undef, length(lat), length(lon))
    lon_grid = Array{Float64}(undef,length(lat), length(lon))
    for i = 1:length(lat)
        for j = 1:length(lon)
            lat_grid[i,j] = lat[i]
            lon_grid[i,j] = lon[j]
        end
    end

    # lat_lon_height grid
    return hcat(reshape(lat_grid, :), reshape(lon_grid, :));
end
