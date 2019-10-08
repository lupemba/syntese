"""
    reproject_dem(input_dem_path, output_dem_path; from_projection="4326+5773", 
                  to_projection="4979", overwrite=false)

Reproject Digital Elevation Model (DEM) from one EPSG projection to another EPSG projection.

# Arguments
- `input_dem_path`: String; input DEM
- `output_dem_path`: String; output DEM
- `from_projection`: String; from EPSG projection, default: 4326+5773 (using EGM96)
- `to_projection`: String; to EPSG projection, default: 4979 (WGS84)
- `overwrite`: Bool; set to `true` overwrite original file.

# Examples:
```jldoctest
julia> f = open("dem_path.txt")
julia> input_dem_path = readlines(f)[1]
julia> tmp = split(input_dem_path, ".")
julia> output_dem_path = tmp[1] * "_reprojected." * tmp[2]
julia> from_projection = "4326+5773"
julia> to_projection = "4979"
julia> reproject_dem(input_dem_path, output_dem_path, from_projection, to_projection)
Process(`gdalwarp input_dem_path/srtm_38_01.tif output_dem_path/srtm_38_01_reprojected.tif -s_srs EPSG:4326+5773 
-t_srs EPSG:4979`, ProcessExited(0))
```

# Returns
- None, new file created in output_dem_path
"""
function reproject_dem(input_dem_path, output_dem_path; from_projection="4326+5773", to_projection="4979", overwrite=false)
    if overwrite
        call = `gdalwarp -overwrite $input_dem_path $output_dem_path -s_srs EPSG:$from_projection -t_srs EPSG:$to_projection`
    else
        call = `gdalwarp $input_dem_path $output_dem_path -s_srs EPSG:$from_projection -t_srs EPSG:$to_projection`
    end
    run(call)
end