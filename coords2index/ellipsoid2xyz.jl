"""
    ellipsoid2xyz(lat::Array{Number,N}, lon::Array{Number,N},
                  height::Array{Number,N}, semi_major_axis=6378137.,
                  flattening=1/298.257223563)

Convert ellipsoid coordinates to cartisian coordinates, based on ellipsoidal
latitudes (rad), longitudes (rad) and heights (m).

# Examples:
```jldoctest
julia> a = 10; e2 = 2; h = [.2, .1, .6];
julia> lat = [pi/4, 3*pi, pi/6]; lon = [pi/2, pi/4, pi/8];
julia> x, y, z = ellipsoid2xyz(lat, lon, h, a, e2);
julia> x
3-element Array{Float64,1}:
  2.90566636126016e-8
 -7.14177848998413
 11.795229079383331
julia> y
3-element Array{Float64,1}:
  4.7453132826267916e8
 -7.141778489984129
  4.885743855978092
julia> z
3-element Array{Float64,1}:
 -4.7453132797983634e8
 -3.637200993467639e-15
 -6.771067811865474
```
"""
function ellipsoid2xyz(lat, lon, height, semi_major_axis=6378137.,
                       flattening=1/298.257223563)
    #=
    TO-DO
    specify arugment types for lat, lon, height.
    =#
    e2 = flattening * (2 - flattening)
    v = semi_major_axis./sqrt.(1 .-e2*sin.(lat).*sin.(lat))
    x=(v+height).*cos.(lat).*cos.(lon)
    y=(v+height).*cos.(lat).*sin.(lon)
    z=(v.*(1-e2)+height).*sin.(lat)
    return x, y, z
end
