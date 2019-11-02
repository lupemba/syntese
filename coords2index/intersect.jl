"""
    intersect(x_sat,los,semi_major_axis=6378137.,flattening=1/298.257223563)

    # Arguments
    - `x_sat::Array{float}(3)`: [ X,Y,Z] of the satellite.
    - `los::Float`: Normalised Line of sight
    # Output
    - `x_0::Array{float}(3)`: intersection between line and elisiod.
"""

function intersect(x_sat,los,semi_major_axis=6378137.,flattening=1/298.257223563)

    semi_minor_axis = semi_major_axis*(1 - flattening)
    epsilon = (semi_major_axis/semi_minor_axis)^2  - 1 # second eccentricity squared
    ecc_squared = flattening*(2-flattening)

    F    = (x_sat'*los + epsilon*x_sat[3]*los[3]) / (1 + epsilon*los[3]^2)
    G    = (x_sat'*x_sat - semi_major_axis^2 + epsilon*x_sat[3]^2) / (1 + epsilon*los[3]^2)
    R    = -F - sqrt(F^2 - G)

    return x_sat + R.* los;
end
