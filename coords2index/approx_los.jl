using LinearAlgebra



"""
    approx_los(osv_fit,s1_annn)

    # Arguments
    - `x_sat::Array{float}(3)`: [ X,Y,Z] of the satellite.
    - `v_sat::Array{float}(3)`: [ V_x,V_y,V_z] of the satellite.
    # Output
    - `los::Array{float}(3)`: Line of sight to mid swath in elipsiodal coordinates
"""
function approx_los(x_sat,v_sat,theta_0)

    #ECEF basis coordinates
    x_hat_ecef = x_sat / sqrt(x_sat'*x_sat) # Towards earth center
    z_hat_ecef = v_sat / sqrt(v_sat'*v_sat) # flight direction
    y_hat_ecef = cross(z_hat_ecef, x_hat_ecef) # Right handed coordinate system

    # Line of sight ECEF basis
    losSat = [-cos(theta_0), sin(theta_0), 0]

    # Basis change matrix from ECEF basis to elipsidal coordinates
    m = hcat(x_hat_ecef, y_hat_ecef,z_hat_ecef)

    # Line of sight in Ellipsoidal coordinates
    los = m*losSat;

    return los
end
