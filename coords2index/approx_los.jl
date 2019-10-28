using LinearAlgebra



"""
    approx_los(osv_fit,s1_annn)ters)

    # Arguments
    - `osv_fit::Array{float}(6)`: [ X,Y,Z,V_x,V_y,V_z] of the satellite.
    - `s1_annn::Dict`: Dict with relevant meta info. See load_s1slc_ann(path)
    # Output
    - `los::Array{float}(3)`: Line of sight to mid swath in elipsiodal coordinates
"""
function approx_los(osv_fit,s1_annn)
    theta_0 = s1_annn["incidence_angle_mid"]*pi/180

    x_sat = osv_fit[1:3]
    v_sat = osv_fit[4:6]

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
