
include("calc_sat_trajectory.jl");
include("approx_los.jl");
include("intersect.jl")
include("xyz2ellh.jl");
include("radar_solve.jl");


"""
    salih2llh(salih, osv, t_sv, sar_parameters)

    Convert from SAR coordinates to latitude and longitude

    # Arguments
    - `sali::Array{float}(Nx3)`: - Array of points line,sample,height
    - `osv::Array{float}(Lx6)`: Array with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `t_sv::Array{float}(L)`: time of each orbit state relative to t_0 in seconds.
    - `sar_parameters::Dict`: Dict with relevant meta info. See load_s1slc_ann(path)
    # Output
    - `llh::Array{float}(Nx3)`: array N points of latitude(deg),longitude(deg),heigt(m)
"""
function salih2llh(salih, osv, t_sv, sar_parameters)
    c = 299792458 # speed of light
    t_start = sar_parameters["t_start"]
    t_stop = sar_parameters["t_stop"]
    sign_angle  = sar_parameters["right_looking"] ? 1 : -1
    theta_0 = sign_angle*abs(sar_parameters["incidence_angle_mid"]*pi/180)
    range_pixel_spacing =  sar_parameters["range_pixel_spacing"]
    inv_azimuth_frequency =  1/sar_parameters["azimuth_frequency"]
    r_near =  sar_parameters["slant_range_time"]  *c/2

    osv_poly, osv_mean, osv_std = calc_sat_trajectory(osv, t_sv, t_start, t_stop)
    rad2deg = 180/pi

    llh = Array{Float64}(undef,size(salih)[1],3)


    for i in 1:size(salih)[1]
        time =  t_start + (salih[i,1]-1)*inv_azimuth_frequency
        range_x = r_near + (salih[i,2] - 1)*range_pixel_spacing

        osv_0 = polyval_sv(osv_poly,time,osv_mean, osv_std)
        x_sat = osv_0[1:3]
        v_sat = osv_0[4:6]

        x_0 = intersect(x_sat, approx_los(x_sat,v_sat,theta_0))
        x = radar_solve(range_x,salih[i,3],x_0,x_sat,v_sat)
        point =  xyz2ellh(x...)
        llh[i,1] = point[1]*rad2deg
        llh[i,2] = point[2]*rad2deg
        llh[i,3] =point[3]

    end

    return llh

end
