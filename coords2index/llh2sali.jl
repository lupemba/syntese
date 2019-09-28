include("zero_doppler_bisect.jl")
include("ellipsoid2xyz.jl");
include("calc_sat_trajectory.jl");
include("load_pod.jl");


"""
    llh2sali(llh, osv, t_sv, sar_parameters)

    # Arguments
    - `llh::Array{float}(Nx3)`: array N points of latitude(deg),longitude(deg),heigt(m)
    - `osv::Array{float}(Lx6)`: Array with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `t_sv::Array{float}(L)`: time of each orbit state relative to t_0 in seconds.
    - `sar_parameters::Dict`: Dict with relevant meta info. See load_s1slc_ann(path)
    # Output
    - `sali::Array{float}(Nx2)`: - Array of points [line,sample]
"""

function llh2sali(llh, osv, t_sv, sar_parameters)

    c = 299792458 # speed of light

    t_0 = sar_parameters["t_0"]
    t_start = sar_parameters["t_start"]
    t_stop = sar_parameters["t_stop"]

    inv_range_pixel_spacing =  1/sar_parameters["range_pixel_spacing"]
    azimuth_frequency =  sar_parameters["azimuth_frequency"]
    r_near =  sar_parameters["slant_range_time"]  *c/2
    deg2rad = pi/180

    osv_poly, osv_mean, osv_std = calc_sat_trajectory(osv, t_sv, t_start, t_stop)
    
    
    sali = Array{Float64}(undef,size(llh)[1],2)

    for i in 1:size(llh)[1]
        point_xyz = ellipsoid2xyz(llh[i,1]* deg2rad, llh[i,2]* deg2rad ,llh[i,3])
        time,range = zero_doppler_bisect(point_xyz, t_start, t_stop,
                                        osv_poly, osv_mean, osv_std)
        sali[i,1] = 1 + (time - t_start) * azimuth_frequency
        sali[i,2] = 1 + (range - r_near) * inv_range_pixel_spacing
    end

    return sali

end
