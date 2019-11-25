module Geometry

import Polynomials
using LinearAlgebra
using Statistics
using Dates

export salih2ll, llh2sali

"""
    salih2llh(sali, height , osv, t_sv, sar_parameters)

    Convert from SAR coordinates to latitude and longitude

    # Arguments
    - `sali::Array{float}(Nx2)`: - Array of points line,sample
    - `height::Array{float}(N)`: - Array with heights (meters) of points
    - `osv::Array{float}(Lx6)`: Array with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `t_sv::Array{float}(L)`: time of each orbit state relative to t_0 in seconds.
    - `sar_parameters::Dict`: Dict with relevant meta info. See load_s1slc_ann(path)
    # Output
    - `llh::Array{float}(Nx3)`: array N points of latitude(deg),longitude(deg),heigt(m)
"""
function salih2ll(sali, height, osv, t_sv, sar_parameters)
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

    ll = Array{Float64}(undef,size(sali)[1],2)


    for i in 1:size(sali)[1]
        time =  t_start + (sali[i,1]-1)*inv_azimuth_frequency
        range_x = r_near + (sali[i,2] - 1)*range_pixel_spacing

        osv_0 = polyval_sv(osv_poly,time,osv_mean, osv_std)
        x_sat = osv_0[1:3]
        v_sat = osv_0[4:6]

        x_0 = intersect(x_sat, approx_los(x_sat,v_sat,theta_0))
        x = radar_solve(range_x,height[i],x_0,x_sat,v_sat)
        point =  xyz2ellh(x...)
        ll[i,1] = point[1]*rad2deg
        ll[i,2] = point[2]*rad2deg

    end

    return ll
end


"""
    llh2sali(llh, osv, t_sv, sar_parameters)

    Convert from latitude longitude and height to SAR coordintes

    # Arguments
    - `llh::Array{float}(Nx3)`: array N points of latitude(deg),longitude(deg)
    - `height::Array{float}(N)`: - Array with heights (meters) of points
    - `osv::Array{float}(Lx6)`: Array with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `t_sv::Array{float}(L)`: time of each orbit state relative to t_0 in seconds.
    - `sar_parameters::Dict`: Dict with relevant meta info. See load_s1slc_ann(path)
    # Output
    - `sali::Array{float}(Nx2)`: - Array of points [line,sample]
"""
function llh2sali(ll, height, osv, t_sv, sar_parameters)

    c = 299792458 # speed of light

    t_0 = sar_parameters["t_0"]
    t_start = sar_parameters["t_start"]
    t_stop = sar_parameters["t_stop"]

    inv_range_pixel_spacing =  1/sar_parameters["range_pixel_spacing"]
    azimuth_frequency =  sar_parameters["azimuth_frequency"]
    r_near =  sar_parameters["slant_range_time"]  *c/2
    deg2rad = pi/180

    osv_poly, osv_mean, osv_std = calc_sat_trajectory(osv, t_sv, t_start, t_stop)


    sali = Array{Float64}(undef,size(ll)[1],2)

    for i in 1:size(ll)[1]
        point_xyz = ellipsoid2xyz(ll[i,1]* deg2rad, ll[i,2]* deg2rad ,height[i])
        time,range = zero_doppler_bisect(point_xyz, t_start, t_stop,
                                        osv_poly, osv_mean, osv_std)
        sali[i,1] = 1 + (time - t_start) * azimuth_frequency
        sali[i,2] = 1 + (range - r_near) * inv_range_pixel_spacing
    end

    return sali
end


"""
    calc_sat_trajectory(osv,t_sv, t_start, t_stop; poly_degree=4,
                        max_time_margin=240)
    Fit normalized polynomials to orbit state vectors (osv) for the time period of
    t_start to t_stop.

    # Arguments
    - `osv::Array(Nx6)`: Array with [X,Y,Z,VX,VY,VZ] data
    - `t_sv::Float`: time in seconds
    - `t_start::Float`: start time
    - `t_stop::Float`: stop time

    # Returns
    - `osv_poly::Array(6)`: Array of normalized fitted polynomials for [X,Y,Z,VX,VY,VZ]
    - `osv_mean::Array(6)`: mean of osv
    - `osv_std::Array(6)`: standard deviation of osv
"""
function calc_sat_trajectory(osv,t_sv, t_start, t_stop; poly_degree=4, max_time_margin=240.)
    dt = t_stop - t_start;
    t_mid = t_start + dt / 2
    t_step = t_sv[2]-t_sv[1]
    max_time_margin = max(max_time_margin, (poly_degree + 2) / 2 * t_step)


    nearby_state_vector_idx = abs.(t_sv .- t_mid) .<= dt/2 + max_time_margin

    # Check if enough point are selected
    if sum(nearby_state_vector_idx) < poly_degree + 1
        println("Insufficient number of state vectors exiting.")
        return
    end

    # slice data to close times
    nearby_state_vector_t_sv = t_sv[nearby_state_vector_idx];
    nearby_state_vector_osv = osv[nearby_state_vector_idx,:];

    osv_poly, osv_mean, osv_std =polyfit_sv(nearby_state_vector_osv,nearby_state_vector_t_sv,poly_degree=poly_degree)

    return osv_poly, osv_mean, osv_std
end


"""
    polyfit_sv(y, t; poly_degree=4, do_scale=1)

    Fit polynomials to each array in y_list using t as the argument

    # Arguments
    - `y:: Array(NxM)`: Array with N observations of M dependd varibles
    - `t::Array(N)`: Array of N times

    # Returns
    - `poly_list::Array(M)`: M normalized polynomials
    - `y_mean::Array(M)`: Array of mean (only if, do_scale=1)
    - `y_std::Array(M)`: Array of standard deviation (only if, do_scale=1)
"""
function polyfit_sv(y, t; poly_degree=4, do_scale=1)
    y_mean = mean(y,dims=1)
    y_std = std(y,dims=1)

    if do_scale==1
        y_norm = (y.-y_mean)./y_std
        poly_list = [Polynomials.polyfit(t, y_norm[:,i],poly_degree) for i in 1:size(y)[2]]
        return poly_list, dropdims(y_mean, dims=1), dropdims(y_std, dims=1)
    elseif do_scale==0
        poly_list = [Polynomials.polyfit(t, y[:,i],poly_degree) for i in 1:size(y)[2]]
        return poly_list
    else
        println("Error: do_scal can be 1 or 0")
        return
    end
end


"""
    function polyval_sv(poly_list,t,y_mean,y_std)

    evaluate normalized polynomials for time, t

    # Arguments
    - `poly_list::Array(M)`: Array of normalized fitted polynomials for [X,Y,Z,VX,VY,VZ]
    - `t`: Number of seconds after t_0
    - `y_std::Array(M)`: standard deviation of the data
    - `y_mean::Array(M)`: mean of the data

    # Returns
    - `y_fit`: Array of Array with each polynomial evaluate for all  t_0 + Seconds(t).
"""
function polyval_sv(poly_list,t,y_mean,y_std)
    y_norm_fit = polyval_sv(poly_list,t)

    # re-normalize the fitted data
    y_fit = y_norm_fit .* y_std .+ y_mean
    return y_fit
end
function polyval_sv(poly_list,t)
    y_fit = [Polynomials.polyval(poly_elem,t) for poly_elem in poly_list]
    return y_fit
end

"""
    radar_solve(range_x,height,x_init,x_sat,v_sat)

    Find the point that is range_x away from the satelite, orthogonal on the flight directions
    and "height" above the elipsiod using Newton_rhapsody method.

    # Arguments
    - `range_x::Float`: distance from satellite to point in meters
    - `height::Float`: height of the point over the referrence elipsoid in meters
    - `x_init::Array{float}(3)`: initial guess of the point .
    - `x_sat::Array{float}(3)`: [ X,Y,Z] of the satellite.
    - `v_sat::Array{float}(3)`: [ V_x,V_y,V_z] of the satellite.
    # Output
    - `x_i::Array{float}(3)`: [ X,Y,Z] of the point.
"""
function radar_solve(range_x,height,x_init,x_sat,v_sat,
        scale_fac = 1e-03,MAX_ITER = 150,eps = 1e-6,
         semi_major_axis=6378137.,flattening=1/298.257223563)


        semi_minor_axis = semi_major_axis*(1 - flattening)

        # inits
        x_i_old = [0.;0.;0.];
        x_i = x_init;
        iter  = 1;
        los = [0.;0.;0.];

        # scale
        v_sat = v_sat.* scale_fac;
        x_i = x_i .*scale_fac;
        x_sat = x_sat.*scale_fac;
        range_x = range_x .*scale_fac;
        a_plus_h = (semi_major_axis + height) .*scale_fac;
        b_plus_h = (semi_minor_axis + height) .*scale_fac;




        while ((x_i - x_i_old)'*(x_i - x_i_old) > eps^2) & (iter < MAX_ITER)

            # Design matrix evaluated at previous solution
            los = x_i - x_sat
            fx_i = [v_sat'* los,
                los' *los - range_x^2,
                ((x_i[1]^2 + x_i[2]^2) / a_plus_h^2 +(x_i[3] / b_plus_h)^2 - 1)];
            # Matrix of partial derivatives
            dfx_i = vcat(v_sat',
                    2*los',
                    [2*x_i[1]/a_plus_h^2, 2*x_i[2]/a_plus_h^2,2*x_i[3]/b_plus_h^2]');

            # Solve linear system
            dx = dfx_i\(-fx_i)
            # Update
            x_i_old = x_i;
            x_i += dx;
            iter += 1;


        end

        if iter == MAX_ITER
                println("Warning Covergens not reached")
        end

        return x_i./scale_fac
end

"""
    zero_doppler_bisect(point_xyz, azimuth_start_time, azimuth_stop_time,
                        osv_poly, osv_mean, osv_std, t0)

    Compute the (slant-range, zero Doppler time)-coordinates for a given point on
    ground.

    # Arguments
    - `point_xyz`: 3-element Array{Float64,1} target point
    - `azimuth_start_time`: Float64, the start time of the subswath in seconds
    - `azimuth_stop_time`: Float64, the end time of the subswath in seconds
    - `osv_poly`: 6-element Array{Polynomials.Poly{Float64},1} see
    calc_sat_trajectory() documentation
    - `osv_mean`: 6-element Array{Float64,1} see calc_sat_trajectory() documentation
    - `osv_std`: 6-element Array{Float64,1} see calc_sat_trajectory() documentation

    # Examples:
    ```jldoctest
    julia> include("load_pod.jl"); include("calc_sat_trajectory.jl");
    julia> f = open("POD_path.txt"); path = readlines(f);
    julia> osv,t_sv = load_pod(path[1]);
    julia> azimuth_start_time = DateTime(2017, 3, 15, 5, 40);
    julia> azimuth_stop_time = DateTime(2017, 3, 15, 5, 42);
    julia> osv_poly, osv_mean, osv_std = calc_sat_trajectory(osv, t_sv, azimuth_start_time, azimuth_stop_time);
    julia> point_xyz = [4.558689828520125e6, 1.3632875482646455e6, 5.716729814351776e6];
    julia> rrange, ttime = zero_doppler_bisect(point_xyz, azimuth_start_time, azimuth_stop_time, osv_poly, osv_mean, osv_std)
    2-element Array{Float64,1}:
    400999.9999999987
    310.3399994522333
    ```

    # Returns
    - `trial_time`: the time after t0 where the satellite is at zero doppler
    wrt to target
    - `sqrt(line_of_sight' * line_of_sight)`: The slant range i.e. the norm of the
    line_of_sight vector.
"""
function zero_doppler_bisect(point_xyz, azimuth_start_time, azimuth_stop_time,
                             osv_poly, osv_mean, osv_std, debug = 0)
    #=
    TO-DO
    Change while loop to more advanced optimization
    =#

    #= Define interval for zero Doppler time search, with
    an interval larger than the data extend for sensitivity
    at the edges =#
    data_duration = azimuth_stop_time - azimuth_start_time
    search_interval_start = (azimuth_start_time - data_duration)
    search_interval_end = (azimuth_stop_time + data_duration)

    # Initial bisection parameters, start iteration with entire data interval
    search_interval_duration = search_interval_end - search_interval_start
    small_time_interval = 1e-6
    max_iter = 100
    iter = 0
    line_of_sight = 0
    trial_time = 0
    osv = [0.,0.,0.,0.,0.,0.]

    while (search_interval_duration > small_time_interval) && (iter < max_iter)

        search_interval_duration = (search_interval_end - search_interval_start) / 2
        trial_time = search_interval_start + search_interval_duration

        osv = polyval_sv(osv_poly, trial_time, osv_mean, osv_std)
        trial_sat_position = osv[1:3]
        line_of_sight = point_xyz .- trial_sat_position

        # Compute the squint angle
        trial_sat_velocity = osv[4:6]
        sin_squint_angle = line_of_sight' * trial_sat_velocity

        if (sin_squint_angle < 0)
            search_interval_end = trial_time
        else
            search_interval_start = trial_time
        end
        iter += 1;
    end

    if iter >= max_iter
        print("Error, max iteration reached")
    end

    if debug == 0
        return [trial_time,sqrt(line_of_sight' * line_of_sight)]
    else
        return [trial_time,sqrt(line_of_sight' * line_of_sight), osv,line_of_sight, iter]
    end
end

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

"""
    ellipsoid2xyz(lat::Array{Number,N}, lon::Array{Number,N},
                  height::Array{Number,N}; semi_major_axis=6378137.,
                  flattening=1/298.257223563)

    Convert ellipsoid coordinates to cartisian coordinates, based on ellipsoidal
    latitudes (rad), longitudes (rad) and heights (m), deafult values are wrt to GRS80.

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
function ellipsoid2xyz(lat, lon, height; semi_major_axis=6378137.,
                       flattening=1/298.257223563)
    #=
    TO-DO
    specify arugment types for lat, lon, height.
    =#
    e2 = flattening * (2 - flattening)
    v = semi_major_axis./sqrt.(1 .- e2*sin.(lat).*sin.(lat))
    x=(v+height).*cos.(lat).*cos.(lon)
    y=(v+height).*cos.(lat).*sin.(lon)
    z=(v.*(1-e2)+height).*sin.(lat)
    return x, y, z
end

"""
    xyz2ellh(X,Y,Z,semi_major_axis=6378137.,flattening=1/298.257223563)

    Go from xyz elipsoidal coordinates to latitude, longitude, height
    (see B.R. Bowring, "The accuracy of geodetic latitude and height equations",
    Survey Review, v28 #218, October 1985 pp.202-206).
"""
function xyz2ellh(X,Y,Z,semi_major_axis=6378137.,flattening=1/298.257223563)
    e2 = flattening*(2-flattening)


    elat=1.e-12
    eht=1.e-5
    p=sqrt(X^2+Y^2)
    lat=atan(Z,p./(1-e2))
    h=0
    dh=1
    dlat=1


    while (dlat>elat) | (dh>eht)
      lat0=lat
      h0=h
      v=semi_major_axis/sqrt(1-e2*sin(lat)*sin(lat))
      h=p*cos(lat)+Z*sin(lat)-(semi_major_axis^2)/v  # Bowring formula
      lat=atan(Z, p*(1-e2*v/(v+h)))
      dlat=abs(lat-lat0)
      dh=abs(h-h0)
    end
    lon=atan(Y,X)
    return lat, lon, h
end

end
