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
