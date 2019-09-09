"""
    zero_doppler_bisect()



# Examples:
```jldoctest
```
"""
using Dates

function zero_doppler_bisect()
    #=
    TO-DO
    Change while loop to more advanced optimization
    =#
    ## SIMON
    
    data_duration = azimuth_stop_time - azimuth_start_time
    search_interval_start = (azimuth_start_time - data_duration) - t0
    search_interval_end = (azimuth_stop_time + data_duration) - t0
    
    
    data_duration = Millisecond(data_duration).value/1000
    search_interval_start = Millisecond(search_interval_start).value/1000
    search_interval_end = Millisecond(search_interval_end).value/1000
    # Initial bisection parameters, start iteration with entire data interval
    search_interval_duration = search_interval_end - search_interval_start
    
    
    
    #####
    #= Define interval for zero Doppler time search, with
    an interval larger than the data extend for sensitivity
    at the edges =#
    data_duration = azimuth_stop_time - azimuth_start_time
    search_interval_start = azimuth_start_time - data_duration
    search_interval_end = azimuth_stop_time + data_duration

    # Initial bisection parameters, start iteration with entire data interval
    search_interval_duration = search_interval_end - search_interval_start
    
    
    small_time_interval = 1e-6
    max_iter = 100
    iter = 0

    while (search_interval_duration > small_time_interval) && (iter < max_iter)
        search_interval_duration = (search_interval_end - search_interval_start) / 2
        trial_time = search_interval_start + search_interval_duration
        
        osv = polyval_sv(osv_poly,trial_time,osv_mean, osv_std)
        
        trial_sat_position = osv[1:3]

        line_of_sight = point_xyz - trial_sat_position

        # Compute the squint angle
        trial_sat_velocity = = osv[4:6]
        sin_squint_angle = line_of_sight * trial_sat_velocity

        if (sin_squint_angle < 0)
            search_interval_end = trial_time
        else
            search_interval_start = trial_time
        end
        iter += 1;
    end

    if iter == max_iter
        print("Error, max iteration reached")
    end

    return [norm(line_of_sight), t0+ Second(trial_time)]
end
