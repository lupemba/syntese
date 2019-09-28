import Polynomials
using Statistics
using Dates


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

##############


function polyval_sv(poly_list,t)
    y_fit = [Polynomials.polyval(poly_elem,t) for poly_elem in poly_list]
    return y_fit
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



###############
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

