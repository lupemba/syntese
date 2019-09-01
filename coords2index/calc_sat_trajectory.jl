import Polynomials
using Statistics
using Dates


"""
    polyfit_sv(y_list, t; poly_degree=4, do_scale=1)

Fit polynomials to each array in y_list using t as the argument

# Arguments
- `y_list`: Array of Arrays with dependent variables
- `t`: Array of independent variable
# returns
- `poly_list`: Array of fitted polynomials
- `y_mean`: Array of mean (only if, do_scale=1)
- `y_std`: Array of standard deviation (only if, do_scale=1)
"""
function polyfit_sv(y_list, t; poly_degree=4, do_scale=1)
    y_mean = mean.(y_list)
    y_std = std.(y_list)

    if do_scale==1
        y_norm = [(y_list[i].-y_mean[i])./y_std[i] for i in range(1,stop=length(y_list))]
        poly_list = [Polynomials.polyfit(t, column,poly_degree) for column in y_norm]
        return poly_list, y_mean, y_std
    elseif do_scale==0
        poly_list = [Polynomials.polyfit(t, column, poly_degree) for column in y_list]
        return poly_list
    else
        println("Error: do_scal can be 1 or 0")
        return
    end
end

"""
    polyfit_sv(y_list, t::Array{DateTime,1}; step=Second, poly_degree=4, do_scale=1)

Convert t to "step" from first time, eg. seconds from t[0] and calls
polyfit_sv(y_list, t, poly_degree, do_scale)

"""
function polyfit_sv(y_list, t::Array{DateTime,1}; step=Second, poly_degree=4, do_scale=1)

    # Convert t to "step" from first time, eg. seconds from t[0]
    dt = t.-t[1]
    dt = [convert(step,elem).value for elem in dt]

    # call standard function
    return polyfit_sv(y_list, dt, poly_degree=poly_degree, do_scale=do_scale)
end





function polyval_sv(poly_list,t)
    y_fit = [Polynomials.polyval(poly_elem,t) for poly_elem in poly_list]
    return y_fit
end

function polyval_sv(poly_list,t,y_mean,y_std)
    y_norm_fit = polyval_sv(poly_list,t)

    # re-normalize the fitted data
    y_fit = [y_norm_fit[i].*y_std[i].+y_mean[i] for i in range(1,stop=length(y_norm_fit))];
    return y_fit
end

"""
    function polyval_sv(poly_list,t::Array{DateTime,1},y_mean,y_std,t0::DateTime,step=Second)

evaluate normalized polynomials for time, t

# Arguments
- `poly_list`: List of normalized fitted polynomials for [X,Y,Z,VX,VY,VZ]
- `y_mean`: mean of the data
- `t::Array{DateTime,1}: times to evaluate the polynomials
- `y_std`: standard deviation of the data
- `t0`: zero time point for polynomial fit

# Returns
- `y_fit`: Array of Array with each polynomial evaluate for all the times.

"""
function polyval_sv(poly_list,t::Array{DateTime,1},y_mean,y_std,t0::DateTime,step=Second)
    # Convert t to "step" from first time, eg. seconds from t0
    dt = t.-t0
    dt = [convert(step,elem).value for elem in dt]
    # call standard function
    return polyval_sv(poly_list,dt,y_mean,y_std)
end



"""
    calc_sat_trajectory(osv,t_sv, t_start, t_stop; poly_degree=4, max_time_margin=Second(240))

Fit normalized polynomials to orbit state vectors for the time period of t_start to t_stop.

# Arguments
- `osv`: Array with [X,Y,Z,VX,VY,VZ] data
- `t_sv::Array{DateTime,N}`: DateTimes associated with osv
- `t_start::Datetime`: start time
- `t_stop::Datetime`: stop time

# Returns
- `osv_poly`: List of normalized fitted polynomials for [X,Y,Z,VX,VY,VZ]
- `osv_mean`: mean of osv
- `osv_std`: standard deviation of osv
- `t0`: zero time point for polynomial fit

"""
function calc_sat_trajectory(osv,t_sv, t_start, t_stop; poly_degree=4, max_time_margin=Second(240))
    dt = t_stop - t_start;
    t_mid = t_start + dt / 2
    t_step = t_sv[2]-t_sv[1]
    max_time_margin = max(convert(typeof(t_step),max_time_margin), (poly_degree + 2) / 2 * t_step)


    nearby_state_vector_idx = abs.(t_sv .- t_mid) .<= dt/2 + max_time_margin
    if sum(nearby_state_vector_idx) < poly_degree + 1
        println("Insufficient number of state vectors exiting.")
        return
    end

    nearby_state_vector_t_sv = t_sv[nearby_state_vector_idx];
    nearby_state_vector_osv = [elem[nearby_state_vector_idx] for elem in osv];
    t0 = nearby_state_vector_t_sv[1]

    osv_poly, osv_mean, osv_std =polyfit_sv(nearby_state_vector_osv,nearby_state_vector_t_sv,poly_degree=poly_degree)

    return osv_poly, osv_mean, osv_std, t0

end
