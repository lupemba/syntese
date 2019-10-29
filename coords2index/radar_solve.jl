
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
