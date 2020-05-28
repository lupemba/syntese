module UnitTest

include("Misc.jl")
include("Geometry.jl")
include("Load.jl")
import .Misc
import .Load
import .Geometry

import EzXML
import XMLDict
using Dates
using Plots
using JSON




# should be included in meta loader
function meta_start_datetime(meta_path)
    # Define gamma meta
    doc = EzXML.readxml(meta_path)
    meta_dict = XMLDict.xml_dict(doc)
    start_time = meta_dict["product"]["adsHeader"]["startTime"]
    stop_time = meta_dict["product"]["adsHeader"]["stopTime"]
    return start_time, stop_time
end


# !!! FOR SOME FUCKED UP REASON THIS DOESN'T WORK INSIDE MODULE !!!
# maybe move to .Load module
# function gamma_dem_json(json_path)
#     srdem_subset = Dict()
#     open(json_path, "r") do f
#         global srdem_subset
#         dicttxt = read(f)  # file information to string
#         srdem_subset = JSON.parse(String(dicttxt))  # parse and transform data
#     end
# 
#     gamma_dem = Array{Float64}(undef, size(srdem_subset)..., size(srdem_subset[1])...)
#     for i in range(1, size(srdem_subset)...)
#         gamma_dem[i, :] = transpose(srdem_subset[i])
#     end
#     return gamma_dem
# end


# maybe move to misc module
function seconds_since_midnight(date_time)
    if typeof(date_time) == DateTime
        date_time_string = Dates.format(date_time, "yyyy-mm-ddTHH:MM:SS")
    elseif typeof(date_time) == String
        date_time_string = date_time
    end
    # Compute time from DateTime input.
    split_time_string = split(split(date_time_string, "T")[end], ":")
    hour_digits, minute_digits = parse.(Int, split_time_string[1:2])
    second_decimals_string = split(split_time_string[3], ".")
    if length(second_decimals_string) == 1
        sec_digits  = parse.(Int, second_decimals_string[1])
        time_now = Time(hour_digits, minute_digits, sec_digits)
    elseif length(second_decimals_string) == 2
        sec_digits, milisec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3]])
        time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits)
    elseif length(second_decimals_string) == 3
        sec_digits, milisec_digits, microsec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3], second_decimals_string[2][4:end]])
        time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits, microsec_digits)
    elseif length(second_decimals_string) == 4
        sec_digits, milisec_digits, microsec_digits, nanosec_digits = parse.(Int, [second_decimals_string[1], second_decimals_string[2][1:3], second_decimals_string[2][4:6], second_decimals_stringsecond_decimals_string[2][6:end]])
        time_now = Time(hour_digits, minute_digits, sec_digits, milisec_digits, microsec_digits, nanosec_digits)
    else
        println("Only down to nanoseconds supported")
    end

    # Convert midnight of date to Time object
    time_midnight = Time(DateTime(split(date_time_string, "T")[1]))

    # get integer number of nanoseconds between input time and midnight
    time_since_midnight = time_now - time_midnight
    time_since_midnight = time_since_midnight / Nanosecond(1) * 10^(-9)
    return time_since_midnight
end

# hardcoded state vectors from gamma .par parameter file.
function gamma_state_vectors()
    # The state-vectors in both master and slave in gamma starts at same time and has intervals of 10 sec.
    gamma_state_vector_time = collect(20342:10:((14*10)+20342));
    
    state_vector_position_1 = [3517345.9262,    1192277.2110,    6013022.4307] 
    state_vector_velocity_1 = [  6568.41746,      -207.88353,     -3792.11336] 
    state_vector_position_2 = [3582830.3513,    1190083.8479,    5974762.8150] 
    state_vector_velocity_2 = [  6528.33957,      -230.77535,     -3859.73816] 
    state_vector_position_3 = [3647910.1657,    1187661.9861,    5935829.1124] 
    state_vector_velocity_3 = [  6487.49607,      -253.58267,     -3926.92946] 
    state_vector_position_4 = [3712577.7365,    1185012.4887,    5896225.6959] 
    state_vector_velocity_4 = [  6445.89162,      -276.30184,     -3993.67967] 
    state_vector_position_5 = [3776825.4776,    1182136.2556,    5855957.0142] 
    state_vector_velocity_5 = [  6403.53094,      -298.92921,     -4059.98125] 
    state_vector_position_6 = [3840645.8508,    1179034.2230,    5815027.5910] 
    state_vector_velocity_6 = [  6360.41887,      -321.46114,     -4125.82670] 
    state_vector_position_7 = [3904031.3668,    1175707.3634,    5773442.0251] 
    state_vector_velocity_7 = [  6316.56033,      -343.89401,     -4191.20858] 
    state_vector_position_8 = [3966974.5859,    1172156.6854,    5731204.9893] 
    state_vector_velocity_8 = [  6271.96032,      -366.22423,     -4256.11948] 
    state_vector_position_9 = [4029468.1188,    1168383.2333,    5688321.2299] 
    state_vector_velocity_9 = [  6226.62394,      -388.44822,     -4320.55207] 
    state_vector_position_10 = [4091504.6276,    1164388.0873,    5644795.5666]
    state_vector_velocity_10 = [  6180.55637,      -410.56242,     -4384.49904]
    state_vector_position_11 = [4153076.8268,    1160172.3631,    5600632.8917]
    state_vector_velocity_11 = [  6133.76289,      -432.56328,     -4447.95315]
    state_vector_position_12 = [4214177.4839,    1155737.2116,    5555838.1700]
    state_vector_velocity_12 = [  6086.24886,      -454.44728,     -4510.90721]
    state_vector_position_13 = [4274799.4204,    1151083.8190,    5510416.4376]
    state_vector_velocity_13 = [  6038.01972,      -476.21093,     -4573.35409]
    state_vector_position_14 = [4334935.5130,    1146213.4062,    5464372.8021]
    state_vector_velocity_14 = [  5989.08103,      -497.85074,     -4635.28668]
    state_vector_position_15 = [4394578.6942,    1141127.2288,    5417712.4415]
    state_vector_velocity_15 = [  5939.43838,      -519.36326,     -4696.69797];

    m_state_vector_position = vcat(transpose(state_vector_position_1), transpose(state_vector_position_2),
                                    transpose(state_vector_position_3), transpose(state_vector_position_4),
                                    transpose(state_vector_position_5), transpose(state_vector_position_6),
                                    transpose(state_vector_position_7), transpose(state_vector_position_8),
                                    transpose(state_vector_position_9), transpose(state_vector_position_10),
                                    transpose(state_vector_position_11), transpose(state_vector_position_12),
                                    transpose(state_vector_position_13), transpose(state_vector_position_14), 
                                    transpose(state_vector_position_15));

    m_state_vector_velocity = vcat(transpose(state_vector_velocity_1), transpose(state_vector_velocity_2),
                                    transpose(state_vector_velocity_3), transpose(state_vector_velocity_4),
                                    transpose(state_vector_velocity_5), transpose(state_vector_velocity_6),
                                    transpose(state_vector_velocity_7), transpose(state_vector_velocity_8),
                                    transpose(state_vector_velocity_9), transpose(state_vector_velocity_10),
                                    transpose(state_vector_velocity_11), transpose(state_vector_velocity_12),
                                    transpose(state_vector_velocity_13), transpose(state_vector_velocity_14), 
                                    transpose(state_vector_velocity_15));

    m_state_vector_gamma = hcat(m_state_vector_position, m_state_vector_velocity);
        
    state_vector_position_1 = [3513295.9759, 1192470.7396, 6015347.3541]
    state_vector_velocity_1 = [  6570.86980,   -206.52190,  -3787.92515]
    state_vector_position_2 = [3578805.1612, 1190290.9638, 5977129.4895]
    state_vector_velocity_2 = [  6530.83920,   -229.41954,  -3855.57626]
    state_vector_position_3 = [3643910.2072, 1187882.6299, 5938237.2725]
    state_vector_velocity_3 = [  6490.04272,   -252.23291,  -3922.79433]
    state_vector_position_4 = [3708603.4783, 1185246.5988, 5898675.0715]
    state_vector_velocity_4 = [  6448.48498,   -274.95835,  -3989.57179]
    state_vector_position_5 = [3772877.3854, 1182383.7682, 5858447.3304]
    state_vector_velocity_5 = [  6406.17073,   -297.59220,  -4055.90108]
    state_vector_position_6 = [3836724.3875, 1179295.0722, 5817558.5686]
    state_vector_velocity_6 = [  6363.10479,   -320.13083,  -4121.77470]
    state_vector_position_7 = [3900136.9919, 1175981.4811, 5776013.3800]
    state_vector_velocity_7 = [  6319.29206,   -342.57063,  -4187.18520]
    state_vector_position_8 = [3963107.7560, 1172444.0012, 5733816.4329]
    state_vector_velocity_8 = [  6274.73756,   -364.90798,  -4252.12519]
    state_vector_position_9 = [4025629.2875, 1168683.6749, 5690972.4692]
    state_vector_velocity_9 = [  6229.44636,   -387.13932,  -4316.58731]
    state_vector_position_10 = [ 4087694.2451, 1164701.5803, 5647486.3040]
    state_vector_velocity_10 = [   6183.42366,   -409.26107,  -4380.56426]
    state_vector_position_11 = [ 4149295.3402, 1160498.8307, 5603362.8252]
    state_vector_velocity_11 = [   6136.67472,   -431.26970,  -4444.04880]
    state_vector_position_12 = [ 4210425.3368, 1156076.5751, 5558606.9930]
    state_vector_velocity_12 = [   6089.20489,   -453.16169,  -4507.03373]
    state_vector_position_13 = [ 4271077.0532, 1151435.9975, 5513223.8393]
    state_vector_velocity_13 = [   6041.01962,   -474.93352,  -4569.51192]
    state_vector_position_14 = [ 4331243.3626, 1146578.3168, 5467218.4672]
    state_vector_velocity_14 = [   5992.12444,   -496.58173,  -4631.47625]
    state_vector_position_15 = [ 4390917.1941, 1141504.7866, 5420596.0504]
    state_vector_velocity_15 = [   5942.52497,   -518.10285,  -4692.91972]

    s_state_vector_position = vcat(transpose(state_vector_position_1), transpose(state_vector_position_2),
                                    transpose(state_vector_position_3), transpose(state_vector_position_4),
                                    transpose(state_vector_position_5), transpose(state_vector_position_6),
                                    transpose(state_vector_position_7), transpose(state_vector_position_8),
                                    transpose(state_vector_position_9), transpose(state_vector_position_10),
                                    transpose(state_vector_position_11), transpose(state_vector_position_12),
                                    transpose(state_vector_position_13), transpose(state_vector_position_14), 
                                    transpose(state_vector_position_15));

    s_state_vector_velocity = vcat(transpose(state_vector_velocity_1), transpose(state_vector_velocity_2),
                                    transpose(state_vector_velocity_3), transpose(state_vector_velocity_4),
                                    transpose(state_vector_velocity_5), transpose(state_vector_velocity_6),
                                    transpose(state_vector_velocity_7), transpose(state_vector_velocity_8),
                                    transpose(state_vector_velocity_9), transpose(state_vector_velocity_10),
                                    transpose(state_vector_velocity_11), transpose(state_vector_velocity_12),
                                    transpose(state_vector_velocity_13), transpose(state_vector_velocity_14), 
                                    transpose(state_vector_velocity_15));

    s_state_vector_gamma = hcat(s_state_vector_position, s_state_vector_velocity);
    
    return m_state_vector_gamma, s_state_vector_gamma, gamma_state_vector_time
end
                    
                    
end  # end module
