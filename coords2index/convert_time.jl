import Dates

"""
    str2date(time)

Convert a string of format "2017-03-15T05:39:50.703105" to a DateTime with an
accuarcy of minutes
"""
function str2date(time)
    y = parse(Int,time[1:4])
    m = parse(Int,time[6:7])
    d = parse(Int,time[9:10])
    h = parse(Int,time[12:13])
    mi = parse(Int,time[15:16])

    return Dates.DateTime(y, m, d, h, mi)
end

"""
    str_date2float(time,t_0)

    # Arguments
    - `time::String`: date of format of format "2017-03-15T05:39:50.703105"
    - `t_0::DateTime`: Reference time

    # Output
    - `dt::Float64`: time relative to t_0 in seconds.
"""
function str_date2float(time,t_0)
    # find time diff with date time down to min
    t_i = str2date(time)
    dt_i = Dates.Second(t_i-t_0).value

    # convert seconds seperately to get float pression
    s = parse(Float64,time[18:end])

    return dt_i + s
end
