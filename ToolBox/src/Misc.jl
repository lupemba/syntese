module Misc


"""
    print2maps_co(lat,lon,name="Corner",color="#FF0000")

    Print lat,lon in format there can be copy pasted to maps.co.
    maps.co is an easy way to plot points on a interactive map.
"""
function print2maps_co(lat,lon,name="Corner",color="#FF0000")
    for i=1:length(lat)
        println(lat[i],",",lon[i],",",name,",",color)
    end
end


end
