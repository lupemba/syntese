import EzXML
import XMLDict
import Dates

function load_pod(path,t_0)

    # Load data as dict
    doc = EzXML.readxml(path)
    pod_dict = XMLDict.xml_dict(doc)

    # Acces orbit state vectors
    osv_dict = pod_dict["Earth_Explorer_File"]["Data_Block"]["List_of_OSVs"]["OSV"];

    # get vectors
    tags = ["X","Y","Z","VX","VY","VZ"]
    osv = [[parse(Float64,elem[tag][""]) for elem in osv_dict] for tag in tags];

    # get times
    t_sv = Array{Float64}(undef, length(osv_dict))
    for i in range(1,stop=length(osv_dict))
        elem = osv_dict[i]["UTC"]
        y = parse(Int,elem[5:8]);
        m = parse(Int,elem[10:11]);
        d = parse(Int,elem[13:14]);
        h = parse(Int,elem[16:17]);
        mi = parse(Int,elem[19:20]);
        
        # find time diff with date time down to min
        t_i = Dates.DateTime(y, m, d, h, mi)
        dt_i = Dates.Second(t_i-t_0).value
        
        # convert seconds seperately to get float pression
        s = parse(Float64,elem[22:end])
        
        t_sv[i] = dt_i + s
    end

    return osv,t_sv
end
