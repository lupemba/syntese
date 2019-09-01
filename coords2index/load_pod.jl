import EzXML
import XMLDict
using Dates

function load_pod(path)

    # Load data as dict
    doc = EzXML.readxml(path)
    pod_dict = XMLDict.xml_dict(doc)

    # Acces orbit state vectors
    osv_dict = pod_dict["Earth_Explorer_File"]["Data_Block"]["List_of_OSVs"]["OSV"];

    # get vectors
    tags = ["X","Y","Z","VX","VY","VZ"]
    osv = [[parse(Float64,elem[tag][""]) for elem in osv_dict] for tag in tags];

    # get times
    t_sv = Array{DateTime}(undef, length(osv_dict))
    for i in range(1,stop=length(osv_dict))
        elem = osv_dict[i]["UTC"]
        y = parse(Int,elem[5:8]);
        m = parse(Int,elem[10:11]);
        d = parse(Int,elem[13:14]);
        h = parse(Int,elem[16:17]);
        mi = parse(Int,elem[19:20]);
        s = parse(Int,elem[22:23]);
        ms = parse(Int,elem[25:27]);
        t_sv[i] = DateTime(y, m, d, h, mi, s, ms)
    end

    return osv,t_sv
end
