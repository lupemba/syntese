import EzXML
import XMLDict
import Dates

include("convert_time.jl")

"""
    load_pod(path,t_0)

    Load preciese orbit files

    # Arguments
    - `path::String`: path to precise orbits file
    - `t_0::DateTime`: Reference time

    # Output
    - `osv:: Array{float}(Nx6)`: six arrays with respectvely X,Y,Z,V_x,V_y,V_z observationer.
    - `t_sv::Array{float}(N)`: time of each orbit state relative to t_0 in seconds.
"""


function load_pod(path,t_0)

    # Load data as dict
    doc = EzXML.readxml(path)
    pod_dict = XMLDict.xml_dict(doc)

    # Acces orbit state vectors
    osv_dict = pod_dict["Earth_Explorer_File"]["Data_Block"]["List_of_OSVs"]["OSV"];

    # get vectors
    tags = ["X","Y","Z","VX","VY","VZ"]
    osv = [[parse(Float64,elem[tag][""]) for elem in osv_dict] for tag in tags];
    osv = hcat(osv[1],osv[2],osv[3],osv[4],osv[5],osv[6])
    # get times
    t_sv = [str_date2float(elem["UTC"][5:end],t_0) for elem in osv_dict]

    return osv,t_sv
end

