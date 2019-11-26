

file = open("ann_path.txt")
path = readlines(file)
sar_parameters = load_s1slc_ann(path[1])

file = open("POD_path.txt")
path = readlines(file)
osv, t_sv = load_pod(path[1],sar_parameters["t_0"]);

sali = llh2sali(llh, osv, t_sv, sar_parameters);
