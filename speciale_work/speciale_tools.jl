include("../ToolBox/ToolBox.jl")
using .ToolBox
using Colors
import PyCall
import StatsBase;
import LsqFit
import Statistics
import FileIO
import Optim
import Distributions
import LinearAlgebra
using Dates

# constant for image visualisations 
min_vv = -19
max_vv = 4

min_vh = -26
max_vh = -3



function loss_rg_fit(VV,VH, seed_mask, rg_mask_change,target_mu,target_sigma, distance_list, d)

    rg_mask_int = ((VV .< target_mu[1]) .& (VH .< target_mu[2])) .| (distance_list .< d)
    rg_mask = rg_mask_change .& rg_mask_int
    flood_mask,steps = region_growing(seed_mask,rg_mask_int .& rg_mask_change);
    loss = loss_2d_gauss(vec(VV[flood_mask]),vec(VH[flood_mask]),target_mu,target_sigma;n_levels=5)
    return loss
end

function loss_2d_gauss(VV,VH,target_mu,target_sigma;n_levels=5, max_q = 1)
    freq = freq_2d_gauss(VV,VH,target_mu,target_sigma;n_levels=n_levels, max_q =max_q)
    ture_freq = max_q/(n_levels*4)
    loss = sum((freq .- ture_freq).^2)
    return loss
end


function freq_2d_gauss(VV,VH,target_mu,target_sigma;n_levels=5,max_q = 1)
    
    N = length(VV)
    distri = Distributions.Chisq(2)
    chi_limits = Distributions.quantile.(distri, collect(range(0,length=n_levels+1,stop=max_q))[2:end])
    
    line1_mask, line2_mask = get_pca_mask(VV,VH,target_mu,target_sigma)
    malhob = get_malhob(VV,VH,target_mu,target_sigma)
    
    frequency = zeros(Float64,(4,n_levels))
    
    frequency[1,:] .= count_chi_limit(malhob[   line1_mask .&   line2_mask],chi_limits)./N
    frequency[2,:] .= count_chi_limit(malhob[ .!line1_mask .&   line2_mask],chi_limits)./N
    frequency[3,:] .= count_chi_limit(malhob[ .!line1_mask .& .!line2_mask],chi_limits)./N
    frequency[4,:] .= count_chi_limit(malhob[   line1_mask .& .!line2_mask],chi_limits)./N
    
    return frequency
end

function get_pca_mask(VV,VH,target_mu,target_sigma)
    pc = LinearAlgebra.eigvecs(target_sigma)
    
    a1 = pc[1,1] /pc[2,1]
    b1 = target_mu[2] - target_mu[1] * a1
    line1_mask = VH .>  (VV.* a1 .+ b1) ;
    
    a2 = pc[1,2] /pc[2,2]
    b2 = target_mu[2] - target_mu[1] * a2
    line2_mask = VH .>  (VV.* a2 .+ b2) ;
    
    return line1_mask, line2_mask
end


function get_malhob(VV,VH,mu,sigma)
    malhob = cat(vec(VV) .- mu[1],vec(VH).- mu[2],dims=2)
    
    inv_sigma = sigma^-1;
    malhob = [transpose(malhob[i,:])*inv_sigma*(malhob[i,:]) for i in 1:length(VV)];
    malhob = reshape(malhob,size(VV));
    return malhob
end


function count_chi_limit(malhob,chi_limits)
    counts = zeros(Int64, length(chi_limits))
    
    counts[1] = sum(malhob .< chi_limits[1])
    
    for i in 2:(length(chi_limits))
        counts[i] = sum(chi_limits[i-1] .< malhob .< chi_limits[i])
    end
    
    return counts
end




################
function moasic2normal_view(moasic_view,meta)
    @assert moasic_view[1].start >= meta["burst_meta"]["first_line_mosaic"][1]
    @assert moasic_view[1].stop <= (meta["burst_meta"]["first_line_mosaic"][end] + meta["lines_per_burst"])

    
    start_burst = findlast(
    meta["burst_meta"]["first_line_mosaic"].<= moasic_view[1].start)

    lines_in_first_burst = moasic_view[1].start -meta["burst_meta"]["first_line_mosaic"][start_burst] +1
    line_start = lines_in_first_burst + (start_burst-1)* meta["lines_per_burst"]

    end_burst = 0
    if moasic_view[1].stop < meta["lines_per_burst"]
        end_burst = 1
    else
        end_burst = findlast(
            (meta["burst_meta"]["first_line_mosaic"] .+meta["lines_per_burst"] ).< moasic_view[1].stop)+1
    end

    lines_in_last_burst =  moasic_view[1].stop - meta["burst_meta"]["first_line_mosaic"][end_burst] +1
    line_end = lines_in_last_burst+ (end_burst-1)* meta["lines_per_burst"]
    return [line_start:line_end,moasic_view[2]]
end


function round_2_int(float_number,tol=0.1)
    res = round(Int64,float_number)
    @assert abs(res-float_number)<tol "Rounding not within tolerance"
return res
end


function split_large_view(large_view,meta_first,meta_second)
    ## TODO- Must be burst_wise
    c = 299792458
    azimuth_frequency =  meta_first["azimuth_frequency"]
    inv_range_pixel_spacing = (2*meta_first["range_sampling_rate"])/c
    delta_r0 = (meta_second["slant_range_time"]-meta_first["slant_range_time"])*c/2
    delta_sample = inv_range_pixel_spacing*delta_r0

    delta_sample = round_2_int(inv_range_pixel_spacing*delta_r0)

    over_lap = meta_first["samples_per_burst"]-delta_sample +1
    
    start_burst = ceil(Int64,large_view[1].start/meta_first["lines_per_burst"])
    last_valid_pixel = maximum(meta_first["burst_meta"]["last_valid_pixel"][start_burst])
    delta_minus = 1 + meta_first["samples_per_burst"] -last_valid_pixel 
    delta_plus = over_lap - delta_minus


    view_first = large_view[1],large_view[2].start:(meta_first["samples_per_burst"]-delta_minus)

    mosaic_view_first = SlcUtil.get_mosaic_view(meta_first,view_first)


    detlta_t0 = convert(Second,(meta_first["t_0"]-meta_second["t_0"]) ).value

    line_t_first_start =  (mosaic_view_first[1].start-1)/azimuth_frequency+meta_first["t_start"]
    line_second_start = round_2_int(1 + ((line_t_first_start+detlta_t0)-meta_second["t_start"])*azimuth_frequency)

    line_t_first_stop =  (mosaic_view_first[1].stop-1)/azimuth_frequency+meta_first["t_start"]
    line_second_stop= round_2_int(1 + ((line_t_first_stop+detlta_t0)-meta_second["t_start"])*azimuth_frequency)


    mosaic_view_second = [line_second_start:line_second_stop,
                        (1+delta_plus): (large_view[2].stop - delta_sample+1)]
    
    delta_line =  round_2_int( mosaic_view_first[1].start-mosaic_view_second[1].start)

    view_second = moasic2normal_view(mosaic_view_second,meta_second)
    return view_first,view_second,delta_line,delta_sample
end




function _meta_cal_datapath_pod(product_folders,polarization,subswath,pod_folder)
    
    file_paths = [Load.slc_paths(folder, polarization, subswath) for folder in product_folders]
    meta = [Load.slc_meta(path[2]) for path in file_paths]
    calibration = [Load.slc_calibration(file_paths[i][3],meta[i]["t_0"]) 
                            for i in 1:length(product_folders)];
    data_path = [path[1] for path in file_paths]
    
    pod_paths = [Load.pod_path(folder, pod_folder) for folder in product_folders]
    
    precise_orbit = [Load.precise_orbit(pod_paths[i],meta[i]["t_0"]) 
                            for i in 1:length(product_folders)]
    
    return meta, calibration, data_path ,precise_orbit
end


function _get_dem(meta,view,dem_path,nan_fill)
    footprint = SlcUtil.footprint(meta[1], view)
    latlon_window = ((minimum(footprint[1]), maximum(footprint[1])), (minimum(footprint[2]), maximum(footprint[2])))
    dem = 0 
    if split(dem_path,".")[end] == "jld"
        pad = 0.1
        dem = JLD.load(dem_path);
        idx1 =(latlon_window[1][1] - pad) .<dem["lat"].< (latlon_window[1][2] + pad)
        idx2 =(latlon_window[2][1] - pad) .<dem["lon"].< (latlon_window[2][2] + pad)
        dem = (dem["lat"][idx1], dem["lon"][idx2], dem["height"][idx1,idx2]);
        dem[3][dem[3].== -32768] .= nan_fill
    elseif split(dem_path,".")[end] == "tiff"
        dem = Load.dem(dem_path, latlon_window; nan_fill = nan_fill, padding=[90,90]);
    else
        @assert 1==2 "DEM format not reconised"
    end
    return dem 
end


function _slave_load_data(master_view,data_path,meta, precise_orbit,dem,calibration)
    #println("Load slave: $(data_path[2])")
    slave_data,flat,lut =  coregister_slave(master_view,
                                                data_path[2],
                                                meta, precise_orbit, dem)  # meta should be 2 Dict array
    slave_data, mosaic_view = SlcUtil.mosaic(slave_data, master_view, meta[1]);
    flat, mosaic_view = SlcUtil.mosaic(flat, master_view, meta[1]);
    slave_data = SlcUtil.calibrate_slave_data(slave_data, mosaic_view, lut, calibration[2]);
    return slave_data,flat,lut, mosaic_view
end

function _master_load_data(master_view,data_path,meta,calibration)
    #println("Load master: $(data_path[1])")
    master_data = Load.slc_data(data_path[1], master_view);
    master_data ,mosaic_view = SlcUtil.mosaic(master_data, master_view, meta[1]);
    master_data = SlcUtil.calibrate_data(master_data, Misc.flatten(mosaic_view...)..., calibration[1]);
    return master_data,mosaic_view
end


function _geo_lut_combi(lut_first,lut_second,delta_line,delta_sample,mosaic_view)
    lut_sample_strid = lut_first["master_sample"][2]-lut_first["master_sample"][1]
    sample_new_start = lut_first["master_sample"][end]+lut_sample_strid
    sample_new_temp = collect(sample_new_start :lut_sample_strid:  mosaic_view[2].stop);
    
    
    geo_lut  = Dict{String,Any}()
    geo_lut["master_sample"] =  cat(lut_first["master_sample"],sample_new_temp,dims=1)
    geo_lut["master_line"] =  lut_first["master_line"]

    row_second = lut_first["master_line"].-delta_line;
    col_second  = sample_new_temp.-delta_sample;
    
    dims_lut_first = (length(lut_first["master_line"]),length(lut_first["master_sample"]))
    dims_lut_second = (length(lut_second["master_line"]),length(lut_second["master_sample"]))
    dims_lut_new = (length(row_second),length(col_second))
    
    key_list = ["latitude","longitude","heights"]
    
    for key in key_list
        old_second = reshape(lut_second[key],dims_lut_second)
        new_second = Misc.interp_grid(lut_second["master_line"] ,lut_second["master_sample"],
                                                        old_second,row_second, col_second)
        commbined = hcat(reshape(lut_first[key],dims_lut_first),new_second)
        geo_lut[key] = commbined
    
    end
        
    return geo_lut
end

function _geo_lut!(lut)
    
    dims_lut = (length(lut["master_line"]),length(lut["master_sample"]))
    key_list = ["latitude","longitude","heights"]
    
    for key in key_list
        lut[key] = reshape(lut[key],dims_lut)
    end
    delete!(lut, "slave_line")
    delete!(lut, "slave_sample")
    
    return lut
end

function _get_data_mosaic_swaths(meta, calibration, data_path ,precise_orbit,
                                    product_folders, master_view, dem_path, 
                                  subswath, dem_nan, polarization)
    
    @assert subswath<3 "The view is outside the image"
    sec_meta,sec_calibration,sec_data_path, sec_precise_orbit = _meta_cal_datapath_pod(
                                                product_folders,polarization,subswath+1,pod_folder);

    view_first,view_second,delta_line,delta_sample = split_large_view(master_view,meta[1],sec_meta[1])

    dem_first = _get_dem(meta,view_first,dem_path,dem_nan)
    dem_second = _get_dem(sec_meta,view_second,dem_path,dem_nan);
    
    data_slave_first, flat_first, 
    lut_first, mosaic_view_first = _slave_load_data(view_first,data_path,meta, 
                                            precise_orbit,dem_first,calibration);

    data_slave_second, flat_second, 
    lut_second, mosaic_view_second = _slave_load_data(view_second,sec_data_path,sec_meta, 
                                            sec_precise_orbit,dem_second ,sec_calibration);


    master_data_first,temp = _master_load_data(view_first,data_path,meta,calibration)
    master_data_second,temp = _master_load_data(view_second,sec_data_path,sec_meta,sec_calibration)

    master_data = hcat(master_data_first,master_data_second);
    slave_data = hcat(data_slave_first,data_slave_second);
    flat =  hcat(flat_first,flat_second);

    mosaic_view = [mosaic_view_first[1],
                master_view[2]]

    geo_lut = _geo_lut_combi(lut_first,lut_second,delta_line,delta_sample,mosaic_view)


    return master_data, slave_data, flat, mosaic_view, geo_lut
end


function _get_data_single_swaths(meta, calibration, data_path ,precise_orbit,
                                    product_folders, master_view, dem_path, 
                                    dem_nan)

    dem = _get_dem(meta,master_view,dem_path,dem_nan)

    slave_data, flat, 
    lut, mosaic_view = _slave_load_data(master_view,data_path,meta, 
                                        precise_orbit,dem,calibration);


    master_data,temp = _master_load_data(master_view,data_path,meta,calibration)

    geo_lut = _geo_lut!(lut)

    return master_data, slave_data, flat, mosaic_view, geo_lut
end


function coherence_worker(product_folders, master_view, dem_path, 
                          subswath, dem_nan, pod_folder, polarization; 
                          kernel = ones(4,14))
    
    stride_line = floor(Int,size(kernel)[1]/2)
    stride_sample = floor(Int,size(kernel)[2]/2)

    meta, calibration, data_path ,precise_orbit = _meta_cal_datapath_pod(
                                                product_folders,polarization,subswath,pod_folder);
    #println("Data path master: $(data_path[1])")
    #println("Data path slave: $(data_path[2])")
    # initiliazise varibles
    master_data, slave_data,
    flat, mosaic_view, geo_lut = Nothing,Nothing,Nothing,Nothing,Nothing


    if meta[1]["samples_per_burst"] < master_view[2].stop 
        master_data, slave_data,
        flat, mosaic_view, geo_lut = _get_data_mosaic_swaths(meta, calibration, data_path 
            ,precise_orbit,product_folders, master_view, dem_path,
            subswath, dem_nan, polarization)
    else
        master_data, slave_data,
        flat, mosaic_view, geo_lut = _get_data_single_swaths(meta, calibration, data_path ,precise_orbit,
                                        product_folders, master_view, dem_path, 
                                        dem_nan)
    end;


    complex_coherence, master_intensity, slave_intensity, lines, samples = SlcUtil.complex_coherence(
                        master_data, slave_data, flat, kernel, mosaic_view);


    # Subsample 
    lines = lines[1:stride_line:end]
    samples = samples[1:stride_sample:end]
    master_intensity = master_intensity[1:stride_line:end,1:stride_sample:end]
    slave_intensity = slave_intensity[1:stride_line:end,1:stride_sample:end]
    complex_coherence = complex_coherence[1:stride_line:end,1:stride_sample:end]

    # Check that the size matches
    @assert length(master_intensity) > 0
    @assert sum(size(master_intensity) .== size(complex_coherence)) ==2
    @assert sum(size(slave_intensity) .== size(complex_coherence)) == 2
    @assert size(master_intensity)[1] == length(lines)
    @assert size(master_intensity)[2] == length(samples)

    return complex_coherence, master_intensity, slave_intensity, lines, samples, geo_lut
end



function _save_to_jld(file_name,data,folder)
    path = joinpath(folder,file_name)*".jld"
    #println("Eltype : $(eltype(data))")
    if isa(data,Dict)
        JLD.save(path,data)
    else
        JLD.save(path, "data",data)
    end
end




function _get_llh(lut,lines,samples,result_folder="")
    llh = Dict{String,Any}()
    key_list = ["latitude","longitude","heights"]
    
    for key in key_list
        llh[key] =  Misc.interp_grid(lut["master_line"],
                                lut["master_sample"],
                                lut[key],
                                collect(lines), 
                                collect(samples));
    end
    
    if length(result_folder) >0
        _save_to_jld("coordinates",llh,result_folder)
    end

    return hcat(reshape(llh["latitude"],:),
                reshape(llh["longitude"],:),
                reshape(llh["heights"],:)) 
end


function _get_target_line_sample(super_llh, meta, safe_folder, pod_folder)
    pod_path = Load.pod_path(safe_folder, pod_folder)
    pod = Load.precise_orbit(pod_path,meta["t_0"])
    line_sample = Geometry.to_line_sample(super_llh[:,1:2], super_llh[:,3], pod..., meta)
    return line_sample[:,1],line_sample[:,2]
end

function _get_meta(safe_folder, polarization, subswath)
    paths = Load.slc_paths(safe_folder, polarization, subswath)
    return Load.slc_meta(paths[2])
end

function _get_master_view(target_line,target_sample,master_meta,pad=50)
    
    start_line = floor(Int64,minimum(target_line)-pad)
    end_line = ceil(Int64,maximum(target_line)+pad)

    start_sample = floor(Int64,minimum(target_sample)-pad)
    end_sample = ceil(Int64,maximum(target_sample)+pad)

    master_view_mosaic = [start_line:end_line,start_sample:end_sample]
    master_view = moasic2normal_view(master_view_mosaic,master_meta)
    return master_view
end

function _resamples_to_target(coherrence,master,slave,lines,samples,target_line,target_sample,super_dims)

    start_lines,step_lines,stop =  start_step_stop(lines)
    index1 = (target_line .- start_lines)./step_lines

    start_samples,step_samples,stop =  start_step_stop(samples)
    index2 = (target_sample .- start_samples)./step_samples
    
    dummy_view = [0:2,0:2]
    
    master_re= Misc.resample(dummy_view, master, index1, index2) 
    master_re = reshape(master_re,super_dims)
    
    slave_re= Misc.resample(dummy_view, slave, index1, index2) 
    slave_re = reshape(slave_re,super_dims)
    
    coherrence_re= Misc.resample(dummy_view, coherrence, index1, index2) 
    coherrence_re = reshape(coherrence_re,super_dims)
    
    return coherrence_re,master_re,slave_re
end

function save_parameters(result_folder,files, super_master_view, 
        subswath,  dem_nan, super_master_index, kernel )
    folder_names_path = joinpath(result_folder,"function_parameters.txt")   
    open(folder_names_path, "w") do io
        # Save info about inputs
        write(io, "super master index: $super_master_index\n")
        write(io, "master_view: $super_master_view\n")
        write(io, "subswath: $subswath\n")
        write(io, "dem_nan: $dem_nan\n")
        write(io, "Kernel size: $(size(kernel))\n")
        write(io, "\n")
        write(io, "files: \n")
        for elem in files
            write(io,"$(splitdir(elem)[2])\n")
        end
    end;
end



function coherence_factory_seq(files, super_master_view, dem_path, subswath, result_folder, dem_nan, 
        pod_folder; super_master_index = 1, polarization = ["VV","VH"],save_intensity = true, kernel = ones(4,14))

    n_files = length(files)

    if super_master_index == "end"
        super_master_index = n_files
    end

    @assert super_master_index<=n_files
    save_parameters(result_folder,files, super_master_view, 
            subswath,  dem_nan, super_master_index, kernel )

    master_list= collect(1:(n_files-1))
    slave_list = collect(master_list).+1
    intensity_saved = zeros(Bool,n_files)

    # flip first coherrence pair is need
    if super_master_index == n_files
        master_list[end]=n_files
        slave_list[end]=n_files-1
    end

    pair_index = findfirst(master_list.==super_master_index)

    ## process first coherence
    master_index = master_list[pair_index]
    slave_index = slave_list[pair_index]
    product_folders = [files[master_index],files[slave_index]]
    println("coherence start")
    println("M:$(splitdir(product_folders[1])[2])")
    println("S:$(splitdir(product_folders[2])[2])")

    super_dims = 0

    lines_super,samples_super, super_geo_lut = [Nothing,Nothing,Nothing]

    for pol in polarization
        println(pol)
        complex_coherence, master_intensity, 
        slave_intensity,lines_super,
        samples_super, super_geo_lut  = coherence_worker(product_folders, super_master_view, dem_path, 
                                                      subswath, dem_nan, pod_folder, pol);
        super_dims = size(master_intensity)
        # save the results
        file_names = _file_name_generater(product_folders, pol)
         _save_to_jld(file_names[1] ,complex_coherence,result_folder)
        if save_intensity
             _save_to_jld(file_names[2] ,master_intensity,result_folder)
             _save_to_jld(file_names[3] ,slave_intensity,result_folder)
        end

    end

    ## update intensity_saved 
    intensity_saved[master_index] = true
    intensity_saved[slave_index] = true
    println("")

    ## remove first coherence pair from list
    deleteat!(master_list,pair_index)
    deleteat!(slave_list,pair_index)


    println("Get super llh")
    super_llh = _get_llh(super_geo_lut,lines_super,samples_super,result_folder);
    println("")


    for i in 1:length(master_list)

        ## get preoducts
        master_index = master_list[i]
        slave_index = slave_list[i]
        product_folders = [files[master_index],files[slave_index]]

        println("Coherence loop: $i/$(length(master_list))")
        println("M:$(splitdir(product_folders[1])[2])")
        println("S:$(splitdir(product_folders[2])[2])")

        println("Get taget line sample")
        #find target line and sample
        master_meta = _get_meta(product_folders[1], polarization[1], subswath)
        target_line,target_sample = _get_target_line_sample(super_llh,master_meta,
                                                        product_folders[1], pod_folder);
        ## get master view
        master_view = _get_master_view(target_line,target_sample,master_meta);

        for pol in polarization
            println(pol)
            # compute coherence
            complex_coherence, master_intensity, 
            slave_intensity,lines,samples, geo_lut  = coherence_worker(product_folders, master_view, dem_path, 
                                                      subswath, dem_nan, pod_folder, pol);
            # resample
            complex_coherence, master_intensity, 
                                    slave_intensity = _resamples_to_target(complex_coherence, master_intensity, 
                                                    slave_intensity, lines,samples,target_line,target_sample,super_dims)

            # save files
            file_names = _file_name_generater(product_folders, pol)
             _save_to_jld(file_names[1] ,complex_coherence,result_folder)
            if save_intensity
                if !intensity_saved[master_index]
                     _save_to_jld(file_names[2] ,master_intensity,result_folder)
                end
                if !intensity_saved[slave_index]
                     _save_to_jld(file_names[3] ,slave_intensity,result_folder)
                end
            end
        end

        ## update intensity_saved 
        intensity_saved[master_index] = true
        intensity_saved[slave_index] = true
        println("")
    end
    return 1
end



#################




function grad_eval(f,x,delta)
    dfdx = similar(x)
    f_0 = f(x...)
    for i in 1:length(x)
        xp = copy(x)
        xb = copy(x)
        xp[i] += delta[i]
        xb[i] -= delta[i]
        dfdx[i] = (f(xp...)-f(xb...))/(2*delta[i])
    end
    return dfdx, f_0
end


function grad_decent(f, x_0 ; no_step_max =3, iter_max = 50,df_tol = 10^-14 ,dfdx_tol = 10^-16, 
        step_size = 10.0 .^(2:6),no_step_count = -2,delta = 10^-4,debug=false)

    x_i = x_0
    iter_count = 0
    n_steps = 0
    delta = ones(length(x_i))*delta
    
    #test if better or the first
    while (no_step_count<no_step_max)
        iter_count += 1;
        dfdx, f_i =  grad_eval(f,x_i,delta)
        
        if debug
            println("")
            println("Iterration : $iter_count")
            println("x_i: $x_i, f_i: $f_i")
             println("dfdx: $dfdx")
        end
        
        
        if sum(abs.(dfdx)) > dfdx_tol
            x_old = x_i
            
            #check different steps sizes and choo the best
            f_steps = [f( (x_i .- dfdx.*elem)...) for elem in step_size]
            min_idx = argmin(f_steps)
            
            # Chec decrease in loose function
            df = f_i - f_steps[min_idx];
            
            if df < df_tol
                #decreasing  delta becuse no imporvment found in the direction of the gradiant
                delta .*= 0.1
                no_step_count += 1    
                if debug
                    println("No change, Decrease delta, no_step: $no_step_count")
                end
            else
                # update x_i with optima step and reset no_step_count
                x_i = x_i.- dfdx.* step_size[min_idx]
                no_step_count = minimum([no_step_count,0])
                 n_steps +=1
                if debug
                    println("Step x, step_idx:$min_idx,  no_step: $no_step_count")
                end
            end
        end
        
        # increase delta
        delta[dfdx .< dfdx_tol] *= 10
       

        if (iter_count > iter_max)
            println("MAX NUMBER OF ITERATIONS REACHED")
            break
        end

    end
    
    dx = x_i .- x_0
    
    if transpose(dx)*dx < 10^-5
        println("WARNING no step taken")
    end

    return x_i, n_steps

end




function _sort_jld_files(file_list,date_position )
    name_list = [split(elem,".")[1] for elem in file_list]
    date_string = [split(elem,"_")[date_position] for elem in name_list]
    date = [Dates.Date(parse.(Int, [elem[1:4], elem[5:6], elem[7:8]])...) for elem in date_string]
    return file_list[sortperm(date)][end:-1:1]
end

function _sort_prossed_files(data_folder,sort_master = true)
    files = readdir(data_folder)


    files = [elem for elem in files if length(elem)>1]
    files = [elem for elem in files if length(split(elem,"_"))>3]

    coherence_idx = [split(elem,"_")[2]=="coh" for elem in files]
    coherence_files = files[coherence_idx ]  
    if sort_master 
        coherence_files = _sort_jld_files(coherence_files,6)
    else
        _sort_jld_files(coherence_files,8)
    end

    coherence_VV_files = [ elem for elem in coherence_files if split(elem,"_")[4]=="VV"]
    coherence_VH_files = [ elem for elem in coherence_files if split(elem,"_")[4]=="VH"]                    

    VV_idx = [(split(elem,"_")[1]=="sigma") & (split(elem,"_")[3]=="VV") for elem in files];
    VV_files = files[VV_idx] 
    VV_files = _sort_jld_files(VV_files,4)

    VH_idx = [(split(elem,"_")[1]=="sigma") & (split(elem,"_")[3]=="VH") for elem in files];
    VH_files= files[VH_idx];                 
    VH_files = _sort_jld_files(VH_files,4)
                                                
    return VV_files,VH_files,coherence_VV_files,coherence_VH_files 
end
                                                
function _load_jld(file_list,data_folder)
    path_names = [joinpath(data_folder,elem) for elem in file_list]
    return [JLD.load(elem,"data") for elem in path_names]
end



shapely_geometry = PyCall.pyimport("shapely.geometry")
skimage_meas = PyCall.pyimport("skimage.measure");
pandas = PyCall.pyimport("pandas")
geopandas = PyCall.pyimport("geopandas")
ndimage = PyCall.pyimport("scipy.ndimage")
skimage_morph = PyCall.pyimport("skimage.morphology");


zip_folder(dir_path) = run(`zip -q -j -r $(dir_path*(".zip")) $dir_path`)# zip files

flat(band,test_area) = reshape(band[test_area...],:)

LsqFit.@. bimodal_gauss_model(x, p) = p[1]*exp(-0.5*((x-p[3])/p[5])^2) +  p[2]*exp(-0.5*((x-p[4])/p[6])^2)

LsqFit.@. gauss_model(x, p) = p[1]*exp(-0.5*((x-p[2])/p[3])^2) 


function get_edges(mask,line_width=2)
    
    h_kernel = [1 1 1 ; 0 0 0 ; -1 -1 -1]
    v_kernel = [1 0 -1 ; 1 0 -1 ; 1 0 -1]
    
    res = Misc.fastconv(convert.(Int64,mask),h_kernel).^2 .+ Misc.fastconv(convert.(Int64,mask),v_kernel).^2;
    res =  res[2:end-1, 2:end-1];
    res = res.>0
    
    for i in 1:line_width
        res = skimage_morph.binary_dilation(res)
    end
    return res
end

function sse_water_fit(flood_band,change_band,seed_mask,bm_mask, p_water,w_sum,edges,y, rg_thresholds)
    #use region growinf
    rg_mask = (flood_band .<rg_thresholds[1]) .& (change_band.<rg_thresholds[2]) .| seed_mask 
    rg_result, steps = region_growing(seed_mask,rg_mask);
    
    # select the water pixelss in the selected tiles
    data = reshape(flood_band,:)[reshape(rg_mask.&bm_mask,:)];
    
    # Compare histogram with the emepircal in p_water
    h = StatsBase.fit(StatsBase.Histogram, data,edges)
    w_sel = h.weights./w_sum;
    w_water_est = gauss_model(y,p_water)
    
    # return sum of squred errors
    SSE = (w_sel .- w_water_est)'*(w_sel .- w_water_est)
    return SSE
end


function sse_water_fit2(flood_band,change_band,ref_band,y_seed,bm_mask, p_water,w_sum,edges,y, rg_thresholds)
    #use region growinf
    
    rg_mask = (flood_band .<rg_thresholds[1]) .& (change_band.<rg_thresholds[2])
    rg_flood, steps = region_growing(flood_band .< y_seed,rg_mask);
    
    
    rg_mask_ref = (ref_band .<rg_thresholds[1])
    rg_ref, steps = region_growing(ref_band .<y_seed,rg_mask_ref);
    
    # Remove flase positives and permant water.
    rg_mask= rg_flood .& (rg_ref .!=true);
    
    # select the water pixelss in the selected tiles
    data = reshape(flood_band,:)[reshape(rg_mask.&bm_mask,:)];
    
    # Compare histogram with the emepircal in p_water
    h = StatsBase.fit(StatsBase.Histogram, data,edges)
    w_sel = h.weights./w_sum;
    w_water_est = gauss_model(y,p_water)
    
    # return sum of squred errors
    SSE = (w_sel .- w_water_est)'*(w_sel .- w_water_est)
    return SSE
end


"""
find_y_seed(p_fit,y,min_ratio = 0.99)

find the y_seed

"""
function find_y_seed(p_fit,y,min_ratio = 0.99)
    p_water = p_fit[[1,3,5]]
    p_nonwater = p_fit[[2,4,6]]
    ratio = gauss_model(y,p_water)./bimodal_gauss_model(y, p_fit)
    # Select y_seed where the single gauuss divege from the bimodal_gauss_model
    y_seed = y[y.>p_water[2]][ findfirst(ratio[y.>p_water[2]].<min_ratio)]
    return y_seed
end



function HSBA_floodmask(flood_band,ref_band)
    change_band = flood_band.-ref_band;
    
    # Find the tiles to fit bimodal
    bm_mask_flood = find_bimodal_tiles(flood_band);
    bm_mask_change = find_bimodal_tiles(change_band);
    bm_mask = bm_mask_change.&bm_mask_flood
    
    # fit bimodel
    data = reshape(flood_band,:)[reshape(bm_mask,:)]
    p_fit,y,w,edges, w_sum = fit_bimodal_gauss(data,round(Int64,length(data)/50))
    
    # Find seed pixels
    y_seed =find_y_seed(p_fit,y)
    seed_mask = flood_band .<y_seed
    
    # optimize y_RG and delta_sigma tresholds
    t_0 = [y_seed+1, -1]
    res = Optim.optimize(
        t -> sse_water_fit(flood_band,change_band,seed_mask,bm_mask, p_fit[[1,3,5]],w_sum,edges,y, t), 
        t_0; autodiff = :forward)
    
    # use optimize values to find flood mask
    rg_thresholds = res.minimizer
    rg_mask = (flood_band .<rg_thresholds[1]) .& (change_band.<rg_thresholds[2]) .| seed_mask 
    flood_mask, steps = region_growing(seed_mask,rg_mask);
    
    # Make region growing og refference image to find Permant water and false positives
    seed_mask_ref = ref_band .<y_seed
    rg_mask_ref = (ref_band .<rg_thresholds[1])
    ref_mask, steps = region_growing(seed_mask_ref,rg_mask_ref);
    
    # Remove flase positives and permant water.
    final_mask = flood_mask .& (ref_mask .!=true);
    
    return final_mask, [y_seed,rg_thresholds...]
end


"""
Test function to use the HSBA algorithmn but with given thresholds

"""
function No_HSBA_floodmask(flood_band,ref_band,thresholds)
    change_band = flood_band.-ref_band;
    seed_mask = flood_band .<thresholds[1]
    
    # use optimize values to find flood mask
    rg_mask = (flood_band .<thresholds[2]) .& (change_band.<thresholds[3]) .| seed_mask 
    flood_mask, steps = region_growing(seed_mask,rg_mask);
    
    # Make region growing og refference image to find Permant water and false positives
    seed_mask_ref = ref_band .<thresholds[1]
    rg_mask_ref = (ref_band .<thresholds[2])
    ref_mask, steps = region_growing(seed_mask_ref,rg_mask_ref);
    
    # Remove flase positives and permant water.
    final_mask = flood_mask .& (ref_mask .!=true);
    
    return final_mask, flood_mask, ref_mask
end


function region_growing(seed_mask,region_grow_mask ; max_inter=100, tol=0)
    res = seed_mask
    region_grow_mask = region_grow_mask .| seed_mask
    steps = 1
    diff = tol +5
    while diff> tol
        
        new_mask = _region_grow_step(res, region_grow_mask)
        diff = sum(new_mask .!= res)
        res= new_mask
        
        if steps >= max_inter
            println("Max Iter reached:  region_growing()")
            break
        end
        steps += 1
    end
    return res,steps
end


function _region_grow_step(seed_mask, region_grow_mask)
    res = skimage_morph.binary_dilation(seed_mask)
    return res.& region_grow_mask 
end



function rg_gif(file_path, seed_mask, region_grow_mask,img,n)
    res = copy(seed_mask);
    gif_file = Array{RGB{Float32}}(undef, size(res)..., n+1);
    gif_file[:,:,1] .= add_mask(img, res ,(0,0,1))

    for i = 1:n
        res = _region_grow_step(res, region_grow_mask)
        gif_file[:,:,i+1] .= add_mask(img, res ,(0,0,1));
        end;

    FileIO.save(file_path, gif_file)
end

function find_bimodal_tiles(data; N_limit = (10^3, 3*10^5), bin_size=50,conditions=(2,0.99,0.1),max_iter=10000)
    
    # initialize mask
    dims = size(data)
    N =  dims[1]*dims[2]
    BM_mask = zeros(Bool,dims...)
    
    # If area is to small return false and exit function
    # serves as a lower limit for the recursive function
    if N < N_limit[1]
        return BM_mask
    end
    
    # If N is under max limit check for bimodal 
    bm_tile = false
    if N < N_limit[2]
        n_bins = round(Int64,N/bin_size)
        bm_tile, cons = bimodal_tile(reshape(data,:),n_bins,conditions=conditions,max_iter=max_iter) 
    end
    
    
    if bm_tile
        BM_mask .= true
    else
        # If tile not bimodal, Split it in 4 and try again
        row_split = floor(Int64,dims[1]/2)
        col_split = floor(Int64,dims[2]/2)
        
        tiles = [(1:row_split,1:col_split),((row_split+1):dims[1],1:col_split),
            (1:row_split,(col_split+1):dims[2]),((row_split+1):dims[1],(col_split+1):dims[2])]
        
        for elem in tiles
            BM_mask[elem...] .= find_bimodal_tiles(data[elem...], 
                N_limit =N_limit, bin_size=bin_size,conditions=conditions,max_iter=max_iter)
        end
        
    end
    
    return BM_mask
end
    


function bimodal_tile(data,n_bins;max_iter=10000,conditions=(2,0.99,0.1))
    p_fit,y,w,edges,w_sum = fit_bimodal_gauss(data,n_bins,max_iter)
    
    ad = sqrt(2)*abs(p_fit[3]-p_fit[4])/sqrt(p_fit[5]^2+p_fit[6]^2)
    bc = sum(sqrt.(w[w .!=0]).*sqrt.(bimodal_gauss_model(y[w .!=0], p_fit)))
    sr = minimum([p_fit[1]*p_fit[5],p_fit[2]*p_fit[6]])/maximum([p_fit[1]*p_fit[5],p_fit[2]*p_fit[6]])

    BM_tile = (ad>conditions[1])&(bc>conditions[2])&(sr>conditions[3])
    return BM_tile, [ad,bc,sr]
end


function fit_bimodal_gauss(data,n_bins,maxIter=1000)
    # Get histogram
    h = StatsBase.fit(StatsBase.Histogram, data,nbins=n_bins)
    w = h.weights
    edges = collect(h.edges[1])
    y = (edges[2:end] +  edges[1:end-1])/2;
    w_sum = sum(w)
    # normalise
    w = w./w_sum
    
    # find  OTSU threshold
    # see https://ieeexplore-ieee-org.proxy.findit.dtu.dk/document/4310076
    mu_k = cumsum(w.*y)
    mu_t = mu_k[end]
    omega_k = cumsum(w)
    sigma2_b = (mu_t.*omega_k .- mu_k).^2 ./(omega_k .*(1 .- omega_k ))
    sigma2_b[end] = 0  # last element is NaN but should be zero
    y_ot = y[argmax(sigma2_b)]

    # initialize parameters using OTSU threshold
    p0 = zeros(Float64,6)
    
    p0[3] = Statistics.mean(data[data.<y_ot])
    p0[5] = Statistics.std(data[data.<y_ot])
    p0[1] = w[findfirst(y.>p0[3])]


    p0[4] = Statistics.mean(data[data.>y_ot])
    p0[6] = Statistics.std(data[data.>y_ot])
    p0[2] = w[findfirst(y.>p0[4])];
    
    # Fit function
    p_fit = LsqFit.curve_fit(bimodal_gauss_model, y, w, p0; autodiff=:forwarddiff,maxIter=maxIter).param
    
    return p_fit, y, w,edges,w_sum

end



function db_scale_img(img , min , max)
    log_img = (10*log10.(img).-min)./(max-min)
    log_img[log_img.>1] .= 0.999
    log_img[log_img.<0] .= 0.001
    return log_img
end


function db_scale_img_diff(img1,img2 , min , max)
    log_img = 10*log10.(img1).-10*log10.(img2)
    log_img = (log_img.-min)./(max-min)
    log_img[log_img.>1] .= 0.999
    log_img[log_img.<0] .= 0.001
    return log_img
end

function scale_img(img , min , max)
    scale_img = (img.-min)./(max-min)
    scale_img[scale_img.>1] .= 0.999
    scale_img[scale_img.<0] .= 0.001
    return scale_img
end


function scale_img(img)
    min = minimum(reshape(img,:))
    max = maximum(reshape(img,:))
    return scale_img(img , min , max)
end



"""
    add_mask(img,mask,color=(1,0,0))

    Add mask to a 2d color image. The mask is given the color stated in the input.
    the mask is a 2d bool mask

"""
function add_mask(img,mask,color=(1,0,0))
    color = Colors.RGB{Float32}(color...)
    img_cop = copy(img)
    img_cop[mask] .= color
    return img_cop
end


function pretty_img(bands,min,max,k=1/1.4)
    
    pre = 0
    if length(bands) ==3
        pre = db_scale_img((bands[2] .+bands[3])./2,min,max) 
    elseif length(bands) ==2
        pre = db_scale_img(bands[2],min,max) 
    end
    co = db_scale_img(bands[1],min,max) 
    co = co.^k
    pre = pre.^k
    return Colors.RGB{Float32}.(pre,co,co);
end



"""
    find_contours(img,threshold = 0.5)

    Find contours in the image. The boundaries are padded with zero to give a har boundary.
    
    Output is a list of contours. Each contour is a 2d array with each row being a point

"""
function find_contours(img,threshold = 0.5)
    temp = convert.(Float32,img)
    # Add hard boundary
    temp = hcat(zeros(size(temp)[1],1),temp,zeros(size(temp)[1],1))
    temp = vcat(zeros(1,size(temp)[2]),temp,zeros(1,size(temp)[2]))
    # Find contours
    contours = skimage_meas.find_contours(temp, 0.5)
    # The paddin os zeros make the python index of temp macth julia index om img
    return contours
end


"""
    close_contours!(contours)

    Make sure that first and last point in the contour is the same. 
    If that is not allready the case the first point is appended to the countour 
     so it is also the new end point
"""
function close_contours!(contours)
    not_closed = [sum(elem[1,:].==elem[end,:]) !=2 for elem in contours]
    if sum(not_closed) != 0
        println("$(sum(not_closed)) Polygons have been closed")
        contours[not_closed] = [vcat(elem,elem[1,:]') for elem in contours[not_closed]]
    end
end


"""
    start_step_stop(range)
    
    get the start, step and the stop value of a range
"""
function start_step_stop(range)
    if isa(range,StepRangeLen)
        start = range.ref.hi
        step = range.step.hi
        stop = range.step.hi* (range.len-1) + range.ref.hi
    elseif isa(range,StepRange)
        start = range.start
        step = range.step
        stop = range.stop
    end
    return start,step,stop
end


"""
    contours_to_coordinates(contours,line_sample,lut)
    
    Convert contours from row and column to coordinates
"""
function contours_to_coordinates(contours,line_sample,lut)
    contours_coords =  deepcopy(contours)
    
    start_line_contour ,step_line_contour ,stop_line_contour  = start_step_stop(line_sample["lines"])
    start_sample_contour ,step_sample_contour ,stop_sample_contour  = start_step_stop(line_sample["samples"])
    
    start_line_lut= lut["master_line"][1]
    step_line_lut= lut["master_line"][2]-lut["master_line"][1]
    
    start_sample_lut= lut["master_sample"][1]
    step_sample_lut= lut["master_sample"][2]-lut["master_sample"][1]
    
    start_sample_contour ,step_sample_contour ,stop_sample_contour  = start_step_stop(line_sample["samples"])

    lut_dims = (length(lut["master_line"]),length(lut["master_sample"]))
    lat_grid = reshape(lut["latitude"],lut_dims)
    lon_grid = reshape(lut["longitude"],lut_dims);
    
    for i in 1:length(contours)
        ## contours in line and sample
        index_1 = (contours[i][:,1].-1) .*step_line_contour .+start_line_contour
        index_2 = (contours[i][:,2].-1) .*step_sample_contour .+start_sample_contour
        
        ## countours in row, coulmn of lut grid
        index_1 = (index_1 .- start_line_lut)./step_line_lut
        index_2 = (index_2 .- start_sample_lut)./step_sample_lut

        contours_coords[i][:,1] = ndimage.map_coordinates(lat_grid,[index_1,index_2])
        contours_coords[i][:,2] = ndimage.map_coordinates(lon_grid,[index_1,index_2])
    end
    
    return contours_coords
end


shapely_polygon(contour) = shapely_geometry.Polygon([(contour[i,2],contour[i,1]) for i in 1:size(contour)[1]])


"""
    geopandas_df(data_dict,geometry_name,crs = "EPSG:4326")

    Makes a geopandas data frame.
    
"""
function geopandas_df(data_dict,geometry_name,crs = "EPSG:4326")
    df = pandas.DataFrame(data_dict)
    gdf = geopandas.GeoDataFrame(df, geometry=geometry_name)
    gdf.crs = Dict("init" => crs)
    return gdf
end



function quicksave_shape(name,folder,mask,line_sample,lut,close=true)
    contours = find_contours(mask)
    if close
        close_contours!(contours)
    end
    contours = contours_to_coordinates(contours,line_sample,lut)
    contours = shapely_polygon.(contours)
    
    number = collect(1:length(contours))
    data_dict = Dict("Number" => number, "Polygon"=> contours)
    gdf = geopandas_df(data_dict,"Polygon")
    save_shape_zip(gdf,name,folder)
    return  gdf
end


"""
    save_shape_zip(gdf,name,folder)

    Saves a geopandas data frame as shape files and then compres it to a .zip
    
"""
function save_shape_zip(gdf,name,folder)
    
    dir_path = joinpath(folder,name)
    zip_path = dir_path*".zip"
    shapepath = joinpath(dir_path,name)*".shp"
    
    run(`mkdir $dir_path`) #Make folder
    gdf.to_file(shapepath) # save files
    run(`zip -q -j -r $zip_path $dir_path`)# zip files
    run(`rm -r $dir_path`) # remod files
end



"""
    poly_slize(data,polygon)

    data is a 2d array with values and polygon is a 2d array where every row is a point.
    Returns a list with all the values from data inside the polygon
"""
function poly_slize(data,polygon)

    dims = size(data)

    index_in, window, window_dims = _get_windows_nindex(polygon,dims) 
    
    data_slice = data[window...]
    
    return reshape(data_slice,:)[index_in]
end


"""
    poly_to_mask(polygon_list,dims)

    Vector to raster. Makes a 2d array with oall the polygons.
    The polygon number is used as fill value
"""
function poly_to_mask(polygon_list,dims)
    
    mask = zeros(Int64,dims...)
    
    for i = 1 : length(polygon_list)
        index_in, window, window_dims = _get_windows_nindex(polygon_list[i],dims) 
        mask[window...] .= reshape(index_in,window_dims).*i
    end

    return mask
end


"""
    _get_windows_nindex(polygon,dims) 

    help function
"""
function _get_windows_nindex(polygon,dims) 

    row_max = minimum([ceil(Int64,maximum(polygon[:,1])),dims[1]])
    row_min = maximum([floor(Int64,minimum(polygon[:,1])),1])

    col_max = minimum([ceil(Int64,maximum(polygon[:,2])),dims[2]])
    col_min = maximum([floor(Int64,minimum(polygon[:,2])),1])


    window = (row_min:row_max,col_min:col_max)
    window_dims = (row_max-row_min+1,col_max-col_min+1)
    
    col, row = Misc.flatten(window...)
    index_in = [inpolygon([col[i],row[i]],polygon) for i =1:(window_dims[1]*window_dims[2])]
    return index_in, window, window_dims
end





######################## Function not made by Simon Kok Lupemba #################################



### Functions by Eigil Lippert ###


function extract_datetime(SAFE_path; start_date=true)
    extract_SAFE_name = split(SAFE_path, "/")[end]
    if start_date
        date_string = split(extract_SAFE_name, "_")[6]
    else
        date_string = split(extract_SAFE_name, "_")[7]
    end
    year = date_string[1:4]
    month = date_string[5:6]
    day = date_string[7:8]
    hour = date_string[10:11]
    minute = date_string[12:13]
    second = date_string[14:end]
    date_int = parse.(Int, [year, month, day, hour, minute, second])
    return DateTime(date_int...)
end

function days_between_acquisitions(date1, date2)
    return Dates.value(Date(date1) - Date(date2))
end

function _file_name_generater(product_folders, polarization)
    ID = Dict{String,String}()
    master_satellite = string(split(split(product_folders[1], "/")[end], "_")[1][end])
    slave_satellite = string(split(split(product_folders[2], "/")[end], "_")[1][end])

   
    # compute days between acquisitions:
    master_date = replace.(string(Date(extract_datetime(product_folders[1]))),  "-" => "")
    slave_date = replace.(string(Date(extract_datetime(product_folders[2]))),  "-" => "")
    days_between_acq = string(abs(days_between_acquisitions(extract_datetime(product_folders[1]), extract_datetime(product_folders[2]))))*"d"

    # define id strings on format:            
    master_name = "sigma_"*"S1"*master_satellite*"_"*polarization*"_"*master_date
    slave_name = "sigma_"*"S1"*slave_satellite*"_"*polarization*"_"*slave_date
    coherence_name = days_between_acq*"_"*"coh_"*master_satellite*slave_satellite*"_"*polarization*"_"*"M_"*master_date*"_"*"S_"*slave_date
    return coherence_name, master_name, slave_name
end



###### Function from the internet   ####


"""
inpolygon(p, poly)
Modified from https://github.com/JuliaGeometry/PolygonOps.jl/blob/master/src/inpolygon.jl
check the membership of `p` in `poly`. Works invariant of winding order.
Returns:
- in = 1
- on = 1
- out = 0
Based on the algorithm by Hao and Sun :
https://www.researchgate.net/publication/328261689_Optimal_Reliable_Point-in-Polygon_Test_and_Differential_Coding_Boolean_Operations_on_Polygons

Modified from 
Copyright (c) 2019 steve <kd2cca@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
function inpolygon(p, poly)
    k = 0

    xp = p[1]
    yp = p[2]

    PT = eltype(p)

    for i = 1:size(poly)[1]-1
        v1 = poly[i,2] - yp
        v2 = poly[i+1,2] - yp

        if v1 < zero(PT) && v2 < zero(PT) || v1 > zero(PT) && v2 > zero(PT)
            continue
        end

        u1 = poly[i,1] - xp
        u2 = poly[i+1,1] - xp

        f = (u1 * v2) - (u2 * v1)

        if v2 > zero(PT) && v1 <= zero(PT)
            if f > zero(PT)
                k += 1
            elseif iszero(f)
                return true
            end
        elseif v1 > zero(PT) && v2 <= zero(PT)
            if f < zero(PT)
                k += 1
            elseif iszero(f)
                return true
            end
        elseif iszero(v2) && v1 < zero(PT)
            iszero(f) && return true
        elseif iszero(v1) && v2 < zero(PT)
            iszero(f) && return true
        elseif iszero(v1) && iszero(v2)
            if u2 <= zero(PT) && u1 >= zero(PT)
                return true
            elseif u1 <= zero(PT) && u2 >= zero(PT)
                return 
            end
        end
    end

    iszero(k % 2) && return false
    return true
end

