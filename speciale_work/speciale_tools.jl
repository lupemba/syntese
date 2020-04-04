include("../ToolBox/ToolBox.jl")
using .ToolBox
using Colors
import PyCall
import StatsBase;
import LsqFit
import Statistics
import FileIO
import Optim

shapely_geometry = PyCall.pyimport("shapely.geometry")
skimage_meas = PyCall.pyimport("skimage.measure");
pandas = PyCall.pyimport("pandas")
geopandas = PyCall.pyimport("geopandas")
ndimage = PyCall.pyimport("scipy.ndimage")
skimage_morph = PyCall.pyimport("skimage.morphology");


flat(band,test_area) = reshape(band[test_area...],:)

LsqFit.@. bimodal_gauss_model(x, p) = p[1]*exp(-0.5*((x-p[3])/p[5])^2) +  p[2]*exp(-0.5*((x-p[4])/p[6])^2)

LsqFit.@. gauss_model(x, p) = p[1]*exp(-0.5*((x-p[2])/p[3])^2) 


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


function region_growing(seed_mask,region_grow_mask ; max_inter=1000, tol=0)
    res = seed_mask
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
    p_fit = LsqFit.curve_fit(bimodal_gauss_model, y, w, p0; autodiff=:forwarddiff,maxIter=1000).param
    
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
    pre = db_scale_img((bands[2] .+bands[3])./2,min,max) 
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

