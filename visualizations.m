% ------------------------------------------------------------------ %
% --------------------- Main function ------------------------------ %
% ------------------------------------------------------------------ %

% Call the execution
visualization_main_exec("2007_esacci/");

function [] = visualization_main_exec(data_src)
    out_dir = "./dtvm_outputs/";
    
    work_dir = strcat(out_dir,data_src);
    % Create maps of freezeup/breakup dates
    %create_maps_of_frbr_dates(work_dir);
    
    % Visualize underperforming points within region
    %visualize_region_points(work_dir,"breakup");
    
    % Visualize equally spaced grid points in region
    %visualize_sampled_points(work_dir, "Foxe_Basin", 0);
    
    % Create map of optimal window for NRC freezeup/breakup calc
    create_map_of_optimal_frbr_windows(work_dir);
end

% ------------------------------------------------------------------ %
% --------------------- Analysis functions ------------------------- %
% ------------------------------------------------------------------ %

function [] = create_map_of_optimal_frbr_windows(work_dir, binfilt)
    
    load(work_dir + "out/coords", "coords");
    
    if binfilt
        load(work_dir + "dtvm/DTVM_frbr_dates_binfilt",...
             "fr_days_DTVM", "br_days_DTVM");
        NRC_varied_frbr_dir =... 
        dir(work_dir + "dtvm/NRC_frbr_dates_varied_periods_binfilt/*.mat");

    else
        load(work_dir + "dtvm/DTVM_frbr_dates",...
             "fr_days_DTVM", "br_days_DTVM");
        NRC_varied_frbr_dir =... 
        dir(work_dir + "dtvm/NRC_frbr_dates_varied_periods/*.mat");
    end
    
    num_of_locs = size(coords, 1);                      
                          
    fr_dates_over_all_periods = nan(length(NRC_varied_frbr_dir),num_of_locs);
    br_dates_over_all_periods = nan(length(NRC_varied_frbr_dir),num_of_locs);
    periods_vec = nan(1,length(NRC_varied_frbr_dir));
    
    for ifile = 1:length(NRC_varied_frbr_dir)
        fname = NRC_varied_frbr_dir(ifile).name;
        str_cell = split(fname, "_");
        
        period = str_cell{end}(1:end-4);
        periods_vec(ifile) = str2double(period);
        
        load(NRC_varied_frbr_dir(ifile).folder+"\"+NRC_varied_frbr_dir(ifile).name,...
             "fr_days_NRC", "br_days_NRC");
         
        fr_dates_over_all_periods(ifile,:) = fr_days_NRC;
        br_dates_over_all_periods(ifile,:) = br_days_NRC;
    end
    
    optimal_NRC_br_periods_vec = nan(1,num_of_locs);
    optimal_NRC_fr_periods_vec = nan(1,num_of_locs);
    
    for k = 1:num_of_locs
        fr_diffs = abs(fr_dates_over_all_periods(:,k) - fr_days_DTVM(k));
        br_diffs = abs(br_dates_over_all_periods(:,k) - br_days_DTVM(k));
        
        [~,I_fr] = min(fr_diffs);
        [~,I_br] = min(br_diffs);
        
        optimal_NRC_br_periods_vec(k) = periods_vec(I_br);
        optimal_NRC_fr_periods_vec(k) = periods_vec(I_fr);
    end
    
    save_dir = work_dir+"maps/Binarized_filtered_maps/";
    create_frbr_dates_map(save_dir, coords, optimal_NRC_br_periods_vec,...
                          "Breakup optimal NRC windows", [5 30]);
    create_frbr_dates_map(save_dir, coords, optimal_NRC_fr_periods_vec,...
                          "Freezeup optimal NRC windows", [5 30]);
end

function [] = create_maps_of_frbr_dates(work_dir, binfilt)
    
    load(work_dir+"out/coords", "coords");
    
    if binfilt
        load(work_dir+"dtvm/DTVM_frbr_dates_binfilt", "br_days_DTVM", "fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates_binfilt", "br_days_NRC", "fr_days_NRC");
        save_dir = work_dir + "maps/Binarized_filtered_maps/";
    else
        load(work_dir+"dtvm/DTVM_frbr_dates", "br_days_DTVM", "fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates", "br_days_NRC", "fr_days_NRC");
        save_dir = work_dir + "maps/Raw_maps/";
    end
        
    create_frbr_dates_map(save_dir, coords, fr_days_NRC, "NRC Freeze-up days", [1 365]);
    create_frbr_dates_map(save_dir, coords, br_days_NRC, "NRC Breakup days", [1 365]);
    create_frbr_dates_map(save_dir, coords, fr_days_DTVM, "DTVM Freeze-up days", [1 365]);
    create_frbr_dates_map(save_dir, coords, br_days_DTVM, "DTVM Breakup days", [1 365]);
    
    create_frbr_dates_map(save_dir, coords, br_days_DTVM-br_days_NRC,...
                          "DTVM-NRC Breakup Difference", [-100 100]);
    create_frbr_dates_map(save_dir, coords, fr_days_DTVM-fr_days_NRC,...
                          "DTVM-NRC Freeze-up Difference", [-100 100]);
                      
    disp("Done creating freeze-up/breakup maps");
end

function [] = visualize_sampled_points(work_dir, region_name, histograms, binfilt)
    
    load(work_dir+"out/coords","coords");
    
    if binfilt
        load(work_dir+"out/sic_mats_binarized_filtered",...
                      "sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates_binfilt","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates_binfilt","br_days_NRC","fr_days_NRC");
        load(work_dir+"dtvm/DTVM_frbr_indexes_binfilt","BR_index","FR_index");
    else
        load(work_dir+"out/sic_mats","sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates","br_days_NRC","fr_days_NRC");
        load(work_dir+"dtvm/DTVM_frbr_indexes","BR_index","FR_index");
    end
    
    [lon_bounds, lat_bounds] = get_region_bounds(region_name);
    
    sampled_pts = sample_points_from_grid(lon_bounds, lat_bounds, 5);
    
    % Find closest coords and create visualization for each point
    for k = 1:length(sampled_pts)
        
        pt = sampled_pts(k,:);
        % Use the index of the sampled point
        indx = find_closest_coords_index(coords, pt);
        
        location = coords(indx,:);
        lon = num2str(location(1));
        lat = num2str(location(2));
        
        if histograms
            % Create histograms instead of time series visualization
            potential_br_dates = BR_index(indx,:);
            potential_fr_dates = FR_index(indx,:);
            
            br_date = br_days_DTVM(indx);
            fr_date = fr_days_DTVM(indx);
            
            plot_histogram(work_dir, potential_br_dates, lon, lat,...
                           br_date, "Breakup");
            plot_histogram(work_dir, potential_fr_dates, lon, lat,...
                           fr_date, "Freezeup");
            
        else
            % Extracting timeseries, freezeup/breakup dates, coord using
            % index
            location = coords(indx,:);
            sic_ts = sic_mat(indx,:);
            sic_mean_ts = sic_mean_mat(indx,:);
            sic_std_ts = sic_std_mat(indx,:);

            frbr_days = [br_days_NRC(indx),fr_days_NRC(indx),...
                         br_days_DTVM(indx),fr_days_DTVM(indx)];
            
            save_dir = work_dir+"points/"+region_name+"_binfilt/";
            if not(isfolder(save_dir, "dir"))
                mkdir(save_dir);
            end
            
            savename = save_dir+"visualization_at_"+lon+"_"+lat+".png";

            % Visualize timeseries and freezeup/breakup dates at chosen
            % location
            visualize_location(sic_ts, sic_mean_ts, sic_std_ts, frbr_days,...
                               location, savename);
        end
    end
end

% Not recommended for use
function [] = visualize_region_points(work_dir, day_type, region_name, binfilt)
    
    load(work_dir+"out/coords","coords");
    
    if binfilt
        load(work_dir+"out/sic_mats_binarized_filtered",...
                      "sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates_binfilt","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates_binfilt","br_days_NRC","fr_days_NRC");
    else
        load(work_dir+"out/sic_mats","sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates","br_days_NRC","fr_days_NRC");
    end
    
    if strcmp(day_type,"freezeup")
        diffs=fr_days_DTVM-fr_days_NRC;
    elseif strcmp(day_type,"breakup")
        diffs=br_days_DTVM-br_days_NRC;
    end
    
    num_of_locations = length(coords);
    
    [lon_bounds, lat_bounds] = get_region_bounds(region_name);
    
    % Controlling parameters
    created = 0;
    creation_limit = 3;
    day_difference_cutoff = 50;
    
    % Randomize selection of coordinates within region
    s = RandStream("mlfg6331_64");
    inds = randsample(s,num_of_locations,num_of_locations);
    
    % Now loop over potential locations
    for k = 1:num_of_locations
        ind = inds(k);
        coord = coords(ind,:);
        % Determine if location falls within bounds of region
        if coord(1) > lon_bounds(1) && coord(1) < lon_bounds(2)...
        && coord(2) > lat_bounds(1) && coord(2) < lat_bounds(2)
            day_diff = diffs(ind);
            if ~isnan(day_diff)
                if abs(day_diff) > day_difference_cutoff
                    % Time-series + map
                    sic_ts = sic_mat(ind,:);
                    sic_mean_ts = sic_mean_mat(ind,:);
                    sic_std_ts = sic_std_mat(ind,:);
                    frbr_days = [br_days_NRC(ind),fr_days_NRC(ind),...
                                 br_days_DTVM(ind),fr_days_DTVM(ind)];
                    lon = num2str(coord(1));
                    lat = num2str(coord(2));
                    
                    savename = work_dir+"points/visualization_at_"+lon+"_"+lat+".png";
                
                    visualize_location(sic_ts, sic_mean_ts, sic_std_ts,...
                                       frbr_days, coord, savename);
                    % Histogram
                    
                    created = created + 1;
                    if created >= creation_limit
                        break
                    end
                end
            end
        end
    end
end

% ------------------------------------------------------------------ %
% --------------------- Plotting functions ------------------------- %
% ------------------------------------------------------------------ %

function [] = visualize_location(sic_ts, sic_mean_ts, sic_std_ts, frbr_days, location, savename)
    
    lon = location(1);
    lat = location(2);
    
    figure("visible","off");
                       
    % SIC  
    sic_ax = subplot(3,2,2);
    hold(sic_ax,"on");
    plot(sic_ax,1:365,sic_ts,"Color","b","DisplayName","Signal");
    plot(sic_ax,[frbr_days(1),frbr_days(1)],[0 1],"Color","r","DisplayName","NRC Breakup");
    plot(sic_ax,[frbr_days(2),frbr_days(2)],[0 1],"Color","g","DisplayName","NRC Freezeup");
    plot(sic_ax,[frbr_days(3),frbr_days(3)],[0 1],"Color","c","DisplayName","DTVM Breakup");
    plot(sic_ax,[frbr_days(4),frbr_days(4)],[0 1],"Color","m","DisplayName","DTVM Freezeup");
    xlabel(sic_ax,"Days of year");
    ylabel(sic_ax,"SIC");
    xlim(sic_ax, [0 365]);
    
    % SIC mean 
    sic_mean_ax = subplot(3,2,4);
    hold(sic_mean_ax,"on");
    plot(sic_mean_ax,1:365,sic_mean_ts,"Color","b","DisplayName","Signal");
    plot(sic_mean_ax,[frbr_days(1),frbr_days(1)],[0 1],"Color","r","DisplayName","NRC Breakup");
    plot(sic_mean_ax,[frbr_days(2),frbr_days(2)],[0 1],"Color","g","DisplayName","NRC Freezeup");
    plot(sic_mean_ax,[frbr_days(3),frbr_days(3)],[0 1],"Color","c","DisplayName","DTVM Breakup");
    plot(sic_mean_ax,[frbr_days(4),frbr_days(4)],[0 1],"Color","m","DisplayName","DTVM Freezeup");
    xlabel(sic_mean_ax,"Days of year");
    ylabel(sic_mean_ax,"SIC Mean");
    xlim(sic_mean_ax, [0 365]);
    
    % SIC std deviation
    sic_std_ax = subplot(3,2,[5 6]);
    hold(sic_std_ax,"on");
    plot(sic_std_ax,1:365,sic_std_ts,"Color","b","DisplayName","Signal");
    plot(sic_std_ax,[frbr_days(1),frbr_days(1)],[0 1],"Color","r","DisplayName","NRC Breakup");
    plot(sic_std_ax,[frbr_days(2),frbr_days(2)],[0 1],"Color","g","DisplayName","NRC Freezeup");
    plot(sic_std_ax,[frbr_days(3),frbr_days(3)],[0 1],"Color","c","DisplayName","DTVM Breakup");
    plot(sic_std_ax,[frbr_days(4),frbr_days(4)],[0 1],"Color","m","DisplayName","DTVM Freezeup");
    xlabel(sic_std_ax,"Days of year");
    ylabel(sic_std_ax,"SIC std deviation");
    xlim(sic_std_ax, [0 365]);
    
    legend(sic_std_ax,"Location","best");

    % Create map
    map_ax = subplot(3,2,[1 3]);

    axesm("mercator","MapLatLimit",[50 70],"MapLonLimit",[-100 -50]);
    setm(map_ax,"mlabellocation",15,"plabellocation",5,"MLabelParallel","south");
    mlabel on;
    plabel on;
    framem on;
    axis off;
    provinces = shaperead("./data/ShapeFiles/PROVINCE.SHP","UseGeoCoords",...
                          true,"BoundingBox",[-100,50;-50,70]);
    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];
    setm(map_ax,"FFaceColor",OceanColor);
    geoshow(provinces,"FaceColor",LandColor);
    geoshow(lat, lon, "DisplayType","Point","Marker","x","Color","red");
    
    % Title the whole visualization
    lon = num2str(lon);
    lat = num2str(lat);
    title(["SIC at (",lon,",",lat,")"]);
    
    saveas(map_ax, savename);
    clf;
end

function [] = create_frbr_dates_map(save_dir, coords, dates, plot_title, color_lims)
    
    figure("visible","off");

    lons = coords(:,1);
    lats = coords(:,2);

    scatter(lons,lats,10,dates,"filled");
    colorbar;
    caxis(color_lims);
        
    title(plot_title);
    xlabel("Longitude");
    ylabel("Latitude");
    save_fname = save_dir + plot_title + ".png";
    
    disp("Saving freezeup/breakup map at path: " + save_fname);
    saveas(gca, save_fname);
    
    clf;
end

function [] = plot_histogram(save_dir, dates_index, lon, lat, flagged_date, day_type_str)

    figure("visible","off");

    histogram(dates_index, 1:365, "Normalization", "probability");
    hold on;
    plot(gca,[flagged_date flagged_date],[0 1]);
    ylim([0 1]);

    fig_title = day_type_str+" histogram"+" at"+" ("+lon+","+lat+")";
    title(fig_title);

    xlabel("Day of year");
    ylabel("Normalized frequency");

    padding=10;
    if ~isnan(dates_index)
        xlim([min(dates_index)-padding max(dates_index)+padding]);
    end
    
    save_fname = save_dir+"histograms/"+day_type_str+"_histogram_at_"+lon+"_"+lat+".png";
    saveas(gca, save_fname);
    
    clf;
end

% ------------------------------------------------------------------ %
% --------------------- Helper functions --------------------------- %
% ------------------------------------------------------------------ %

function [indx] = find_closest_coords_index(coords, chosen_coord)
    
    diffs_norms = vecnorm(chosen_coord - coords, 2, 2);
    [~,indx] = min(diffs_norms,[],"all","linear");
end

function [lon_bounds, lat_bounds] = get_region_bounds(name)
    
    region_names = ["James_Bay", "Hudson_Bay", "Foxe_Basin",...
                    "Hudson_Strait", "Baffin_Sea", "East_Coast"];
                
              %Lat_bounds Lon_bounds
    regions = {{[51 55], [-83 -78]},... % James Bay
               {[55 64], [-96 -77]},... % Hudson Bay
               {[64 70], [-85 -75]},... % Foxe Basin
               {[58 65], [-77 -66]},... % Hudson Strait
               {[63 70], [-65 -50]},... % Baffin Sea
               {[50 62], [-64 -50]}};   % East Coast
    
    regions_dict = containers.Map(region_names, regions);

    select_region = regions_dict(name);
    lat_bounds = select_region(1);
    lon_bounds = select_region(2);
end

function [points] = sample_points_from_grid(lon_bounds, lat_bounds, grid_res)
    
    % Sample lons and lats from bounds
    lon_sample = linspace(lon_bounds(1),lon_bounds(2),grid_res);
    lat_sample = linspace(lat_bounds(1),lat_bounds(2),grid_res);
    points = cell(grid_res, grid_res);
    
    % Combine sampled lons/lats for new coordinates
    for m = 1:grid_res
        for n = 1:grid_res
            lon = lon_sample(m);
            lat = lat_sample(n);
            points{m,n} = [lon lat];
        end
    end
    
    % Flatten to 1D array
    points = vertcat(points{:});
end
