% ------------------------------------------------------------------ %
% --------------------- Main function ------------------------------ %
% ------------------------------------------------------------------ %

% Call the execution
tic;
visualization_main_exec("2007_esacci/");
toc;

function [] = visualization_main_exec(data_src, process_type)
    % Entry point of execution of DTVM method
    %
    % arguments:
    %   data_src - string describing which data source to use
    %       allowed: (2007_esacci/, 2008_esacci/, 
    %                 2009_esacci/, 2010_esacci/)
    %
    % return: None
    %
    % loaded variables:
    %   sic_mat - 2D matrix of sea ice concentrations (SIC)
    %   sic_std_mat - 2D matrix of moving std deviation of SIC
    %
    % saved variables: None
    
    out_dir = "./out/";
    work_dir = out_dir+data_src;
    
    % Create maps of freezeup/breakup dates
    create_maps_of_frbr_dates(work_dir, process_type);
    
    %create_landmask(work_dir);
    % Visualize underperforming points within region
    %visualize_region_points(work_dir,"breakup");
    
    % Visualize equally spaced grid points in region
    %visualize_sampled_points(work_dir, "Foxe_Basin", 0, 1);
    
    % Create map of optimal window for NRC freezeup/breakup calc
    %create_map_of_optimal_frbr_windows(work_dir);
end

% ------------------------------------------------------------------ %
% --------------------- Analysis functions ------------------------- %
% ------------------------------------------------------------------ %

function [] = create_map_of_optimal_frbr_windows(work_dir, binfilt)
    % Create maps of NRC window values that yield NRC freeze-up/breakup
    % dates closest to DTVM calculated dates
    %
    % arguments:
    %   work_dir -  
    %   binfilt -
    %
    % return: None
    %
    % loaded variables:
    %   DTVM_frbr_dates(_binfilt)
    %       fr_days_DTVM - vector of DTVM freeze-up dates
    %       br_days_DTVM - vector of DTVM breakup dates
    %   NRC_frbr_p_{window}
    %       fr_days_NRC - vector of NRC freeze-up dates for {window}
    %       br_days_NRC - vector of NRC breakup dates for {window}
    %
    % saved variables: None
    
    load(work_dir + "mats/coords", "coords");
    
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
    % Create maps of freeze-up/breakup dates
    %
    % arguments:
    %   work_dir - path to directory with dtvm outputs
    %   binfilt - boolean for whether processed/raw data is used
    %
    % return: None
    %
    % loaded variables:
    %   coords
    %       coords - 2D matrix of coordinates for every location
    %   DTVM_frbr_dates(_binfilt)
    %       br_days_DTVM - vector of DTVM breakup dates
    %       fr_days_DTVM - vector of DTVM freeze-up dates
    %   NRC_frbr_dates(_binfilt)
    %       br_days_NRC - vector of NRC breakup dates
    %       fr_days_NRC - vector of NRC freezeup dates
    %
    % saved variables: None
    
    load(work_dir+"mats/coords", "coords");
    
    if binfilt
        load(work_dir+"dtvm/DTVM_frbr_dates_binfilt", "br_days_DTVM", "fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates_binfilt", "br_days_NRC", "fr_days_NRC");
        save_dir = work_dir + "maps/Binarized_filtered_maps/";
    else
        load(work_dir+"dtvm/DTVM_frbr_dates", "br_days_DTVM", "fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates", "br_days_NRC", "fr_days_NRC");
        save_dir = work_dir + "maps/Raw_maps/";
    end
    
    % make the save directory if it doesn't exist
    if not(isfolder(save_dir))
        mkdir(save_dir);
    end
    
    pad = 2;
        
    create_frbr_dates_map(save_dir, coords, fr_days_NRC, "NRC Freeze-up days",...
                          [min(fr_days_NRC)-pad, max(fr_days_NRC)+pad]);
                      
    create_frbr_dates_map(save_dir, coords, br_days_NRC, "NRC Breakup days",...
                          [min(br_days_NRC)-pad, max(br_days_NRC)+pad]);
                      
    create_frbr_dates_map(save_dir, coords, fr_days_DTVM, "DTVM Freeze-up days",...
                          [min(fr_days_DTVM)-pad, max(fr_days_DTVM)+pad]);
                      
    create_frbr_dates_map(save_dir, coords, br_days_DTVM, "DTVM Breakup days",...
                          [min(br_days_DTVM)-pad, max(br_days_DTVM)+pad]);
    
    create_frbr_dates_map(save_dir, coords, br_days_DTVM-br_days_NRC,...
                          "DTVM-NRC Breakup Difference", [-30 30]);
                      
    create_frbr_dates_map(save_dir, coords, fr_days_DTVM-fr_days_NRC,...
                          "DTVM-NRC Freeze-up Difference", [-30 30]);
                      
    disp("Done creating freeze-up/breakup maps");
end

function [] = visualize_sampled_points(work_dir, region_name, histograms, binfilt)
    % Create visualizations for regularly spaced points in region
    %
    % arguments:
    %   work_dir - working directory path
    %   region_name - region name string
    %   histograms - boolean for whether histograms should be gen.
    %   binfilt - whether to use raw or processed data
    %
    % return: None
    %
    % loaded variables:
    %   sic_mats(_binarized_filtered)
    %       sic_mat
    %       sic_mean_mat
    %       sic_std_mat
    %   NRC_frbr_dates(_binfilt)
    %       br_days_NRC
    %       fr_days_NRC
    %   NRC_frbr_indexes(_binfilt)
    %       BR_index
    %       FR_index
    %
    % saved variables: None
    
    load(work_dir+"mats/coords","coords");
    
    [lon_bounds, lat_bounds] = get_region_bounds(region_name);
    
    sampled_pts = sample_points_from_grid(lon_bounds, lat_bounds, 5);
    
    if binfilt
        load(work_dir+"mats/sic_mats_binarized_filtered",...
                      "sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates_binfilt","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates_binfilt","br_days_NRC","fr_days_NRC");
        load(work_dir+"dtvm/DTVM_frbr_indexes_binfilt","BR_index","FR_index");
        save_dir = work_dir+"points/"+region_name+"_binfilt/";
    else
        load(work_dir+"mats/sic_mats","sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates","br_days_NRC","fr_days_NRC");
        load(work_dir+"dtvm/DTVM_frbr_indexes","BR_index","FR_index");
        save_dir = work_dir+"points/"+region_name+"/";
    end
    
    if not(isfolder(save_dir))
        mkdir(save_dir);
    end
    
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
            
            savename = save_dir+"visualization_at_"+lon+"_"+lat+".png";

            % Visualize timeseries and freezeup/breakup dates at chosen
            % location
            visualize_location(sic_ts, sic_mean_ts, sic_std_ts, frbr_days,...
                               location, savename);
        end
    end
end

function [] = visualize_region_points(work_dir, day_type, region_name, binfilt)
    % Not recommended for use
    
    load(work_dir+"mats/coords","coords");
    
    if binfilt
        load(work_dir+"mats/sic_mats_binarized_filtered",...
                      "sic_mat","sic_mean_mat","sic_std_mat");
        load(work_dir+"dtvm/DTVM_frbr_dates_binfilt","br_days_DTVM","fr_days_DTVM");
        load(work_dir+"dtvm/NRC_frbr_dates_binfilt","br_days_NRC","fr_days_NRC");
    else
        load(work_dir+"mats/sic_mats","sic_mat","sic_mean_mat","sic_std_mat");
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
    % Visualize time-series and freeze-up dates at a location
    %
    % arguments:
    %   sic_ts - vector describing time series of SIC values
    %   sic_mean_ts - vector desribing time series of mean SIC values
    %   sic_std_ts - vector describing time series of std SIC values
    %   frbr_days - vector of DTVM/NRC freezeup/breakup days at loc.
    %   location - vector describing the location
    %   savename - save name of the plot file
    %
    % return: None
    % saved variables: None
    % loaded variables: None
    
    lon = location(1);
    lat = location(2);
    
    figure("visible","off");
                       
    % SIC  
    sic_ax = subplot(3,2,2);
    hold(sic_ax,"on");
    plot(sic_ax,1:length(sic_ts),sic_ts,"Color","b","DisplayName","Signal");
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
    plot(sic_mean_ax,1:length(sic_mean_ts),sic_mean_ts,"Color","b","DisplayName","Signal");
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
    plot(sic_std_ax,1:length(sic_std_ts),sic_std_ts,"Color","b","DisplayName","Signal");
    plot(sic_std_ax,[frbr_days(1),frbr_days(1)],[0 1],"Color","r","DisplayName","NRC Breakup");
    plot(sic_std_ax,[frbr_days(2),frbr_days(2)],[0 1],"Color","g","DisplayName","NRC Freezeup");
    plot(sic_std_ax,[frbr_days(3),frbr_days(3)],[0 1],"Color","c","DisplayName","DTVM Breakup");
    plot(sic_std_ax,[frbr_days(4),frbr_days(4)],[0 1],"Color","m","DisplayName","DTVM Freezeup");
    xlabel(sic_std_ax,"Days of year");
    ylabel(sic_std_ax,"SIC std deviation");
    xlim(sic_std_ax, [0 365]);
    
    legend(sic_std_ax,"Location","best");

    % Create map of region
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
    title("SIC at ("+lon+","+lat+")");
    
    saveas(map_ax, savename);
    clf;
end

function [] = create_hist_and_plot(sic_std, dates_index, location, savename)
    % create histogram + plot of variability signal
    % 
    % arguments:
    %   sic_std - 
    %   dates -
    %   coords -
    %   savename -
    %
    % return: None
    % saved variables: None
    % loaded variables: None
    
    figure("visible","off");
end

function [] = create_frbr_dates_map(save_dir, coords, dates, plot_title,...
                                    color_lims)
    % Create map of freezeup or breakup dates
    %
    % arguments:
    %   save_dir - save directory for plot .png
    %   coords - 2D matrix of all coordinates
    %   dates - freeze-up or breakup dates to plot
    %   plot_title - string for title of plot
    %   cmap_name - variable for plot color map
    %   color_lims - vector for limits of colorbar
    %
    % return: None
    % saved variables: None
    % loaded variables: None
    
    figure("visible","off");

    lons = coords(:,1);
    lats = coords(:,2);
    
    % hard coded marker size for now
    if size(coords, 1) > 20000
        marker_size = 3;
    else
        marker_size = 10;
    end
    
    geoscatter(lats,lons,marker_size,dates,"filled");
    geobasemap colorterrain;
    
    colorbar;
    caxis(color_lims);
        
    title(plot_title);
    save_fname = save_dir + plot_title + ".png";
    
    disp("Saving freezeup/breakup map at path: " + save_fname);
    saveas(gca, save_fname);
    
    clf;
end

function [] = plot_histogram(save_dir, dates_index, lon, lat, flagged_date, day_type_str)
    % Plot histogram of DTVM potential freezeup/breakup dates
    %
    % arguments:
    %   save_dir - directory path to save at
    %   dates_index - potential freezeup/breakup dates at location
    %   lon - longitude of location considered
    %   lat - latitude of location considered
    %   flagged_date - special date at which to draw vertical line, usually
    %   the precise freeze-up/breakup date determined by DTVM
    %   day_type_str - string desribing whether freezeup/breakup is
    %   considered
    %
    % return: None
    % saved variables: None
    % loaded variables: None
    
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
    % Find the index of the closest available coordinate
    % arguments:
    %   coords - 2D matrix of coordinates for each location
    %   chosen_coord - provided coordinate tuple
    %
    % return:
    %   indx - index of closest available coordinate
    %
    % saved variables: None
    %
    % loaded variables: None
    
    diffs_norms = vecnorm(chosen_coord - coords, 2, 2);
    [~,indx] = min(diffs_norms,[],"all","linear");
end

function [lon_bounds, lat_bounds] = get_region_bounds(name)
    % Converts region name string to set of longitude/latitude bounds
    % arguments:
    %   name - string representing region name
    %   
    % return:
    %   lon_bounds - vector of longitude bounds of the region
    %   lat_bounds - vector of latitude bounds of the region
    %
    % saved variables: None
    %
    % loaded variables: None
    
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
    lat_bounds = select_region{1};
    lon_bounds = select_region{2};
end

function [points] = sample_points_from_grid(lon_bounds, lat_bounds, grid_res)
    % Sample coordinates from region at regularly spaced points
    %
    % arguments:
    %   lon_bounds - vector of longitude bounds of the region
    %   lat_bounds - vector of latitude bounds of the region
    %   grid_res - number of points wanted between each bound
    %
    % return:
    %   points - 2D matrix of sampled lon,lat coordinates
    %
    % saved variables: None
    %
    % loaded variables: None
    
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

function [] = create_landmask(work_dir)
    % create landmask array
    %
    % arguments:
    %   work_dir - 
    %   grid_spacing - scalar
    %
    % return:
    %   landmask - 2D logical matrix of 0 - ocean and 1 - land
    
    load(work_dir+"mats/coords", "coords");
    
    avg_lat_step = mean(diff(coords(:,1)));
    avg_lon_step = mean(diff(coords(:,2)));
    
    lon_min = min(coords(:,1));
    lon_max = max(coords(:,1));
    lat_min = min(coords(:,2));
    lat_max = max(coords(:,2));
    
    lon_sample = lon_min:0.2:lon_max;
    lat_sample = lat_min:0.2:lat_max;
    
    [lons, lats] = meshgrid(lat_sample, lon_sample);
    land = landmask(lats, lons, 'North and South America');
    
    keyboard;
    
end
