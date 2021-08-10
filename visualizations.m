% ------------------------------------------------------------------ %
% --------------------- Main function ------------------------------ %
% ------------------------------------------------------------------ %

% Call the execution
tic;
visualization_main_exec("2007_esacci/", "raw");
%visualize_batch()
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
    
    % set working directory to generate visualizations inside
    work_dir = './out/'+data_src+process_type+'/';
    
    % Create maps of freezeup/breakup dates
    %create_maps_of_frbr_dates(work_dir);
    
    %create_landmask(work_dir);
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

function [] = create_map_of_optimal_frbr_windows(work_dir)
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
    
    % load coordinates and DTVM frbr data
    load(work_dir + "mats/coords", "coords");
    load(work_dir + "dtvm/DTVM_frbr_dates","fr_days_DTVM", "br_days_DTVM");
    load(work_dir + "dtvm/NRC_frbr_cubes","NRC_br_cube", "NRC_fr_cube");
    
    % find the index corresponding to the minimal difference from DTVM
    [~,br_optimal] = min(abs(bsxfun(@minus, NRC_br_cube, br_days_DTVM)));
    [~,fr_optimal] = min(abs(bsxfun(@minus, NRC_fr_cube, fr_days_DTVM)));
    
    % since window sizes < 5 weren't considered, set those to NaN
    br_optimal(br_optimal < 5) = nan;
    fr_optimal(fr_optimal < 5) = nan;
    
    % create optimal windows maps
    br_optimal_map = create_frbr_dates_map(coords, br_optimal,...
                                           "Breakup optimal NRC windows",...
                                           [5 30]);
    fr_optimal_map = create_frbr_dates_map(coords, fr_optimal,...
                                           "Freezeup optimal NRC windows",...
                                           [5 30]);
    % saving
    save_dir = work_dir+"maps/";
    saveas(br_optimal_map,save_dir+"Breakup_optimal_NRC_windows_map.png");
    saveas(fr_optimal_map,save_dir+"Freezeup_optimal_NRC_windows_map.png");
end

function [] = create_maps_of_frbr_dates(work_dir)
    % Create maps of freeze-up/breakup dates and differences
    %
    % arguments:
    %   work_dir - path to directory with dtvm outputs
    %   process_type - boolean for whether processed/raw data is used
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
    
    % load coordinates
    load(work_dir+"mats/coords", "coords");
    
    % load freeze-up/breakup data
    load(work_dir+"dtvm/DTVM_frbr_dates", "br_days_DTVM", "fr_days_DTVM");
    load(work_dir+"dtvm/NRC_frbr_dates", "br_days_NRC", "fr_days_NRC");
    
    % define and make the save directory if it doesn't exist
    save_dir = work_dir + "maps/";
    if not(isfolder(save_dir))
        mkdir(save_dir);
    end
    
    pad = 2;
    
    % NRC freeze-up dates map
    NRC_fr_dates_map = create_frbr_dates_map(coords, fr_days_NRC,... 
        "NRC Freeze-up dates",[min(fr_days_NRC)-pad, max(fr_days_NRC)+pad]);
    
    % NRC breakup dates map
    NRC_br_dates_map = create_frbr_dates_map(coords, br_days_NRC,... 
        "NRC Breakup dates",[min(br_days_NRC)-pad, max(br_days_NRC)+pad]);
    
    % DTVM freeze-up dates map                   
    DTVM_fr_dates_map = create_frbr_dates_map(coords, fr_days_DTVM,...
        "DTVM Freeze-up days",[min(fr_days_DTVM)-pad, max(fr_days_DTVM)+pad]);
    
    % DTVM breakup dates map
    DTVM_br_dates_map = create_frbr_dates_map(coords, br_days_DTVM,...
        "DTVM Breakup days",[min(br_days_DTVM)-pad, max(br_days_DTVM)+pad]);
    
    % color bar limit for plotting frbr differences
    c_abs_lim = 20;
    
    % Breakup differences map
    br_diffs_map = create_frbr_dates_map(coords, br_days_DTVM-br_days_NRC,...
        "DTVM-NRC Breakup Difference", [-c_abs_lim c_abs_lim]);
    
    % Freeze-up differences map
    fr_diffs_map = create_frbr_dates_map(coords, fr_days_DTVM-fr_days_NRC,...
        "DTVM-NRC Freeze-up Difference", [-c_abs_lim c_abs_lim]);
    
    % Saving
    saveas(NRC_fr_dates_map,save_dir+"NRC_freezeup_dates.png");
    saveas(NRC_br_dates_map,save_dir+"NRC_breakup_dates.png");
    saveas(DTVM_fr_dates_map,save_dir+"DTVM_freezeup_dates.png");
    saveas(DTVM_br_dates_map,save_dir+"DTVM_breakup_dates.png");
    saveas(br_diffs_map,save_dir+"Breakup_differences.png");
    saveas(fr_diffs_map,save_dir+"Freezeup_difference.png");
                 
    disp("Done creating freeze-up/breakup maps");
end

function [] = visualize_sampled_points(work_dir, region_name, histograms)
    % Create visualizations for regularly spaced points in region
    %
    % arguments:
    %   work_dir - working directory path
    %   region_name - region name string
    %   histograms - boolean for whether histograms should be gen.
    %   process_type - process type
    %
    % return: None
    %
    % loaded variables:
    %   sic_mats
    %       sic_mat
    %       sic_mean_mat
    %       sic_std_mat
    %   NRC_frbr_dates
    %       br_days_NRC
    %       fr_days_NRC
    %   NRC_frbr_indexes
    %       BR_index
    %       FR_index
    %
    % saved variables: None
    
    load(work_dir+"mats/coords","coords");
    
    [lon_bounds, lat_bounds] = get_region_bounds(region_name);
    
    sampled_pts = sample_points_from_grid(lon_bounds, lat_bounds, 5);
    
    % load SIC, SIC mean, SIC variability
    load(work_dir+"mats/sic","sic_mat");
    load(work_dir+"mats/sic_std","sic_std_mat");
    load(work_dir+"mats/sic_mean","sic_mean_mat");
    % load DTVM data
    load(work_dir+"dtvm/DTVM_frbr_indexes","BR_index","FR_index");
    load(work_dir+"dtvm/DTVM_frbr_dates","br_days_DTVM","fr_days_DTVM");
    % load NRC calculated frbr
    load(work_dir+"dtvm/NRC_frbr_dates","br_days_NRC","fr_days_NRC");
    
    % define save directory based on region examined
    save_dir = work_dir+"points/"+region_name+"/";
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
        
        % Create histograms instead of time series visualization
        if histograms
            
            % slice frbr dates and SIC variability at location's index
            potential_br_dates = BR_index(indx,:);
            potential_fr_dates = FR_index(indx,:);
            sic_std_at_loc = sic_std_mat(indx,:);
            
            %keyboard;
            
            create_hist_and_plot(save_dir, sic_std_at_loc,...
                                 potential_br_dates, potential_fr_dates,...
                                 location);
        
        % Create time-series visualization
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

function [] = visualize_location(sic_ts, sic_mean_ts, sic_std_ts,...
                                 frbr_days, location, savename)
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
    % create filled area representing threshold of SIC = 0.15
    a = area(sic_ax,1:length(sic_ts),repmat(0.15,1,length(sic_ts)));
    a.DisplayName = "Open water zone";
    a.FaceColor = "c";
    a.FaceAlpha = 0.1;
    a.LineStyle = "none";
    
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

function [] = create_hist_and_plot(save_dir, sic_std, br_dates, fr_dates, location)
    % create histogram + plot of variability signal
    % 
    % arguments:
    %   sic_std - 
    %   dates -
    %   location - lon,lat coordinates of location, for titling
    %   savename -
    %
    % return: None
    % saved variables: None
    % loaded variables: None
    
    fig = figure("visible","off");
    
    br_date_at_loc = quantile(br_dates, 0.75);
    fr_date_at_loc = quantile(fr_dates, 0.25);
    
    % variability plot with potential freezeup/breakup dates
    sic_ax = subplot(2,1,1);
    hold(sic_ax,'on');
    % draw SIC variability line
    std_line = plot(sic_ax,1:length(sic_std),sic_std);
    std_line.Color = "blue";
    std_line.DisplayName = "SIC variability";
    % draw breakup date
    sic_ax_br_day = plot(sic_ax,[br_date_at_loc, br_date_at_loc],[0, 1.1]);
    sic_ax_br_day.Color = "red";
    sic_ax_br_day.DisplayName = "DTVM Breakup date";
    % draw freeze-up date
    sic_ax_fr_day = plot(sic_ax,[fr_date_at_loc, fr_date_at_loc],[0, 1.1]);
    sic_ax_fr_day.Color = "magenta";
    sic_ax_fr_day.DisplayName = "DTVM Freeze-up date";
    % titling, labels, limits
    sic_title = "SIC variability signal with DTVM flagged freeze-up/breakup";
    title(sic_ax, sic_title);
    xlabel(sic_ax, "Day of year");
    ylabel(sic_ax, "SIC variability");
    xlim(sic_ax, [1 365]);
    ylim(sic_ax, [0 1.1]);
    % draw the legend
    legend(sic_ax,"Location","best");
    
    % histogram
    hist_ax = subplot(2,1,2);
    hold(hist_ax, "on");
    % draw breakup date on histogram
    hist_ax_br_day = plot(hist_ax,[br_date_at_loc, br_date_at_loc],[0, 1.1]);
    hist_ax_br_day.Color = "red";
    hist_ax_br_day.DisplayName = "DTVM Breakup date";
    % draw freezeup date on histogram
    hist_ax_fr_day = plot(hist_ax,[fr_date_at_loc, fr_date_at_loc],[0, 1.1]);
    hist_ax_fr_day.Color = "magenta";
    hist_ax_fr_day.DisplayName = "DTVM Freeze-up date";
    % plot breakup frequency histogram
    hist_ax_br_dates = histogram(hist_ax,br_dates,1:length(br_dates));
    hist_ax_br_dates.DisplayName = "DTVM potential breakup dates";
    hist_ax_br_dates.Normalization = "Probability";
    % plot freezeup frequency histogram
    hist_ax_fr_dates = histogram(hist_ax,fr_dates,1:length(fr_dates));
    hist_ax_fr_dates.DisplayName = "DTVM potential freezeup dates";
    hist_ax_fr_dates.Normalization = "Probability";
    % titling, labels, limits
    hist_title = "Normalized frequency of DTVM potential freeze-up/breakup dates";
    title(hist_ax, hist_title);
    xlabel(hist_ax, "Day of year");
    ylabel(hist_ax, "Normalized frequency");
    xlim(hist_ax, [1 365]);
    ylim(hist_ax, [0 1.1]);
    % draw the legend
    legend(hist_ax,"Location","best");
    
    % finally, save the visualization
    location_str = strjoin(string(location),"_");
    savename = save_dir + "histogram_at_" + location_str + ".png";
    saveas(fig, savename);
   
    clf;
end

function [fig] = create_frbr_dates_map(coords, dates, plot_title, color_lims)
    % Create map of freezeup or breakup dates
    %
    % arguments:
    %   coords - 2D matrix of all coordinates
    %   dates - freeze-up or breakup dates to plot
    %   plot_title - string for title of plot
    %   cmap_name - variable for plot color map
    %   color_lims - vector for limits of colorbar
    %
    % return: fig - figure object of map
    % saved variables: None
    % loaded variables: None
    
    fig = figure("visible","off");

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
end

% ------------------------------------------------------------------ %
% --------------------- Helper functions --------------------------- %
% ------------------------------------------------------------------ %

function [] = visualize_batch()
    % Utility function for visualizing multiple sources of data
    % and multiple processing methods
    %
    % arguments: None
    
    data_srcs = ["2007_esacci/","2008_esacci/"];
    process_types = ["hysteresis","binfilt","raw"];

    for j = 1:length(data_srcs)
        for k = 1:length(process_types)        
            data_src = data_srcs(j);
            process_type = process_types(k);
            fprintf("Source: %s, Process type: %s\n", data_src, process_type);
            visualization_main_exec(data_src, process_type);
            fprintf("\n");
        end
    end
end

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
