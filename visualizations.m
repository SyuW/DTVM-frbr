% ------------------------------------------------------------------ %
% --------------------- Main function ------------------------------ %
% ------------------------------------------------------------------ %

% Call the execution
visualization_main_exec();

function [] = visualization_main_exec()
    out_directory = './dtvm_outputs/';
    
    % Create maps of freezeup/breakup dates
    % create_maps_of_frbr_dates(out_directory);
    
    % Visualize underperforming points within region
    % visualize_region_points(out_directory,'breakup');
    
    % Visualize equally spaced grid points in region
    visualize_sampled_points(out_directory);
end

% ------------------------------------------------------------------ %
% --------------------- Analysis functions ------------------------- %
% ------------------------------------------------------------------ %

function [] = create_maps_of_frbr_dates(out_dir)

    load(strcat(out_dir,'dtvm/','DTVM_frbr_dates'),...
         'br_days_DTVM', 'fr_days_DTVM');
    load(strcat(out_dir,'dtvm/','NRC_frbr_dates'),...
         'br_days_NRC', 'fr_days_NRC');
     
    create_frbr_dates_map(out_dir, fr_days_NRC, 'NRC Freeze-up days', [1 365]);
    create_frbr_dates_map(out_dir, br_days_NRC, 'NRC Breakup days', [1 365]);
    create_frbr_dates_map(out_dir, fr_days_DTVM, 'DTVM Freeze-up days', [1 365]);
    create_frbr_dates_map(out_dir, br_days_DTVM, 'DTVM Breakup days', [1 365]);
    
    create_frbr_dates_map(out_dir, br_days_DTVM-br_days_NRC,...
                          'DTVM-NRC Breakup Difference', [-100 100]);
    create_frbr_dates_map(out_dir, fr_days_DTVM-fr_days_NRC,...
                          'DTVM-NRC Freeze-up Difference', [-100 100]);
                      
    disp('Done creating freeze-up/breakup maps');
end

function [] = visualize_sampled_points(out_dir)

    load(strcat(out_dir,'dtvm/','DTVM_frbr_dates'),'br_days_DTVM','fr_days_DTVM');
    load(strcat(out_dir,'dtvm/','NRC_frbr_dates'),'br_days_NRC','fr_days_NRC');
    load(strcat(out_dir,'out','sic_mats'),'sic_mat','sic_mean_mat','sic_std_mat');
    
    load(strcat(out_dir,'out/','coords'),'coords');
    
    [lon_bounds, lat_bounds] = get_region_bounds('Foxe Basin');
    
    sampled_pts = sample_points_from_grid(lon_bounds, lat_bounds, 5);
    
    % Find closest coords and create visualization for each point
    for k = 1:length(sampled_pts)
        
        pt = sampled_pts(k,:);
        indx = find_closest_coords_index(coords, pt);
        
        % Extracting timeseries, freezeup/breakup dates, coord using index
        location = coords(indx,:);
        sic_ts = sic_mat(indx,:);
        sic_mean_ts = sic_mean_mat(indx,:);
        sic_std_ts = sic_std_mat(indx,:);
        
        frbr_days = [br_days_NRC(indx),fr_days_NRC(indx),...
                     br_days_DTVM(indx),fr_days_DTVM(indx)];
        
        lon = num2str(location(1));
        lat = num2str(location(2));
        savename = strcat(out_dir,...
                   'points/visualization_at_',lon,'_',lat,'.png');
                             
        visualize_location(sic_ts, sic_mean_ts, sic_std_ts, frbr_days,...
                           location, savename);
    end
end

function [] = visualize_region_points(out_dir, day_type)
    
    load(strcat(out_dir,'dtvm/','DTVM_frbr_dates'),'br_days_DTVM','fr_days_DTVM');
    load(strcat(out_dir,'dtvm/','NRC_frbr_dates'),'br_days_NRC','fr_days_NRC');
    
    if strcmp(day_type,'freezeup')
        diffs=fr_days_DTVM-fr_days_NRC;
    elseif strcmp(day_type,'breakup')
        diffs=br_days_DTVM-br_days_NRC;
    end
    
    load(strcat(out_dir,'out/','coords'),'coords');
    num_of_locations = length(coords);
    
    load(strcat(out_dir,'out/','sic_mats_binarized_filtered'),...
                'sic_mat','sic_mean_mat','sic_std_mat');
    
    [lon_bounds, lat_bounds] = get_region_bounds('Hudson Bay');
    
    % Controlling parameters
    created = 0;
    creation_limit = 3;
    day_difference_cutoff = 50;
    
    % Randomize selection of coordinates within region
    s = RandStream('mlfg6331_64');
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
                    
                    savename = strcat(out_dir,...
                    'points/visualization_at_',lon,'_',lat,'.png');
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
    
    lon = location(1);
    lat = location(2);
    
    figure('visible','off');
                       
    % SIC  
    sic_ax = subplot(3,2,2);
    hold(sic_ax,'on');
    plot(sic_ax,1:365,sic_ts,'Color','b','DisplayName','Signal');
    plot(sic_ax,[frbr_days(1),frbr_days(1)],[0 1],'Color','r','DisplayName','NRC Breakup');
    plot(sic_ax,[frbr_days(2),frbr_days(2)],[0 1],'Color','g','DisplayName','NRC Freezeup');
    plot(sic_ax,[frbr_days(3),frbr_days(3)],[0 1],'Color','c','DisplayName','DTVM Breakup');
    plot(sic_ax,[frbr_days(4),frbr_days(4)],[0 1],'Color','m','DisplayName','DTVM Freezeup');
    xlabel(sic_ax,'Days of year');
    ylabel(sic_ax,'SIC');
    xlim(sic_ax, [0 365]);
    
    % SIC mean 
    sic_mean_ax = subplot(3,2,4);
    hold(sic_mean_ax,'on');
    plot(sic_mean_ax,1:365,sic_mean_ts,'Color','b','DisplayName','Signal');
    plot(sic_mean_ax,[frbr_days(1),frbr_days(1)],[0 1],'Color','r','DisplayName','NRC Breakup');
    plot(sic_mean_ax,[frbr_days(2),frbr_days(2)],[0 1],'Color','g','DisplayName','NRC Freezeup');
    plot(sic_mean_ax,[frbr_days(3),frbr_days(3)],[0 1],'Color','c','DisplayName','DTVM Breakup');
    plot(sic_mean_ax,[frbr_days(4),frbr_days(4)],[0 1],'Color','m','DisplayName','DTVM Freezeup');
    xlabel(sic_mean_ax,'Days of year');
    ylabel(sic_mean_ax,'SIC Mean');
    xlim(sic_mean_ax, [0 365]);
    
    % SIC std deviation
    sic_std_ax = subplot(3,2,[5 6]);
    hold(sic_std_ax,'on');
    plot(sic_std_ax,1:365,sic_std_ts,'Color','b','DisplayName','Signal');
    plot(sic_std_ax,[frbr_days(1),frbr_days(1)],[0 1],'Color','r','DisplayName','NRC Breakup');
    plot(sic_std_ax,[frbr_days(2),frbr_days(2)],[0 1],'Color','g','DisplayName','NRC Freezeup');
    plot(sic_std_ax,[frbr_days(3),frbr_days(3)],[0 1],'Color','c','DisplayName','DTVM Breakup');
    plot(sic_std_ax,[frbr_days(4),frbr_days(4)],[0 1],'Color','m','DisplayName','DTVM Freezeup');
    xlabel(sic_std_ax,'Days of year');
    ylabel(sic_std_ax,'SIC std deviation');
    xlim(sic_std_ax, [0 365]);
    
    legend(sic_std_ax,'Location','best');

    % Create map
    map_ax = subplot(3,2,[1 3]);

    axesm('mercator','MapLatLimit',[50 70],'MapLonLimit',[-100 -50]);
    setm(map_ax,'mlabellocation',15,'plabellocation',5,'MLabelParallel','south');
    mlabel on;
    plabel on;
    framem on;
    axis off;
    provinces = shaperead('./data/ShapeFiles/PROVINCE.SHP','UseGeoCoords',...
                          true,'BoundingBox',[-100,50;-50,70]);
    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];
    setm(map_ax,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);
    geoshow(lat, lon, 'DisplayType','Point','Marker','x','Color','red');
    
    % Title the whole visualization
    lon = num2str(lon);
    lat = num2str(lat);
    title(['SIC at (',lon,',',lat,')']);
    
    saveas(map_ax, savename);
    clf;
end

function [] = create_frbr_dates_map(out_dir, dates, plot_title, color_lims)

    load(strcat(out_dir,'out/','coords'),'coords');
    
    figure('visible','off');

    lons = coords(:,1);
    lats = coords(:,2);

    scatter(lons,lats,10,dates,'filled');
    colorbar;
    caxis(color_lims);
        
    title(plot_title);
    xlabel('Longitude');
    ylabel('Latitude');
    save_fname = strcat(out_dir,'maps/',plot_title,'.png');
    saveas(gca, save_fname);
    
    clf;
end

function [] = plot_histogram(dates_index, location, flagged_date, day_type_str)

    lon = num2str(location(1));
    lat = num2str(location(2));
    
    figure('visible','off');

    histogram(dates_index, 1:365, 'Normalization', 'probability');
    hold on;
    plot(gca,[flagged_date flagged_date],[0 1]);
    ylim([0 1]);

    fig_title = strcat(day_type_str,' histogram',' at',' (',lon,',',lat,')');
    title(fig_title);
    save_fname = strcat(output_directory,'histograms/',day_type_str,...
                        '_histogram','_at',lon,'_',lat,'.png');

    xlabel('Day of year');
    ylabel('Normalized frequency');

    padding=10;
    if ~isnan(dates_index)
        xlim([min(dates_index)-padding max(dates_index)+padding]);
    end
    saveas(gca, save_fname);
    
    clf;
end

% ------------------------------------------------------------------ %
% --------------------- Helper functions --------------------------- %
% ------------------------------------------------------------------ %

function [indx] = find_closest_coords_index(coords, chosen_coord)
    
    diffs_norms = vecnorm(chosen_coord - coords, 2, 2);
    [~,indx] = min(diffs_norms,[],'all','linear');
end

function [lon_bounds, lat_bounds] = get_region_bounds(name)
    
    region_names = {'James Bay', 'Hudson Bay', 'Foxe Basin',...
                    'Hudson Strait', 'Baffin Sea', 'East Coast'};
                
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
