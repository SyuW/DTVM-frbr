% ------------------------------------------------------------------ %
% --------------------- Main function ------------------------------ %
% ------------------------------------------------------------------ %

% Call the execution
visualization_main_exec();

function [] = visualization_main_exec()
    out_directory = './dtvm_outputs/';
    create_maps_of_frbr_dates(out_directory);
end

% ------------------------------------------------------------------ %
% --------------------- Analysis functions ------------------------- %
% ------------------------------------------------------------------ %

function [] = create_maps_of_frbr_dates(out_dir)

    load(strcat(out_dir,'out/','DTVM_frbr_dates'),...
         'br_days_DTVM', 'fr_days_DTVM');
    load(strcat(out_dir,'out/','NRC_frbr_dates'),...
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

function [] = plot_timeseries_at_random_points()
    % Random coord sampling
    num_of_coords = 5;
    inds = randsample(size(coords,1), num_of_coords);
    for k = 1:num_of_coords
        ind = inds(k);
        coord = coords(ind,:);
        choices{k} = coord;
    end
    % Plot for coordinate choices
    for k = 1:length(choices)
        choice = choices{k};
        chosen_index = coord_to_closest_coords_index(coords, choice);
        chosen_coord = coords(chosen_index,:);
        visualize_location({sic_mat sic_mean_mat sic_std_mat},...
                        {'SIC' 'Mean SIC' 'SIC Std deviation'},... 
                        chosen_index,...
                        chosen_coord,...
                        date_vec,...
                        {freezeup_days_NRC,breakup_days_NRC,freezeup_days_DTVM,breakup_days_DTVM},...
                        0);
    end
    disp('Done plotting time-series for random points');
end

function [] = visualize_specific_points_in_region(dtvm_dir, out_dir)

    load(strcat(dtvm_dir, 'DTVM_frbr_dates'),'br_days_DTVM','fr_days_DTVM');
    load(strcat(dtvm_dir, 'NRC_frbr_dates'),'br_days_NRC','fr_days_NRC');
    load(strcat(out_dir, 'coords'),'coords');
    
    num_of_locations = length(coords);
    
    [lon_bounds, lat_bounds] = get_region('Hudson Bay');
    
    % Controlling parameters
    is_breakup = 1;
    created = 0;
    creation_limit = 3;
    day_difference_cutoff = 80;

    s = RandStream('mlfg6331_64');
    inds = randsample(s,num_of_locations,num_of_locations);
    
    for k = 1:num_of_locations
        ind = inds(k);
        coord = coords(ind,:);
        % Determine if location falls within bounds of region of interest
        if coord(1) > lon_bounds(1) && coord(1) < lon_bounds(2)...
        && coord(2) > lat_bounds(1) && coord(2) < lat_bounds(2)
            if is_breakup
                dtvm_day = breakup_days_DTVM(k);
                nrc_day = breakup_days_NRC(k);
            else
                dtvm_day = freezeup_days_DTVM(k);
                nrc_day = freezeup_days_NRC(k);
            end
            day_diff = dtvm_day - nrc_day;
            if ~isnan(day_diff)
                if abs(day_diff) > day_difference_cutoff
                    % Time-series + map
                    visualize_location({sic_mat sic_mean_mat sic_std_mat},...
                                    {'SIC' 'Mean SIC' 'SIC Std deviation'},...
                                    k,...
                                    coord,...
                                    date_vec,...
                                    {freezeup_days_NRC,breakup_days_NRC,freezeup_days_DTVM,breakup_days_DTVM},...
                                    1);
                    % Histogram
                    if is_breakup
                        days_for_histogram = BR_index(k,:);
                    else
                        days_for_histogram = FR_index(k,:);
                    end
                    plot_histogram(days_for_histogram, coord, is_breakup);
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
    provinces = shaperead('canada_shape_files/PROVINCE.SHP','UseGeoCoords',...
                          true,'BoundingBox',[-100,50;-50,70]);
    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];
    setm(map_ax,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);
    geoshow(lat, lon, 'DisplayType','Point','Marker','x','Color','red');
    
    lon = num2str(lon);
    lat = num2str(lat);
    
    % Title the whole visualization
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

function [index] = coord_to_closest_coords_index(all_coords, chosen_coord)
    %disp(['Finding nearest coordinates to choice']);
    coord_diff = chosen_coord - all_coords;
    index = 1;
    min_norm = norm(coord_diff(1,:));
    for k = 1:length(coord_diff)
        n = norm(coord_diff(k,:));
        if n < min_norm
            min_norm = n;
            index = k;
        end
    end
end

function [lon_bounds, lat_bounds] = get_region(name)
    
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
