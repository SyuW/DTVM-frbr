%%%%%% Create maps of breakup/freeze-up days %%%%%%
plot_frbr_date_maps = 1;
if plot_frbr_date_maps
    % Plot names
    names = {'NRC Freeze-up days','NRC Breakup days',...
             'DTVM Freeze-up days','DTVM Breakup days',...
             'DTVM-NRC Breakup Difference',...
             'DTVM-NRC Freeze-up Difference'};
    % Maps
    map_mats = {fr_days_NRC,br_days_NRC,...
                fr_days_DTVM,br_days_DTVM,...
                br_days_DTVM-br_days_NRC,...
                fr_days_DTVM-fr_days_NRC};
                
    create_frbr_dates_maps(map_mats, names, coords);
    disp('Done creating freeze-up/breakup maps');
end



%%%%% Plot time series for randomly chosen locations %%%%%
plot_ts_for_rand_pts = 0;
if plot_ts_for_rand_pts
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
        plot_ts_and_map({sic_mat sic_mean_mat sic_std_mat},...
                        {'SIC' 'Mean SIC' 'SIC Std deviation'},... 
                        chosen_index,...
                        chosen_coord,...
                        date_vec,...
                        {freezeup_days_NRC,breakup_days_NRC,freezeup_days_DTVM,breakup_days_DTVM},...
                        0);
    end
    disp('Done plotting time-series for random points');
end
          %Lat_bounds Lon_bounds
regions = {{[51 55], [-83 -78]},... % James Bay
           {[55 64], [-96 -77]},... % Hudson Bay
           {[64 70], [-85 -75]},... % Foxe Basin
           {[58 65], [-77 -66]},... % Hudson Strait
           {[63 70], [-65 -50]},... % Baffin Sea
           {[50 62], [-64 -50]}};   % East Coast

custom_region = {[-64 -56] [54 60]};
region = custom_region; %regions{6};

lat_bounds = region{1};
lon_bounds = region{2};

%%%%% Plot time-series and histograms for 'spurious points' %%%%%
plot_ts_for_region = 1;
% Controlling parameters
is_breakup = 1;
created = 0;
creation_limit = 3;
day_difference_cutoff = 80;

s = RandStream('mlfg6331_64');
inds = randsample(s,length(coords),length(coords));

if plot_ts_for_region
    for k = 1:length(coords)
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
                    plot_ts_and_map({sic_mat sic_mean_mat sic_std_mat},...
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
    disp('Done plotting time-series for spurious/regular points and generating histograms');
end

% Clear everything to conserve memory
clear;

function [index] = coord_to_closest_coords_index(all_coords, chosen_coord)
    %disp(['Finding nearest coordinates to choice']);
    coord_diff = chosen_coord - all_coords;
    index = 1;
    min_norm = norm(coord_diff(1,:));
    for k = 1:length(coord_diff)
        c = coord_diff(k);
        n = norm(coord_diff(k,:));
        if n < min_norm
            min_norm = n;
            index = k;
        end
    end
end


function [] = plot_ts_and_map(mats, names, loc_index, location, days, dates_cells, spurious)

    lon = location(1);
    lat = location(2);

    colors = {'r', 'g', 'c', 'm'};
    day_type_names = ["NRC Freeze-up", "NRC Breakup", "DTVM Freeze-up", "DTVM Breakup"];
    figure;
    axs = [subplot(3,2,2) subplot(3,2,4) subplot(3,2, [5 6])];
    for k = 1:length(axs)
        hold(axs(k),'on');
        for a = 1:length(dates_cells)
            found_date = dates_cells{a}(loc_index);
            if ~isnan(found_date)
                plot(axs(k),[found_date found_date],[0 1],'Color',colors{a},'DisplayName',day_type_names{a});
            end
        end
        plot(axs(k), days, mats{k}(loc_index,:),'Color','b','DisplayName','Signal');
        xlabel(axs(k),'Days of year'); ylabel(axs(k),names{k});
        xlim(axs(k), [0 365]);
    end

    legend(axs(3),'Location','best');

    % Create map
    ax4 = subplot(3,2,[1 3]);

    axesm('mercator','MapLatLimit',[50 70],'MapLonLimit',[-100 -50]);
    setm(ax4,'mlabellocation',15,'plabellocation',5,'MLabelParallel','south');

    mlabel on;
    plabel on;
    framem on;
    axis off;

    provinces = shaperead('canada_shape_files/PROVINCE.SHP','UseGeoCoords',...
                          true,'BoundingBox',[-100,50;-50,70]);

    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];

    setm(ax4,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);
    geoshow(lat, lon, 'DisplayType','Point','Marker','x','Color','red');
    title(['SIC at (' num2str(lon) ',' num2str(lat) ')']);

    if spurious
        save_fname = strcat(output_directory,'plots/', 'spurious_point_plot_at_',num2str(lon), '_', num2str(lat),'.png');
    else
        save_fname = strcat(output_directory,'plots/', 'plot_at_',num2str(lon), '_', num2str(lat),'.png');
    end

    saveas(ax4, save_fname);
    close all;
end


function [] = create_frbr_dates_maps(frbr_cell_arr, names, coords_mat)
    global output_directory;

    lons = coords_mat(:,1);
    lats = coords_mat(:,2);

    for k = 1:length(frbr_cell_arr)
        scatter(lons,lats,10,frbr_cell_arr{k},'filled');
        colorbar;
        if k <= 4
            caxis([0 365]);
        else
            colormap(jet);
            caxis([-100 100]);
        end
        title(names{k});
        xlabel('Longitude');
        ylabel('Latitude');
        save_fname = strcat(output_directory,'maps/',names{k},'.png');
        saveas(gca, save_fname);
    end
    close all; 
end


function [] = plot_histogram(X, location, is_breakup)
    global output_directory;

    lon = location(1);
    lat = location(2);
    br_range = [60 306];
    fr_range = [245 365];

    if is_breakup
        day_type_str = 'Breakup';
        X = X(X > br_range(1) & X < br_range(2));
        if ~isempty(X)
            flagged_date = quantile(X, 0.75);
        else
            flagged_date = br_range(1);
        end
    else
        day_type_str = 'Freezeup';
        X = X(X > fr_range(1) & X < fr_range(2));
        if ~isempty(X)
            flagged_date = quantile(X, 0.25);
        else
            flagged_date = fr_range(2);
        end
    end

    histogram(X, 1:365, 'Normalization', 'probability');
    hold on;
    plot(gca,[flagged_date flagged_date],[0 1]);

    fig_title = strcat(day_type_str,' histogram',' at',' (',num2str(lon),',',num2str(lat),')');
    title(fig_title);
    save_fname = strcat(output_directory,'histograms/',day_type_str,'_histogram','_at',num2str(lon),'_',num2str(lat),'.png');

    xlabel('Day of year');
    ylabel('Normalized frequency');

    padding=10;
    if ~isnan(X)
        xlim([min(X)-padding max(X)+padding]);
    end
    ylim([0 1]);

    saveas(gca, save_fname);

    close all;
end