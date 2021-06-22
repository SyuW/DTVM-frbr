% sic_mat, sic_mean_mat, sic_std_mat
% date_vec, coords

use_binarized = 1;
if use_binarized
    disp('Using the binarized SIC signal');
    load('~/scratch/dtvm_outputs/out/calc_mats_bin_sic');
else
    disp('Using the unbinarized SIC signal')
    load('~/scratch/dtvm_outputs/out/calc_mats');
end

% Controllable parameters
num_of_thresholds = 500;
iqr_lim = 20;
% Calculate freeze-up/breakup days using the DTVM method
[freezeup_days_DTVM breakup_days_DTVM BR_index FR_index] = DTVM_freezeup_breakup(date_vec, sic_std_mat, num_of_thresholds, iqr_lim);

% Calculate freeze-up/breakup days using the NRC method
freezeup_days_NRC = cts_presence_breakup_freezeup(sic_mat, date_vec, 'Freeze-up');
breakup_days_NRC = cts_presence_breakup_freezeup(sic_mat, date_vec, 'Breakup');

%%%%%% Create maps of breakup/freeze-up days %%%%%%
plot_frbr_date_maps = 1;
if plot_frbr_date_maps
    % Also plot the differences as maps
    names = {'NRC Freeze-up days','NRC Breakup days','DTVM Freeze-up days','DTVM Breakup days',...
             'DTVM-NRC Breakup Difference','DTVM-NRC Freeze-up Difference'};
    create_frbr_dates_maps({freezeup_days_NRC,breakup_days_NRC,...
                            freezeup_days_DTVM,breakup_days_DTVM,...
                            breakup_days_DTVM-breakup_days_NRC,...
                            freezeup_days_DTVM-freezeup_days_NRC},...
                           {'NRC Freeze-up days','NRC Breakup days','DTVM Freeze-up days','DTVM Breakup days',...
                            'DTVM-NRC Breakup Difference','DTVM-NRC Freeze-up Difference'},...
                           coords);
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

%%%%% Plot time-series and histograms for 'spurious points' %%%%%
plot_ts_for_region = 1;
% Controlling parameters
is_breakup = 1;
created = 0;
creation_limit = 3;
day_difference_cutoff = 50;
lat_bounds = [52 65]; lon_bounds = [-60 -50];
if plot_ts_for_region
    if is_breakup
        diffs = breakup_days_DTVM-breakup_days_NRC;
        freezeup_or_breakup_index = BR_index;
    else
        diffs = freezeup_days_DTVM-freezeup_days_NRC;
        freezeup_or_breakup_index = FR_index;
    end

    clearvars FR_index BR_index;

    for k = 1:length(coords)
        coord = coords(k,:);
        % Determine if location falls within bounds of region of interest
        if coord(1) > lon_bounds(1) && coord(1) < lon_bounds(2)...
        && coord(2) > lat_bounds(1) && coord(2) < lat_bounds(2)
            day_diff = diffs(k);
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
                    days_for_histogram = freezeup_or_breakup_index(k,:);
                    plot_histogram(days_for_histogram, coord, is_breakup);
                    created = created + 1;
                    if created >= creation_limit
                        break
                    end
                end
            end
        end
    end
    disp('Done plotting time-series for spurious points and generating histograms');
end

% Clear variables to conserve memory
clearvars;

function [day_of_interest] = cts_presence_breakup_freezeup(mat, days, day_type)
    period = 15;
    threshold = 0.15;
    day_of_interest = nan(1,size(mat,1));
    if day_type == "Freeze-up"
        mat = mat > threshold;
        season_start_end = [245 length(days)];
    elseif day_type == "Breakup"
        mat = mat < threshold;
        season_start_end = [60 306];
    else
        error('Error. Day type input is not correct.');
    end
    for loc = 1:size(mat,1)
        binarized_signal_at_loc = mat(loc,:) > threshold;
        for d = 1:length(days)-period+1
            if season_start_end(1) < days(d) && days(d) < season_start_end(2)
                if sum(binarized_signal_at_loc(d:d+period-1))==period
                    day_of_interest(loc) = days(d);
                    break
                end
            end
        end
    end
end
       

% Dynamic Threshold Variability Method algorithm for flagging freeze-up/breakup dates
function [FR,BR,BR_index,FR_index] = DTVM_freezeup_breakup(days, mat, num_of_thresholds, iqr_lim, is_fr)
    br_range = [60 250];
    fr_range = [245 365];

    % Create the thresholds and identify dates based on the thresholds
    BR_index(1:size(mat,1),1:num_of_thresholds)=nan;
    FR_index(1:size(mat,1),1:num_of_thresholds)=nan;

    % th = threshold, loc = location, d = day
    for loc = 1:size(mat, 1)
        max_val = max(mat(loc,:));
        thresholds_vec = linspace(0, max_val, num_of_thresholds);

        % Calculate possible breakup days by counting forwards
        for th = 1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(BR_index(loc, th));
                % Iterate over days to find when threshold exceeded
                for d = 1:days(end)-1
                    threshold_exceeded = (mat(loc,d) >= threshold) ||...
                                         (mat(loc,d) > threshold & mat(loc,d+1) < threshold);
                    if threshold_exceeded
                        BR_index(loc,th)=days(d);
                        break
                    end
                end
            end
        end

        % Filter out NaN dates and dates outside of season
        BR_dates_at_loc = BR_index(loc,:);
        BR_dates_at_loc = BR_dates_at_loc(~isnan(BR_dates_at_loc));
        BR_dates_at_loc = BR_dates_at_loc(BR_dates_at_loc > br_range(1) & BR_dates_at_loc < br_range(2));
        if ~isempty(BR_dates_at_loc)
            % Use 75th percentile for breakup to favor later days
            BR(loc) = quantile(BR_dates_at_loc, 0.75);
        else
            % Set the breakup day to beginning of season
            BR(loc) = br_range(1);
        end

        % Calculate possible freeze-up days by counting backwards
        for th=1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(FR_index(loc, th));
                % Iterate over days starting from freeze-up season start
                for df = days(end):-1:2
                    threshold_exceeded = (mat(loc,df) >= threshold) ||...
                                         (mat(loc,df) > threshold & mat(loc,df-1) < threshold);
                    if threshold_exceeded
                        FR_index(loc,th)=days(df);
                        break
                    end
                end
            end
        end

        % Filter out NaN dates and dates outside of season
        FR_dates_at_loc = FR_index(loc,:);
        FR_dates_at_loc = FR_dates_at_loc(~isnan(FR_dates_at_loc));
        FR_dates_at_loc = FR_dates_at_loc(FR_dates_at_loc > fr_range(1) & FR_dates_at_loc < fr_range(2));
        if ~isempty(FR_dates_at_loc)
            % Use 25th percentile for freeze-up to favor earlier days 
            FR(loc) = quantile(FR_dates_at_loc, 0.25);
        else
            % Set the freeze-up day to the end of the season
            FR(loc) = fr_range(2);
        end
    end
disp('Done DTVM for freeze-up/breakup dates');
end


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
    day_type_names = {"NRC Freeze-up", "NRC Breakup", "DTVM Freeze-up", "DTVM Breakup"};
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

    legend(axs(3),'Location','northwest');

    % Create map
    ax4 = subplot(3,2,[1 3]);

    axesm('mercator','MapLatLimit',[50 70],'MapLonLimit',[-100 -50]);
    setm(ax4,'mlabellocation',15,'plabellocation',5,'MLabelParallel','south');

    mlabel on;
    plabel on;
    framem on;
    axis off;

    provinces = shaperead('canada_shape_files/PROVINCE.SHP','UseGeoCoords',true,'BoundingBox',[-100,50;-50,70]);

    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];

    setm(ax4,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);
    geoshow(lat, lon, 'DisplayType','Point','Marker','x','Color','red');
    title(['SIC at (' num2str(lon) ',' num2str(lat) ')']);

    if spurious
        save_fname = strcat('~/scratch/dtvm_outputs/plots/', 'spurious_point_plot_at_',num2str(lon), '_', num2str(lat),'.png');
    else
        save_fname = strcat('~/scratch/dtvm_outputs/plots/', 'plot_at_',num2str(lon), '_', num2str(lat),'.png');
    end

    saveas(ax4, save_fname);
    close all;
end


function [] = create_frbr_dates_maps(frbr_cell_arr, names, coords_mat)
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
        save_fname = strcat('~/scratch/dtvm_outputs/maps/',names{k},'.png');
        saveas(gca, save_fname);
    end
    close all; 
end


function [] = plot_histogram(X, location, is_breakup)
    lon = location(1);
    lat = location(2);
    br_range = [60 250];
    fr_range = [245 365];

    if is_breakup
        day_type_str = 'Breakup';
        X = X(X > br_range(1) & X < br_range(2));
        flagged_date = quantile(X, 0.75);
    else
        day_type_str = 'Freezeup';
        X = X(X > fr_range(1) & X < fr_range(2));
        flagged_date = quantile(X, 0.25);
    end

    edges = 1:365;
    histogram(X, edges, 'Normalization', 'probability');
    hold on;
    plot(gca,[flagged_date flagged_date],[0 1]);

    fig_title = strcat(day_type_str,' histogram',' at',' (',num2str(lon),',',num2str(lat),')');
    title(fig_title);
    save_fname = strcat('~/scratch/dtvm_outputs/histograms/',day_type_str,'_histogram','_at',num2str(lon),'_',num2str(lat),'.png');

    xlabel('Day of year');
    ylabel('Normalized frequency');
    ylim([0 1]);

    saveas(gca, save_fname);

    close all;
end


