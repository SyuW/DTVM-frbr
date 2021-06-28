% Create maps of breakup/freeze-up days
plot_frbr_date_maps = 0;
if plot_frbr_date_maps
    % Also plot the differences as maps
    flagged_dates_cells{5} = breakup_days_DTVM-breakup_days_NRC;
    flagged_dates_cells{6} = freezeup_days_DTVM-freezeup_days_NRC;
    names = {'NRC Freeze-up days','NRC Breakup days','DTVM Freeze-up days','DTVM Breakup days',...
             'DTVM-NRC Breakup Difference','DTVM-NRC Freeze-up Difference'};
    create_frbr_dates_maps(flagged_dates_cells,names,coords);
end

% Plot time series for randomly chosen locations
plot_ts_for_rand_pts = 0;
if plot_ts_for_rand_pts
    flagged_dates_cells = {freezeup_days_NRC,breakup_days_NRC,...
                           freezeup_days_DTVM,breakup_days_DTVM};
    matrices = {sic_mat sic_mean_mat sic_std_mat};
    names = {'SIC' 'Mean SIC' 'SIC Std deviation'};

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
                        {freezeup_days_NRC,breakup_days_NRC,freezeup_days_DTVM,breakup_days_DTVM});
    end
end


plot_ts_for_region = 1;
breakup_diffs = breakup_days_DTVM-breakup_days_NRC;
freezeup_diffs = freezeup_days_DTVM-freezeup_days_NRC;
created = 0;
if plot_ts_for_region
    lat_bounds = [52 65]; lon_bounds = [-60 -50];
    for k = 1:length(coords)
        coord = coords(k,:);
        if coord(1) > lon_bounds(1) && coord(1) < lon_bounds(2)...
        && coord(2) > lat_bounds(1) && coord(2) < lat_bounds(2)
            day_diff = freezeup_diffs(k);
            if ~isnan(day_diff)
                if abs(day_diff) >= 20
                    %disp(coord); disp(k);
                    plot_ts_and_map({sic_mat sic_mean_mat sic_std_mat},...
                                    {'SIC' 'Mean SIC' 'SIC Std deviation'},...
                                    k,...
                                    coord,...
                                    date_vec,...
                                    {freezeup_days_NRC,breakup_days_NRC,freezeup_days_DTVM,breakup_days_DTVM})
                    created = created + 1
                    if created > 8
                        break
                    end
                end
            end
        end
    end
end

function [] = plot_ts_and_map(mats, names, loc_index, location, days, dates_cells)
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

    save_fname = strcat('~/scratch/dtvm_outputs/plots/', 'freezeup_over_plot_at_',num2str(lon), '_', num2str(lat),'.png');
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
        save_fname = strcat('~/scratch/dtvm_outputs/experiments/',names{k},'.png');
        saveas(gca, save_fname);
    end
    close all; 
end