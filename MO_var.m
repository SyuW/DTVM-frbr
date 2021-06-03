% Data has the form
% longitude latitude SIC ... ...
dir_sic=dir(['~/sic_data_2007/']);
sic_mat = [];
date_vec = [];

% Read data in from directory
for ifile=3:length(dir_sic)
    fname=dir_sic(ifile).name;
    fdate=fname(15:22);
    fyear=fdate(1:4);
    fmonth=fdate(5:6);
    fday=fdate(7:8);
    fdate2=([ fyear '-' fmonth '-' fday]);
    t=datetime(fdate2,'InputFormat','yyyy-MM-dd');
    doy_tmp=day(t,'dayofyear');
    date_vec=[date_vec doy_tmp];
    tmpdata=load(['~/sic_data_2007/' fname]);
    sic_day=tmpdata(:,3);
    coords=tmpdata(:,[1 2]);
    sic_mat=[sic_mat sic_day];
end

days = date_vec;

[sic_std_mat, sic_mean_mat] = create_mean_and_std(date_vec,sic_mat,3);

matrices = {sic_mat sic_mean_mat sic_std_mat};
names = {'SIC' 'Mean SIC' 'SIC Std deviation'};

choice = [-85 60];
chosen_index = coord_to_closest_coords_index(coords, choice);
chosen_coord = coords(chosen_index,:);
plot_ts_and_map(matrices,names,chosen_index,chosen_coord,date_vec);

MO_index = find_possible_MO_dates(days, sic_std_mat);
determine_MO_date(MO_index, sic_std_mat);

% creating an array of possible melt onset dates based on different thresholds
% rd = linspace(0, 1); % 500 used thresholds
% MO_index(1:size(sic_mean_mat,1),1:length(rd))=days(1);
% for factor = 1:length(rd);
%    acc_range = rd(factor);
%    for j = 1:size(sic_mean_mat,1)
%        for k = 1:size(sic_mean_mat,2)
%            if sic_mean_mat(j,k) < acc_range && isnan(MO_index(j,factor));
%                MO_index(j,factor)=days(k);
%            end
%        end
%    end
%end

keyboard;

iqr_lim=20;
for j=1:size(sic_mean_mat,1)
    %MO_index_nonan = MO_index(~isnan(MO_index(j,:)));
    MO_index_nonan = MO_index(~isnan(MO_index(j,:)));
    MO_index_nonan = MO_index(j,:);
    filter = MO_index_nonan > 61 & MO_index_nonan  <300; % 61 corresponds to AHRA algorithm
    MO_filter = MO_index_nonan(filter);
    inter_q(j) = iqr(MO_filter);
    MO(j) = quantile(MO_filter, 0.25);
    %hist_ratio = sum(MO_index(j,:) > 61 & MO_index(j,:)  <200)/sum(MO_index(j,:)  <200);
    %if inter_q > iqr_lim || hist_ratio <0.5
    %    MO(j) = nan;
    %end
    if(mod(j,50)==0)
        keyboard;
    end
end


function [MO_index] = find_possible_MO_dates(days, mat)
    num_of_thresholds = 500;
    % Initialize matrix of MO dates for each location, threshold
    MO_index(1:size(mat,1),1:num_of_thresholds)=nan;
    % th = threshold, l = location, d = day
    for loc = 1:size(mat, 1)
        max_val = max(mat(loc,:));
        thresholds_vec = linspace(0, max_val, num_of_thresholds);
        for th = 1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(MO_index(loc, th));
                for d = 1:size(mat, 2)-1
                    threshold_exceeded = (mat(loc,d) >= threshold) | (mat(loc,d) < threshold & mat(loc,d+1) > threshold);
                    if threshold_exceeded
                        MO_index(loc,th)=days(d);
                        break
                    end
                end
            end
        end
    end
end


function [] = determine_MO_date(MO_index, mat)
    iqr_lim=20;
    for loc=1:size(mat,1)
        MO_index_nonan = MO_index(~isnan(MO_index(loc,:)));
        MO_index_nonan(loc,:);
        filter = MO_index_nonan > 61 & MO_index_nonan < 300;
        MO_filtered = MO_index_nonan(filter);

        keyboard;
    end
end


% creating a time series of daily variance for every day in the series
function [sic_std_mat,sic_mean_mat] = create_mean_and_std(date_vec, sic, window)
    sic_std_mat = [];
    sic_mean_mat = [];
    sic_std= nan(size(date_vec));
    sic_mean= nan(size(date_vec));
    num_of_locations = size(sic,1);
    for k = 1:length(date_vec);
        % Select only the current day and previous two days
        logical = ((date_vec(k)-date_vec)>-1 & (date_vec(k)-date_vec)<=window-1);
        if (sum(logical)>2 & sum(logical)<=5)
            days_sic = sic(:,logical);
            sic_std = std(days_sic');
            sic_mean = mean(days_sic');
            sic_std_mat = [sic_std_mat sic_std'];
            sic_mean_mat = [sic_mean_mat sic_mean'];
        else
            sic_std_mat = [sic_std_mat nan(num_of_locations, 1)];
            sic_mean_mat = [sic_mean_mat nan(num_of_locations, 1)];
        end
    end
end


function [index] = coord_to_closest_coords_index(all_coords, chosen_coord)
    disp(['Finding nearest coordinates to choice']);
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


function [] = plot_ts_and_map(mats, names, loc_index, location, days)

    
    lon = location(1);
    lat = location(2);

    figure;
    axs = [subplot(3,2,2) subplot(3,2,4) subplot(3,2, [5 6])];
    for k = 1:length(axs)
        plot(axs(k), days, mats{k}(loc_index,:));
        xlabel(axs(k),'Days of year'); ylabel(axs(k),names{k});
        xlim(axs(k), [0 365]);
    end

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
    title(['SIC at (' num2str(lon) ',' num2str(lat) ')'])

    saveas(ax4, '~/scratch/MO_outputs/plots.png');
    close;

end