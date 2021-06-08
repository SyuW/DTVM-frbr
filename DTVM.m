% Controllable parameters
calc_window = 5;
num_of_thresholds = 500;
iqr_lim = 20;

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
    %disp(['~/sic_data_2007/' fname])
    tmpdata=load(['~/sic_data_2007/' fname]);
    sic_day=tmpdata(:,3);
    coords=tmpdata(:,[1 2]);
    sic_mat=[sic_mat sic_day];
end

% Binarizing the signal for mean, std deviation calculation
cutoff = 0.15;
if cutoff
    disp('Binarizing the SIC signal before calculating mean and standard deviation');
    sic = sic > cutoff;
end

days = date_vec;
[sic_std_mat, sic_mean_mat] = create_mean_and_std(date_vec,sic_mat,calc_window);

freezeup_days = find_breakup_freezeup_from_continuous_presence(sic_mat, days, 'Freeze-up');
breakup_days = find_breakup_freezeup_from_continuous_presence(sic_mat, days, 'Breakup');
MO_dates_DTVM = DTVM(days, sic_std_mat, num_of_thresholds, iqr_lim);

matrices = {sic_mat sic_mean_mat sic_std_mat};
flagged_dates_cells = {freezeup_days breakup_days MO_dates_DTVM};
names = {'SIC' 'Mean SIC' 'SIC Std deviation'};

choices = {[-85 60], [-80 65], [-90 62], [-80 52], [-70 63]};
%choices = {[-85 60]};

for k = 1:length(choices)
    choice = choices{k};
    chosen_index = coord_to_closest_coords_index(coords, choice);
    chosen_coord = coords(chosen_index,:);
    plot_ts_and_map(matrices,names,chosen_index,chosen_coord,date_vec,flagged_dates_cells);
end

%iqr_lim=20;
%for j=1:size(sic_mean_mat,1)
    %MO_index_nonan = MO_index(~isnan(MO_index(j,:)));
%    MO_index_nonan = MO_index(~isnan(MO_index(j,:)));
%    MO_index_nonan = MO_index(j,:);
%    filter = MO_index_nonan > 61 & MO_index_nonan  <300; % 61 corresponds to AHRA algorithm
%    MO_filter = MO_index_nonan(filter);
%    inter_q(j) = iqr(MO_filter);
%    MO(j) = quantile(MO_filter, 0.25);
%    %hist_ratio = sum(MO_index(j,:) > 61 & MO_index(j,:)  <200)/sum(MO_index(j,:)  <200);
%    %if inter_q > iqr_lim || hist_ratio <0.5
%    %    MO(j) = nan;
%    %end
%    if(mod(j,50)==0)
%        keyboard;
%    end
%end

function [MO] = find_breakup_freezeup_from_continuous_presence(mat, days, day_type)
    period = 15;
    threshold = 0.15;
    MO = nan(1,size(mat,1));
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
                    MO(loc) = days(d);
                    break
                end
            end
        end
    end
end
        

% Dynamic Threshold Variability Method for flagging dates
function [MO] = DTVM(days, mat, num_of_thresholds, iqr_lim)
    % Create the thresholds and identify dates based on the thresholds
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
    for loc=1:size(mat,1)
        MO_dates_at_loc = MO_index(loc,:);
        % Consider only dates that are not NaN
        MO_dates_at_loc = MO_dates_at_loc(~isnan(MO_dates_at_loc));

        melt_range_start = 1; melt_range_end = 364;
        dates_inside_melt_range = MO_dates_at_loc(MO_dates_at_loc > melt_range_start & MO_dates_at_loc < melt_range_end);
        %dates_before_melt_range = MO_dates_at_loc(MO_dates_at_loc < 61);

        %inter_q(loc) = iqr(dates_inside_melt_range)
        if ~isempty(dates_inside_melt_range)
            MO(loc) = quantile(dates_inside_melt_range, 0.25);
        else
            MO(loc) = nan;
        end
    end
end


function [MO_index] = old_MO_dates(days, sic_mean_mat)
    %creating an array of possible melt onset dates based on different thresholds
    rd = linspace(0, 1, 500); % 500 used thresholds
    MO_index(1:size(sic_mean_mat,1),1:length(rd))=nan;
    for factor = 1:length(rd);
        acc_range = rd(factor);
        for j = 1:size(sic_mean_mat,1)
            for k = 1:size(sic_mean_mat,2)
                if sic_mean_mat(j,k) < acc_range && isnan(MO_index(j,factor));
                    MO_index(j,factor)=days(k);
                end
            end
        end
    end
end


% creating a time series of daily variance for every day in the series
function [sic_std_mat,sic_mean_mat] = create_mean_and_std(date_vec, sic, calc_window)

    sic_std_mat = [];
    sic_mean_mat = [];
    sic_std= nan(size(date_vec));
    sic_mean= nan(size(date_vec));
    num_of_locations = size(sic,1);
    for k = 1:length(date_vec);
        % Select only the current day and previous (window-1) days
        logical = ((date_vec(k)-date_vec)>-1 & (date_vec(k)-date_vec)<=calc_window-1);
        if (sum(logical)==calc_window)
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
    %keyboard;
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


function [] = plot_ts_and_map(mats, names, loc_index, location, days, dates_cells)
    lon = location(1);
    lat = location(2);

    colors = {'r', 'g', 'm'};
    day_type_names = {"NRC Freeze-up", "NRC Breakup", "DTVM"};
    figure;
    axs = [subplot(3,2,2) subplot(3,2,4) subplot(3,2, [5 6])];
    for k = 1:length(axs)
        hold(axs(k),'on');
        for a = 1:length(dates_cells)
            MO_date = dates_cells{a}(loc_index);
            if ~isnan(MO_date)
                plot(axs(k),[MO_date MO_date],[0 1],'Color',colors{a},'DisplayName',day_type_names{a});
            end
        end
        plot(axs(k), days, mats{k}(loc_index,:),'Color','b','DisplayName','Signal');
        xlabel(axs(k),'Days of year'); ylabel(axs(k),names{k});
        xlim(axs(k), [0 365]);
    end

    legend(axs(3));

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

    save_fname = strcat('~/scratch/MO_outputs/plots/', 'plot_at_',num2str(lon), '_', num2str(lat),'.png');
    saveas(ax4, save_fname);
    close;

end


% WIP for saving matrix to binary
function [sic_mat, date_vec] = create_sic_mat(load_existing)
    output_dir_str='~/scratch/MO_outputs/';
    output_dir = dir('~/scratch/MO_outputs/');
    mat_fname = "foo";
    for ifile=1:length(output_dir)
        fname = output_dir(ifile).name;
        disp(fname);
        C = strsplit(fname, '_');
        if contains(fname, 'SIC_mat_')
            mat_fname=fname;
            mat_size = [str2num(C{3}) str2num(C{4})];
        elseif contains(fname, 'found_days')
            mat_fname=fname;
            day_size = [str2num(C{3}) str2num(C{4})];
            disp(C);
        end
    end

    keyboard;

    fid = fopen(mat_fname);
    if fid == -1 | ~load_existing
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
            %disp(['~/sic_data_2007/' fname])
            tmpdata=load(['~/sic_data_2007/' fname]);
            sic_day=tmpdata(:,3);
            coords=tmpdata(:,[1 2]);
            sic_mat=[sic_mat sic_day];
        end
        
        sic_mat_size = size(sic_mat);
        day_size = size(date_vec);

        fname1 = strcat("~/scratch/MO_outputs/SIC_mat_",num2str(mat_size(1)),'_',num2str(mat_size(2)),'_.bin');
        fname2 = strcat("~/scratch/MO_outputs/found_days_",num2str(day_size(1)),"_",num2str(day_size(2)),"_.bin");

        save(fname1,'sic_mat');
        save(fname2, 'date_vec');
    end
end

