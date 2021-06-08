% Controllable parameters
function [] = create_all_mats()
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

    days = date_vec;
    [sic_std_mat, sic_mean_mat] = create_mean_and_std(date_vec,sic_mat,calc_window);

    days_size = size(days);
    fid = fopen(strcat('days_',num2str(days_size(1)),'_',num2str(days_size(2)),'_',.mat),'w');
    fwrite(fid,days);

    sic_size = size(sic_mat);
    fid = fopen(strcat('sic_',num2str(sic_size(1)),'_',num2str(sic_size(2)),'_',.mat),'w');
    
end


% creating a time series of daily variance for every day in the series
function [sic_std_mat,sic_mean_mat] = create_mean_and_std(date_vec, sic, calc_window, threshold)
    cutoff = 0.15 % Option for binarizing the signal
    if threshold
        sic = sic > cutoff
    end
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
end