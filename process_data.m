calc_window = 5;

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

bin_sic = 0.15;

% Binarizing the signal for mean, std deviation calculation
% Set to zero if you don't want to binarize
if bin_sic
    disp('Binarizing the SIC signal before calculating mean and standard deviation');
    sic_mat = sic_mat > bin_sic;
end

[sic_std_mat, sic_mean_mat] = create_mean_and_std(date_vec,sic_mat,calc_window);

bin_mean = 0; 
if bin_mean
    sic_mean_mat = sic_mean_mat > bin_mean;
end

bin_std = 0;
if bin_std
    sic_std_mat = sic_std_mat > bin_std;
end

save('~/scratch/dtvm_outputs/out/calc_mats','sic_mat','sic_std_mat','sic_mean_mat','date_vec','coords');

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