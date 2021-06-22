global output_directory
output_directory = '~/scratch/dtvm_outputs/';

calc_window = 5;

dir_sic=dir(['~/sic_data_2007/']);
sic_mat = [];
date_vec = [];

% Read data in from directory
for ifile=3:length(dir_sic)
    fname=dir_sic(ifile).name;
    tmpdata=load(['~/sic_data_2007/' fname]);

    fdate=fname(15:22);
    fyear=fdate(1:4);
    fmonth=fdate(5:6);
    fday=fdate(7:8);
    fdate2=([ fyear '-' fmonth '-' fday]);
    t=datetime(fdate2,'InputFormat','yyyy-MM-dd');
    doy_tmp=day(t,'dayofyear');
    current_sic = tmpdata(:,3);
    date_vec=[date_vec doy_tmp];
    %disp(['~/sic_data_2007/' fname])
    coords=tmpdata(:,[1 2]);
    sic_mat=[sic_mat current_sic];
end

% Linearly interpolate SIC signal
sic_mat = interp1(date_vec, sic_mat', date_vec(1):date_vec(end))';
date_vec = date_vec:date_vec(end);

mats_savename = strcat(output_directory,'out/calc_mats');

% Apply a median filter before binarization
remove_fluctuations_before_bin = 0;
filter_order = 3;
if remove_fluctuations_before_bin
    disp('Removing short term fluctuations in SIC signal before binarization');
    for loc = 1:size(sic_mat, 1)
        sic_mat(loc,:) = medfilt1(double(sic_mat(loc,:)), filter_order, 'truncate');
    end
    mats_savename = strcat(mats_savename, '_filtered');
end

% Binarizing the signal for mean, std deviation calculation
bin_sic = 1;
sic_threshold = 0.15;
if bin_sic
    disp('Binarizing the SIC signal before calculating mean and standard deviation');
    sic_mat = sic_mat > sic_threshold;
    mats_savename = strcat(mats_savename, '_bin_sic');
else
    disp('Continuing without binarizing');
end

% Apply a median filter after binarization
remove_fluctuations_after_bin = 1;
filter_order = 3;
if remove_fluctuations_after_bin
    disp('Removing short term fluctuations in SIC signal after binarization');
    for loc = 1:size(sic_mat, 1)
        sic_mat(loc,:) = medfilt1(double(sic_mat(loc,:)), filter_order, 'truncate');
    end
    mats_savename = strcat(mats_savename, '_filtered');
end

[sic_std_mat, sic_mean_mat] = create_mean_and_std(date_vec,sic_mat,calc_window);

% Binarize mean
bin_mean = 0; 
if bin_mean
    sic_mean_mat = sic_mean_mat > bin_mean;
    mats_savename = strcat(mats_savename, '_mean');
    disp('Binarizing the mean SIC signal');
end

% Binarize std deviation
bin_std = 0;
if bin_std
    sic_std_mat = sic_std_mat > bin_std;
    mats_savename = strcat(mats_savename, '_std')
    disp('Binarizing the std SIC signal');
end

disp(strcat('Writing workspace variables to ',mats_savename))
save(mats_savename,'sic_mat','sic_std_mat','sic_mean_mat','date_vec','coords');
clearvars;

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