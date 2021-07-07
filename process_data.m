% ----------------------------------------------------------- %
% --------------------- Main function ----------------------- %
% ----------------------------------------------------------- %

process_data_main_exec('2007_esacci/');

function [] = process_data_main_exec(data_src)
    % Entry point of execution when process_data.m is run
    % arguments (input):
    %   data_src - string describing which data source to use
    %       allowed: (2007_esacci/, 2008_esacci/, 
    %                 2009_esacci/, 2010_esacci/)
    %
    % arguments (output): None
    %
    % saved variables:
    %   sic_mat - 2D matrix of sea ice concentrations (SIC)
    %   sic_mean_mat - 2D matrix of moving mean of SIC
    %   sic_std_mat - 2D matrix of moving std deviation of SIC
    %   coords - 2D matrix of coordinates (one coord per row)

    % Start of timer
    tic
    
    % output/data paths based on data source choice
    output_directory = strcat('./dtvm_outputs/',data_src,'out/');
    data_directory = strcat('./data/esacci_sic/',data_src);
    
    % Create sea ice concentration matrix
    sic_mat = create_sic_mat(data_directory);
    
    % Additional processing
    sic_mat = binarize_signal(sic_mat, 0.15);
    sic_mat = filter_signal(sic_mat, 5, 2);
    disp('Done constructing SIC matrix from data files');

    % Calculate moving mean/standard deviation of signal
    calc_window = 5;
    sic_mean_mat = movmean(sic_mat, [calc_window-1, 0], 2);
    sic_std_mat = movstd(sic_mat, [calc_window-1, 0], 0, 2);
    disp('Done calculating moving mean/std deviation of signal');
    
    % Change the following line for a different filename
    mats_save_path = strcat(output_directory, 'sic_mats_binarized_filtered');
    
    % Write all sea ice concentration related matrices to file
    disp(strcat('Writing SIC matrices to file'));
    save(mats_save_path,'sic_mat','sic_std_mat','sic_mean_mat');
    
    % Get coordinate locations of all points
    coords = get_all_coordinates(data_directory);
    coords_save_path = strcat(output_directory, 'coords');
    disp('Writing coordinates to file');
    save(coords_save_path,'coords');
    
    % End of timer
    toc
end

% ------------------------------------------------------------- %
% --------------------- Sub functions ------------------------- %
% ------------------------------------------------------------- %

function [sic_mat] = create_sic_mat(data_dir)
    % Create sea ice concentration matrix from .dat files
    % arguments (input):
    %   data_dir - path string to data source
    %       example: './data/esacci_sic/2007_esacci/'
    %
    % arguments (output):
    %   sic_mat - 2D matrix of shape [num_of_locations 365]

    % Create listing of files inside data folder
    dir_sic = dir([data_dir '*sic*']);
    
    % Get number of locations
    num_of_locations = size(load([data_dir dir_sic(1).name]),1);
    
    % prellocating SIC matrix
    sic_mat = nan(num_of_locations, 365);
    
    % Iterate over files to enter SIC values into matrix
    for ifile=1:length(dir_sic)
        fname=dir_sic(ifile).name;
        fdate=fname(15:22);
        
        t=datetime(fdate,'InputFormat','yyyyMMdd');
        doy=day(t,'dayofyear');

        tmpdata=load([data_dir fname]);
        current_day_sic=tmpdata(:,3);
        
        sic_mat(:,doy)=current_day_sic;
    end
    
    % Interpolate NaN values and set negative values to zero
    sic_mat = inpaint_nans(sic_mat,0);
    sic_mat(sic_mat < 0) = 0;
end

function [coords] = get_all_coordinates(data_dir)
    % Get all coordinates
    % arguments (input):
    %   data_dir - path string to data source
    %       example: './data/esacci_sic/2007_esacci/'
    %
    % arguments (output):
    %   coords - matrix with each row representing a coord
    
    dir_sic = dir([data_dir '*sic*']);
    % grabbing coordinates
    first_file = load([data_dir dir_sic(1).name]);
    coords = first_file(:,[1 2]);
end

function [signal] = filter_signal(signal, window_size, dim)
    % Get all coordinates
    % arguments (input):
    %   signal - 1D matrix representing signal to be filtered
    %   window_size - order of the median filter
    %   dim - which dimension of input to filter (should be 1)
    %
    % arguments (output):
    %   signal - 1D matrix representing filtered signal
    
    signal = single(signal);
    signal = medfilt1(signal, window_size, [], dim, 'zeropad');
end

function [signal] = binarize_signal(signal, cutoff)
    % Binarize signal at cutoff
    % arguments (input):
    %   signal - 1D matrix representing signal to be binarized
    %   cutoff - cutoff value
    %
    % arguments (output):
    %   signal - 1D matrix representing binarized signal
    
    signal = signal > cutoff;
    signal = single(signal);
end