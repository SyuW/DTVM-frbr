% ----------------------------------------------------------- %
% --------------------- Main function ----------------------- %
% ----------------------------------------------------------- %

% Call the execution of data processing
process_data_main_exec();

function [] = process_data_main_exec()
    tic
    output_directory = './dtvm_outputs/out/';
    data_directory = './data/esacci_sic_data_2007/';
    
    % Create sea ice concentration matrix
    sic_mat = create_sic_mat(data_directory);
    
    % Additional processing
    sic_mat = binarize_signal(sic_mat, 0.15);
    sic_mat = filter_signal(sic_mat, 5, 2);
    
    disp('Done constructing SIC matrix from data files');

    calc_window = 5;
    % Calculate moving mean/standard deviation of signal
    sic_mean_mat = movmean(sic_mat, [calc_window-1, 0], 2);
    sic_std_mat = movstd(sic_mat, [calc_window-1, 0], 0, 2);
    disp('Done calculating moving mean/std deviation of signal');
    
    % Change the following line for a different filename
    mats_save_path = strcat(output_directory, 'sic_mats_binarized_filtered');
    
    disp(strcat('Writing SIC matrices to file'));
    save(mats_save_path,'sic_mat','sic_std_mat','sic_mean_mat');
    
    % Get coordinate locations of all points
    coords = get_all_coordinates(data_directory);
    coords_save_path = strcat(output_directory, 'coords');
    disp('Writing coordinates to file')
    save(coords_save_path,'coords');
    
    toc
end

% ------------------------------------------------------------- %
% --------------------- Sub functions ------------------------- %
% ------------------------------------------------------------- %

% Create the sea ice concentration matrix from data files
function [sic_mat] = create_sic_mat(data_dir_path)
    % Create listing of files inside data folder
    dir_sic = dir([data_dir_path '*sic*']);
    
    % Get number of locations
    num_of_locations = size(load([data_dir_path dir_sic(1).name]),1);
    
    % prellocating SIC matrix
    sic_mat = nan(num_of_locations, 365);
    
    % Iterate over files to enter SIC values into matrix
    for ifile=1:length(dir_sic)
        fname=dir_sic(ifile).name;
        fdate=fname(15:22);
        
        t=datetime(fdate,'InputFormat','yyyyMMdd');
        doy=day(t,'dayofyear');

        tmpdata=load([data_dir_path fname]);
        current_day_sic=tmpdata(:,3);
        
        sic_mat(:,doy)=current_day_sic;
    end
    
    % Interpolate for NaN values
    sic_mat = inpaint_nans(sic_mat,0);
end

function [coords] = get_all_coordinates(data_dir_path)
    dir_sic = dir([data_dir_path '*sic*']);
    % grabbing coordinates
    first_file = load([data_dir_path dir_sic(1).name]);
    coords = first_file(:,[1 2]);
end

% Apply median filtering to signal
function [signal] = filter_signal(signal, window_size, dim)
    signal = single(signal);
    signal = medfilt1(signal, window_size, [], dim, 'zeropad');
end

% Binarize signal at cutoff
function [signal] = binarize_signal(signal, cutoff)
    signal = signal > cutoff;
end