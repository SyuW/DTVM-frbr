% ----------------------------------------------------------- %
% --------------------- Main function ----------------------- %
% ----------------------------------------------------------- %

tic;
process_batch();
toc;

function [] = process_batch()
    % Utility function for processing multiple sources of data
    % as well as multiple types of ways it's processed
    %
    % arguments: None
    
    data_srcs = ["2007_esacci","2008_esacci","2009_esacci","2010_esacci/"];
    process_types = ["binfilt","raw"];

    for j = 1:length(data_srcs)
        for k = 1:length(process_types)        
            data_src = data_srcs(j);
            process_type = process_types(k);
            fprintf("Processing source, process type: %s, %s\n", data_src, process_type);
            process_data_main_exec(data_src, process_type);
            fprintf("\n");
        end
    end
end

function [] = process_data_main_exec(data_src, process_type)
    % Main entry point of execution
    %
    % arguments:
    %   data_src - string describing which data source to use
    %       allowed: (2007_esacci/, 2008_esacci/, 
    %                 2009_esacci/, 2010_esacci/)
    %
    % return: None
    %
    % saved variables:
    %   sic_mat - 2D matrix of sea ice concentrations (SIC)
    %   sic_mean_mat - 2D matrix of moving mean of SIC
    %   sic_std_mat - 2D matrix of moving std deviation of SIC
    %   coords - 2D matrix of coordinates (one coord per row)
    
    % output/data paths based on data source choice
    out_dir = "./out/" + data_src + process_type + "/mats/";
    data_dir = "./data/esacci_sic/" + data_src;
    
    % Make the output directory if it doesn't exist
    if not(isfolder(out_dir))
        mkdir(out_dir);
    end
    
    % Create sea ice concentration matrix
    sic_mat = create_sic_mat(data_dir);
    disp("Done constructing initial SIC matrix from data files");
    
    % Calculate moving mean/standard deviation of signal
    calc_window = 5;
    sic_mean_mat = movmean(sic_mat, [calc_window-1, 0], 2, "omitnan");
    sic_std_mat = movstd(sic_mat, [calc_window-1, 0], 0, 2, "omitnan");
    disp("Done calculating moving mean/std deviation of signal");
    
    % create second array for filtering
    filt_sic = sic_mat(:,:);
    
    % don't do anything to the signal    
    if process_type == "raw"
        disp("Not doing any processing");
        
    % binarize and apply filtering    
    elseif process_type == "binfilt"
        disp("Binarizing and median filtering SIC signal");
        filt_sic = binarize_signal(filt_sic, 0.15);
        % window size = 5 so that filter includes two neighboring points on
        % either side of five-point stencil
        % dim = 2 so that filtering is along each row
        filt_sic = filter_signal(filt_sic, 5, 2);
        
    % try using hysteresis
    elseif process_type == "hysteresis"
        disp("Applying hysteresis to SIC signal");
        % Apply hysteresis to time series
        for row = 1:size(filt_sic, 1)
            filt_sic(row,:) = hysteresis_binarize(filt_sic(row,:));
        end
        
    % invalid argument
    else
        error("Invalid process type %s", process_type)
 
    end
    
    disp("Done filtering");
    
    % Write all sea ice concentration related matrices to file
    disp(strcat("Writing SIC matrices to file"));
    save(out_dir+"sic", "sic_mat");
    save(out_dir+"sic_std", "sic_std_mat");
    save(out_dir+"sic_mean","sic_mean_mat");
    save(out_dir+"filt_sic", "filt_sic");
    
    % Get coordinate locations of all points
    coords = get_all_coordinates(data_dir);
    disp("Writing coordinates to file");
    save(out_dir+"coords","coords");
end

% ------------------------------------------------------------- %
% --------------------- Sub functions ------------------------- %
% ------------------------------------------------------------- %

function [sic_mat] = create_sic_mat(data_dir)
    % Create sea ice concentration matrix from .dat files
    %
    % arguments:
    %   data_dir - path string to data source
    %       example: './data/esacci_sic/2007_esacci/'
    %
    % return:
    %   sic_mat - 2D matrix of shape [num_of_locations 365]
    
    % Create listing of files inside data folder
    dir_sic = dir(data_dir+"*sic*");
    
    % Get number of locations and year
    first_fname = dir_sic(1).name;
    num_of_locations = size(load(data_dir+first_fname),1);
    year = str2double(first_fname(15:18));
    
    % use year to determine number of days
    division_by_four = mod(year, 4);
    if division_by_four == 0
        total_num_of_days = 366;
    else
        total_num_of_days = 365;
    end
    
    % days of year vector
    doy_vec = 1:total_num_of_days;
    
    % prellocating SIC matrix
    sic_mat = nan(num_of_locations, total_num_of_days);
    
    % Iterate over files to enter SIC values into matrix
    for ifile=1:length(dir_sic)
        fname=dir_sic(ifile).name;
        fdate=fname(15:22);
        
        t=datetime(fdate,"InputFormat","yyyyMMdd");
        doy=day(t,"dayofyear");

        tmpdata=load(data_dir+fname);
        current_day_sic=tmpdata(:,3);
        
        % fill matrix column with SIC data for day
        sic_mat(:,doy) = current_day_sic;
    end
    
    % Linearly interpolate
    for row=1:size(sic_mat,1)
        current_day_sic = sic_mat(row,:);
        
        doy_vec_i = doy_vec(~isnan(current_day_sic));
        current_day_sic_i = current_day_sic(~isnan(current_day_sic));
        
        % replace with interpolated SIC signal at location
        sic_mat(row,:) = interp1(doy_vec_i, current_day_sic_i, doy_vec);
    end
 
end

function [coords] = get_all_coordinates(data_dir)
    % Get all coordinates
    %
    % arguments:
    %   data_dir - path string to data source
    %       example: './data/esacci_sic/2007_esacci/'
    %
    % arguments:
    %   coords - matrix with each row representing a coord
    
    dir_sic = dir(data_dir+"*sic*");
    % grabbing coordinates
    first_file = load(data_dir+dir_sic(1).name);
    coords = first_file(:,[1 2]);
end

function [signal] = filter_signal(signal, window_size, dim)
    % Apply median filter to signal
    %
    % arguments:
    %   signal - 1D matrix representing signal to be filtered
    %   window_size - order of the median filter
    %   dim - which dimension of input to filter (should be 1)
    %
    % return:
    %   signal - 1D matrix representing filtered signal
    
    % convert to single in case signal is logical array
    signal = single(signal);
    
    % 'zeropad' param to consider signal to be zero beyond endpoints when
    % endpoint filtering
    signal = medfilt1(signal, window_size, [], dim, 'zeropad');
end

function [signal] = binarize_signal(signal, cutoff)
    % Binarize signal at cutoff
    %
    % arguments:
    %   signal - vector representing signal to be binarized
    %   cutoff - cutoff value
    %
    % return:
    %   signal - vector representing binarized signal
    
    signal = signal > cutoff;
    signal = single(signal);
end

function [signal] = hysteresis_binarize(signal)
    % Experimenting with hysteresis behavior
    % 
    % arguments:
    %   signal - vector representing signal to apply hysteresis to
    %
    % return:
    %   signal - vector representing 'hysteresized' signal
    
    % set ice state based on SIC on first day
    if signal(1) > 0.15
        ice = 1;
    else
        ice = 0;
    end
    
    % go through every day in the time series signal
    for day = 1:length(signal)
        % if it's ice, convert to water if SIC < 0.15
        if ice
            if signal(day) < 0.15
                ice = 0;
                signal(day) = 0;
            else
                ice = 1;
                signal(day) = 1; 
            end
        % if it became water, only convert back to ice if SIC > 0.75
        else
            if signal(day) < 0.75
                ice = 0;
                signal(day) = 0;
            else
                ice = 1;
                signal(day) = 1;
            end
        end
    end
end
           