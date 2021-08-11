% ----------------------------------------------------------- %
% --------------------- Main function ----------------------- %
% ----------------------------------------------------------- %

% Call the execution of DTVM method for freeze-up/breakup
% Important - choose the data source you want to apply DTVM for
% Options: 2007_esacci/, 2008_esacci/, 2009_esacci/, 2010_esacci/
tic;
%DTVM_main_exec("2007_esacci/", "raw");
process_batch()
toc;

function [] = process_batch()
    % Utility function for processing multiple sources of data
    % as well as multiple types of ways it's processed
    %
    % arguments: None
    
    data_srcs = ["2007_esacci/","2008_esacci/"];
    process_types = ["hysteresis","binfilt","raw"];

    for j = 1:length(data_srcs)
        for k = 1:length(process_types)        
            data_src = data_srcs(j);
            process_type = process_types(k);
            fprintf("Source: %s, Process type: %s\n", data_src, process_type);
            DTVM_main_exec(data_src, process_type);
            fprintf("Done freezeup/breakup calculation\n");
        end
    end
end

function [] = DTVM_main_exec(data_src, process_type)
    % Entry point of execution of DTVM method
    % arguments:
    %   data_src - string describing which data source to use
    %       allowed: (2007_esacci/, 2008_esacci/, 
    %                 2009_esacci/, 2010_esacci/)
    %   process_type - string describing how the data was processed
    %
    % return: None
    %
    % loaded variables:
    %   sic_mat - 2D matrix of sea ice concentrations (SIC)
    %   sic_std_mat - 2D matrix of moving std deviation of SIC
    %
    % saved variables: None
    
    % type checking
    valid_process_types = ["binfilt","hysteresis","raw"];
    if ~any(strcmp(valid_process_types, process_type))
        error(['Error. %s is not a valid processed type.',...
               'Please choose from raw, binfilt, hysteresis'], process_type);
    end
    
    fprintf("Starting freeze-up/breakup calculation\n");
    
    % output directory for dtvm outputs
    out_dir = "./out/"+data_src+process_type+"/dtvm/";
    
    % input directory for loading SIC data
    mats_dir = "./out/"+data_src+process_type+"/mats/";
    
    % Make the output directory if it doesn't exist
    if not(isfolder(out_dir))
        mkdir(out_dir);
    end
    
    % load SIC matrix and also filtered SIC matrix
    load(mats_dir+"sic", "sic_mat");
    load(mats_dir+"filt_sic", "filt_sic");

    % create NRC + DTVM frbr dates -- this is the standard function
    period = 15;
    create_NRC_DTVM_frbr(out_dir, sic_mat, filt_sic, period);
    
    % create NRC frbr dates for varying periods
    disp("Starting NRC freeze-up/breakup gen. for window range");
    window_range = 5:30;
    create_NRC_frbr_for_window_range(out_dir, sic_mat, window_range);
    disp("Done generating NRC freeze-up/breakup data for window range");
end

% ----------------------------------------------------------------------- %
% --------------------- Data processing Functions ----------------------- %
% ----------------------------------------------------------------------- %

function [] = create_NRC_DTVM_frbr(out_dir, sic, filt_sic, period)
    % Entry point of execution of DTVM method
    % arguments (input):
    %   out_dir - base output directory path (e.g.
    %   ./dtvm_outputs/2007_esacci/dtvm/)
    %   sic - matrix containing SIC signals (1D arrays) for each
    %   location
    %   filt_sic - matrix containing filtered SIC signals (1D arrays) for
    %   each location
    %
    % arguments (output): None
    %
    % loaded variables: None
    %
    % saved variables:
    %   NRC_frbr_dates
    %       br_days_NRC - vector of NRC breakup days
    %       fr_days_NRC - vector of NRC freezeup days
    %   DTVM_frbr_dates
    %       br_days_DTVM - vector of DTVM breakup days
    %       fr_days_DTVM - vector of DTVM freezeup days
    %   DTVM_frbr_indexes
    %       BR_index - 2D matrix of potential breakup days per location
    %       FR_index - 2D matrix of potential freezeup days per location
    
    % Number of thresholds to use with DTVM algorithm
    num_of_thresholds = 500;
    
    % Calculate NRC freezeup/breakup dates vector
    fr_days_NRC = NRC_freezeup_breakup(sic, "Freeze-up", period);
    br_days_NRC = NRC_freezeup_breakup(sic, "Breakup", period);
    
    disp("Done creating NRC freeze-up/breakup dates");
    
    % Calculate DTVM freezeup/breakup dates vectors/indexes
    [fr_days_DTVM, br_days_DTVM, BR_index, FR_index]... 
        = DTVM_freezeup_breakup(filt_sic, num_of_thresholds);
    
    % Save freeze-up/breakup dates to output directory
    save(out_dir+"NRC_frbr_dates", "br_days_NRC", "fr_days_NRC")
    save(out_dir+"DTVM_frbr_dates", "br_days_DTVM", "fr_days_DTVM");
    save(out_dir+"DTVM_frbr_indexes", "BR_index", "FR_index");
      
    disp("Done saving DTVM freeze-up/breakup dates");   
end

function [] = create_NRC_frbr_for_window_range(out_dir, sic_mat, window_range)
    % Calculate NRC freeze-up/breakup for a range of window values
    % arguments (input)
    %   out_dir - base output directory path (e.g.
    %   ./dtvm_outputs/2007_esacci/dtvm/)
    %   sic_mat - 2D matrix of SIC values for (location, day of year)
    %   window_range - vector of window values
    %   process_type - process type
    %
    % arguments (output): None
    %
    % loaded variables: None
    %
    % saved variables:
    %   NRC_frbr_cubes
    %       NRC_br_cube
    %       NRC_fr_cube
    %
    
    start_period = window_range(1);
    end_period = window_range(end);
    
    % preallocate cubes to store data
    NRC_br_cube = nan(window_range(end), size(sic_mat,1));
    NRC_fr_cube = nan(window_range(end), size(sic_mat,1));
    
    % iterate over periods 
    for period = window_range
    
        % Calculate freeze-up/breakup days using the NRC method for period
        fr_days_NRC = NRC_freezeup_breakup(sic_mat, "Freeze-up", period);
        br_days_NRC = NRC_freezeup_breakup(sic_mat, "Breakup", period);
        
        % add to cube
        NRC_br_cube(period,:) = br_days_NRC;
        NRC_fr_cube(period,:) = fr_days_NRC;
    
    end
    
    % save the data cubes
    save_fname = "NRC_frbr_cubes";
    save(out_dir+save_fname, "NRC_br_cube", "NRC_fr_cube");
    
end

% --------------------------------------------------------------- %
% --------------------- FrBr Functions -------------------------- %
% --------------------------------------------------------------- %

function [frbr_dates] = NRC_freezeup_breakup(mat, day_type, period)
    % Freeze-up/Breakup dates using NRC definition
    % arguments:
    %   mat - 2D matrix of SIC values for (location, day of year)
    %   day_type - string describes whether freeze-up/breakup considered
    %   period - integer. how many days of continuous presence needed
    %   
    % return:
    %   frbr_dates - vector of freezeup or breakup dates
    %
    % loaded variables: None
    %
    % saved variables: None
    
    % Controllable parameters
    threshold = 0.15;
    
    % Number of coordinates
    num_of_locations = size(mat, 1);
    
    % Pre-allocate dates array
    frbr_dates = nan(1,num_of_locations);
    
    % Binarize for ice presence:
    % {0 = water, 1 = ice} if 'Freeze-up'
    % {1 = water, 0 = ice} if 'Breakup'
    if day_type == "Freeze-up"
        mat = mat > threshold;
        season_start_end = [245 365];
    elseif day_type == "Breakup"
        mat = mat < threshold;
        season_start_end = [60 306];
    else
        error("Error. Day type input is not correct.");
    end
    
    % Find freeze-up/breakup date at each location
    for loc = 1:num_of_locations
        presence_at_loc = mat(loc,:);
        for d = 1:365-period+1
            if season_start_end(1) < d && d < season_start_end(2)
                % If detected ${period} consecutive days of water
                % Set freeze-up/breakup day to beginning of that period
                if sum(presence_at_loc(d:d+period-1))==period
                    frbr_dates(loc) = d;
                    break
                end
            end
        end
    end
end

function [FR,BR,BR_index,FR_index] = DTVM_freezeup_breakup(mat, num_of_thresholds)
    % Freeze-up/breakup dates using the DTVM algorithm
    % arguments:
    %   mat - input matrix, usually should be SIC variability
    %   num_of_thresholds - integer. number of threshold values for DTVM
    %
    % return:
    %   FR - vector of DTVM freezeup dates
    %   BR - vector of DTVM breakup dates
    %   FR_index - 2D matrix of potential DTVM
    %
    % loaded variables: None
    %
    % saved variables: None

    % Breakup season range
    br_range = [60 258];
    
    % Freeze-up season range
    fr_range = [259 365];
    
    % Calculate moving standard deviation of signal
    calc_window = 5;
    std_mat = movstd(mat, [calc_window-1, 0], 0, 2);
    
    % Number of coordinates
    num_of_locations = size(std_mat, 1);
    
    % Pre-allocate freeze-up/breakup mats
    BR_index(1:num_of_locations,1:num_of_thresholds)=nan;
    FR_index(1:num_of_locations,1:num_of_thresholds)=nan;
    BR(1:num_of_locations)=nan;
    FR(1:num_of_locations)=nan;

    % loc = location, th = threshold, d = day
    for loc = 1:num_of_locations
        max_val = max(std_mat(loc,:));
        thresholds_vec = linspace(0, max_val, num_of_thresholds);

        % Calculate possible breakup days by counting forwards
        for th = 1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(BR_index(loc, th))
                % Iterate over days to find when threshold exceeded
                % "Jump over fluctuations for dates earlier than start of breakup season"
                for d = br_range(1):365-1
                    threshold_exceeded = (std_mat(loc,d) >= threshold);
                    if threshold_exceeded
                        BR_index(loc,th)=d;
                        break
                    end
                end
            end
        end

        % Filter out NaN dates and breakup dates outside of breakup season
        BR_dates_at_loc = BR_index(loc,:);
        BR_dates_at_loc = BR_dates_at_loc(~isnan(BR_dates_at_loc));
        BR_dates_at_loc = BR_dates_at_loc(BR_dates_at_loc > br_range(1) &...
                                          BR_dates_at_loc < br_range(2));
        if ~isempty(BR_dates_at_loc)
            % Use 75th percentile for breakup to favor later days
            BR(loc) = quantile(BR_dates_at_loc, 0.75);
        else
            % Set the breakup day to beginning of season
            %BR(loc) = br_range(1);
        end

        % Calculate possible freeze-up days by counting backwards
        for th=1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(FR_index(loc, th))
                % Iterate backwards over days starting from freeze-up season end
                for df = 365:-1:2
                    threshold_exceeded = (std_mat(loc,df) >= threshold);
                    if threshold_exceeded
                        FR_index(loc,th)=df;
                        break
                    end
                end
            end
        end

        % Filter out NaN dates and freeze-up dates outside of freeze-up season
        FR_dates_at_loc = FR_index(loc,:);
        FR_dates_at_loc = FR_dates_at_loc(~isnan(FR_dates_at_loc));
        FR_dates_at_loc = FR_dates_at_loc(FR_dates_at_loc > fr_range(1) &...
                                          FR_dates_at_loc < fr_range(2));
        if ~isempty(FR_dates_at_loc)
            % Use 25th percentile for freeze-up to favor earlier days 
            FR(loc) = quantile(FR_dates_at_loc, 0.25);
        else
            % Set the freeze-up day to the end of the season
            %FR(loc) = fr_range(2);
        end
    end
end