% ----------------------------------------------------------- %
% --------------------- Main function ----------------------- %
% ----------------------------------------------------------- %

% Call the execution of DTVM method for freeze-up/breakup
% Important - choose the data source you want to apply DTVM for
% Options: 2007_esacci/, 2008_esacci/, 2009_esacci/, 2010_esacci/
tic;
DTVM_main_exec("2007_esacci/", 1);
toc;

function [] = DTVM_main_exec(data_src, binfilt)
    % Entry point of execution of DTVM method
    % arguments:
    %   data_src - string describing which data source to use
    %       allowed: (2007_esacci/, 2008_esacci/, 
    %                 2009_esacci/, 2010_esacci/)
    %
    % return: None
    %
    % loaded variables:
    %   sic_mat - 2D matrix of sea ice concentrations (SIC)
    %   sic_std_mat - 2D matrix of moving std deviation of SIC
    %
    % saved variables: None
    
    % output directory for dtvm outputs
    out_dir = "./dtvm_outputs/"+data_src+"dtvm/";
    
    % input directory for loading SIC data
    mats_dir = "./dtvm_outputs/"+data_src+"out/";
    
    % Make the output directory if it doesn't exist
    if not(isfolder(out_dir))
        mkdir(out_dir);
    end
    
    % Load the sea ice concentration matrix created by process_data.m
    if binfilt
        load(mats_dir+"sic_mats_binarized_filtered", "sic_mat", "sic_std_mat");
    else
        load(mats_dir+"sic_mats", "sic_mat", "sic_std_mat");
    end
    
    % create NRC + DTVM frbr dates -- this is the standard function 
    create_NRC_DTVM_frbr(out_dir, sic_mat, sic_std_mat, 15, binfilt);
    
    % create NRC frbr dates for varying periods
    %window_range = 5:30;
    %create_NRC_frbr_for_window_range(out_dir, sic_mat, 5:30, binfilt);
end

% ----------------------------------------------------------------------- %
% --------------------- Data processing Functions ----------------------- %
% ----------------------------------------------------------------------- %

function [] = create_NRC_DTVM_frbr(out_dir, sic_mat, sic_std_mat, period, binfilt)
    % Entry point of execution of DTVM method
    % arguments (input):
    %   out_dir - base output directory path (e.g.
    %   ./dtvm_outputs/2007_esacci/dtvm/)
    %   sic_mat - matrix containing SIC signals (1D arrays) for each
    %   location
    %   sic_std_mat - matrix containing SIC variability signals
    %   binfilt - boolean for whether processed/raw data is being used
    %
    % arguments (output): None
    %
    % loaded variables: None
    %
    % saved variables:
    %   NRC_frbr_dates(_binfilt)
    %       br_days_NRC - vector of NRC breakup days
    %       fr_days_NRC - vector of NRC freezeup days
    %   DTVM_frbr_dates(_binfilt)
    %       br_days_DTVM - vector of DTVM breakup days
    %       fr_days_DTVM - vector of DTVM freezeup days
    %   DTVM_frbr_indexes(_binfilt)
    %       BR_index - 2D matrix of potential breakup days per location
    %       FR_index - 2D matrix of potential freezeup days per location
    
    % Controllable parameters
    num_of_thresholds = 500;
    
    % Calculate NRC freezeup/breakup dates vector
    fr_days_NRC = cts_presence_breakup_freezeup(sic_mat, "Freeze-up", period);
    br_days_NRC = cts_presence_breakup_freezeup(sic_mat, "Breakup", period);
    disp("Done creating NRC freeze-up/breakup dates");
    
    % Calculate DTVM freezeup/breakup dates vectors/indexes
    [fr_days_DTVM, br_days_DTVM, BR_index, FR_index]... 
        = DTVM_freezeup_breakup(sic_std_mat, num_of_thresholds);
    
    disp("Writing freeze-up/breakup data to " + out_dir);
    
    % if using processed signal, save under a different file name
    if binfilt
        disp("Using binarized and filtered signal for freeze-up/breakup date calculation");
        save(out_dir+"NRC_frbr_dates_binfilt", "br_days_NRC", "fr_days_NRC")
        save(out_dir+"DTVM_frbr_dates_binfilt", "br_days_DTVM", "fr_days_DTVM");
        save(out_dir+"DTVM_frbr_indexes_binfilt", "BR_index", "FR_index");
        disp("Done creating DTVM freeze-up/breakup dates");
        
    % using regular signal    
    else
        disp("Using raw signal for freeze-up/breakup date calculation");
        save(out_dir+"NRC_frbr_dates", "br_days_NRC", "fr_days_NRC")
        save(out_dir+"DTVM_frbr_dates", "br_days_DTVM", "fr_days_DTVM");
        save(out_dir+"DTVM_frbr_indexes", "BR_index", "FR_index");
        disp("Done creating DTVM freeze-up/breakup dates");
    end
    
    % Save NRC freezeup/breakup dates to file
    save(out_dir+"NRC_frbr_dates", "br_days_NRC", "fr_days_NRC");
end

function [] = create_NRC_frbr_for_window_range(out_dir, sic_mat, window_range, binfilt)
    % Calculate NRC freeze-up/breakup for a range of window values
    % arguments (input)
    %   out_dir - base output directory path (e.g.
    %   ./dtvm_outputs/2007_esacci/dtvm/)
    %   sic_mat - 2D matrix of SIC values for (location, day of year)
    %   window_range - vector of window values
    %   binfilt - boolean for whether raw/processed data is used
    %
    % arguments (output): None
    %
    % loaded variables: None
    %
    % saved variables:
    %   NRC_frbr_p_{window}:
    %       fr_days_NRC - NRC freeze-up dates for {window} value
    %       br_days_NRC - NRC breakup dates for {window} value
    
    % Iterate over periods
    for period = window_range
    
        % Calculate freeze-up/breakup days using the NRC method for period
        fr_days_NRC = cts_presence_breakup_freezeup(sic_mat, "Freeze-up", period);
        br_days_NRC = cts_presence_breakup_freezeup(sic_mat, "Breakup", period);
        
        % save to a different folder if using processed signal
        if binfilt
            save_folder = out_dir + "NRC_frbr_dates_varied_periods_binfilt/";
        else
            save_folder = out_dir + "NRC_frbr_dates_varied_periods/";
        end
        
        % save the created frbr dates vectors
        save_name = save_folder+"NRC_frbr_p_"+num2str(period);
        save(save_name,"fr_days_NRC","br_days_NRC");
        disp("Done creating NRC freeze-up/breakup dates for period of "...
             +num2str(period));
         
    end
    
end

% --------------------------------------------------------------- %
% --------------------- FrBr Functions -------------------------- %
% --------------------------------------------------------------- %

function [frbr_dates] = cts_presence_breakup_freezeup(mat, day_type, period)
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

    % Freeze-up/Breakup season range
    br_range = [60 306];
    fr_range = [245 365];
    
    % Number of coordinates
    num_of_locations = size(mat, 1);
    
    % Pre-allocate freeze-up/breakup mats
    BR_index(1:num_of_locations,1:num_of_thresholds)=nan;
    FR_index(1:num_of_locations,1:num_of_thresholds)=nan;
    BR(1:num_of_locations)=nan;
    FR(1:num_of_locations)=nan;

    % loc = location, th = threshold, d = day
    for loc = 1:num_of_locations
        max_val = max(mat(loc,:));
        thresholds_vec = linspace(0, max_val, num_of_thresholds);

        % Calculate possible breakup days by counting forwards
        for th = 1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(BR_index(loc, th))
                % Iterate over days to find when threshold exceeded
                % "Jump over fluctuations for dates earlier than start of breakup season"
                for d = br_range(1):365-1
                    threshold_exceeded = (mat(loc,d) >= threshold) ||...
                                         (mat(loc,d) > threshold && mat(loc,d+1) < threshold);
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
            BR(loc) = br_range(1);
        end

        % Calculate possible freeze-up days by counting backwards
        for th=1:length(thresholds_vec)
            threshold = thresholds_vec(th);
            if isnan(FR_index(loc, th))
                % Iterate backwards over days starting from freeze-up season end
                for df = 365:-1:2
                    threshold_exceeded = (mat(loc,df) >= threshold) ||...
                                         (mat(loc,df) > threshold && mat(loc,df-1) < threshold);
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