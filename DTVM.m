% ----------------------------------------------------------- %
% --------------------- Main function ----------------------- %
% ----------------------------------------------------------- %

% Call the execution of DTVM method for freeze-up/breakup
DTVM_main_exec();

function [] = DTVM_main_exec()
    % load variables {sic_mat, sic_mean_mat, sic_std_mat, coords} to workspace
    out_directory = './dtvm_outputs/dtvm';
    mats_directory = './dtvm_outputs/out/';
    load(strcat(mats_directory, 'sic_mats_binarized_filtered'),'sic_mat','sic_std_mat');

    % Controllable parameters
    num_of_thresholds = 500;

    tic;
    % Calculate freeze-up/breakup days using the NRC method
    fr_days_NRC = cts_presence_breakup_freezeup(sic_mat, 'Freeze-up');
    br_days_NRC = cts_presence_breakup_freezeup(sic_mat, 'Breakup');
    disp('Done creating NRC freeze-up/breakup dates');
    toc;

    tic;
    % Calculate freeze-up/breakup days using the DTVM method
    [fr_days_DTVM, br_days_DTVM, BR_index, FR_index] = ...
    DTVM_freezeup_breakup(sic_std_mat, num_of_thresholds);
    disp('Done creating DTVM freeze-up/breakup dates');
    toc;
    
    disp(strcat('Writing freeze-up/breakup data to ', out_directory));
    save(strcat(out_directory,'NRC_frbr_dates'),'br_days_NRC','fr_days_NRC');
    save(strcat(out_directory,'DTVM_frbr_dates'),'br_days_DTVM','fr_days_DTVM');
    save(strcat(out_directory,'DTVM_frbr_indexes'),'BR_index','FR_index');
end

% --------------------------------------------------------- %
% --------------------- Functions ------------------------- %
% --------------------------------------------------------- %

% Freeze-up/Breakup dates using NRC definition
function [frbr_dates] = cts_presence_breakup_freezeup(mat, day_type)
    
    % Controllable parameters
    period = 15;
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
        error('Error. Day type input is not correct.');
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

% Dynamic Threshold Variability Method algorithm for flagging freeze-up/breakup dates
function [FR,BR,BR_index,FR_index] = DTVM_freezeup_breakup(mat, num_of_thresholds)

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
            FR(loc) = fr_range(2);
        end
    end
end