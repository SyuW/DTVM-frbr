tic;
main_exec_sic_gif("2009", 0.15);
toc;


function [] = main_exec_sic_gif(chosen_year, threshold)
    % Main entry point of sea ice concentration (SIC) GIF generation
    %
    % arguments:
    %   chosen_year - chosen year for processing data into GIF
    %   threshold - threshold for binarizing the sea ice concentration
    %
    % return: None
    %
    % loaded variables: None
    % saved variables: None

    data_src = chosen_year+"_esacci/";
    sic_dir = "./out/"+data_src+"raw/mats/";
    
    save_fname = "GIFs/Hudson_ice_presence_in_"+chosen_year+".gif";
    create_sic_gif(chosen_year, sic_dir, save_fname, 0.1, threshold);
end


function [] = create_sic_gif(chosen_year, sic_dir, save_fname,...
                             gif_delay, th)
    % Generates GIFs of sea ice concentration
    % 
    % arguments:
    %   chosen_year
    %   sic_dir
    %   save_fname
    %   gif_delay
    %   th 
    %
    % return: None
    % 
    % loaded variables:
    %   sic_mats
    %       sic_mat - sea ice concentration matrix
    %   coords
    %       coords - matrix containing coordinates
    %
    % saved variables: None
    
    % convert year string to number for rest of calc.
    chosen_year = str2double(chosen_year);
    
    % load sic data and also coordinate information
    load(sic_dir+"sic.mat", "sic_mat");
    load(sic_dir+"coords.mat", "coords");
    
    % Thresholding
    if th ~= 0
        % binarize SIC values at threshold
        sic_mat = sic_mat > th;
        sic_mat = single(sic_mat);
    end
    
    % get total number of days in year
    days_in_year = size(sic_mat,2);
    
    % get longitudes/latitudes from coordinates matrix
    lons = coords(:,1);
    lats = coords(:,2);
    
    % instantiate a new figure
    fig = figure("visible", "off");
    colormap;
    
    % Scatter SIC for first day
    markersize = 10;
    sic_date = datetime(chosen_year, 1, 1);
    s = scatter(lons, lats, markersize, sic_mat(:,1), "filled");
    
    % Set title and x,y axis labels
    title(string(sic_date));
    xlabel("Longitude");
    ylabel("Latitude");
    
    % Get frame for the initial day
    f = getframe(fig);
    height = size(f.cdata, 1);
    width = size(f.cdata, 2);
    
    % preallocate movie matrix using frame dimensions
    mov = zeros(height, width, 1, days_in_year, "uint8");
    
    % update movie matrix and get color map of initial day
    [mov(:,:,1,1), map] = rgb2ind(f.cdata, 256, "nodither");
    
    % Create the rest of the animation by updating color map
    for d=2:days_in_year
        
        % Updating the title
        sic_date = datetime(chosen_year, 1, d);
        title(string(sic_date));
        
        % Update the color map of the scatter object
        s.CData = sic_mat(:,d);
        
        % Get frame of GIF from figure
        f = getframe(fig);
        
        % Create the rest of the frames
        mov(:,:,1,d) = rgb2ind(f.cdata, map, "nodither");
        
    end
    
    % Create animated GIF
    imwrite(mov, map, save_fname, "DelayTime", gif_delay, "LoopCount", inf)
    
    clf;
end

