tic
main_exec_sic_gif(2010);
toc

function [] = main_exec_sic_gif(chosen_year)
    data_src = strcat(num2str(chosen_year),'_esacci/');
    sic_dir = strcat('./dtvm_outputs/', data_src, 'out/');
    
    save_fname = strcat('GIFs/Hudson_SIC_in_',num2str(chosen_year),'.gif');
    create_sic_gif(chosen_year, sic_dir, save_fname, 0.1);
end

function [] = create_sic_gif(chosen_year, sic_dir, save_fname, gif_delay)

    load(strcat(sic_dir, 'sic_mats.mat'), 'sic_mat');
    load(strcat(sic_dir, 'coords.mat'), 'coords');
    
    days_in_year = size(sic_mat,2);
    
    lons = coords(:,1);
    lats = coords(:,2);

    fig = figure('visible','off');
    colormap;
    hold on;
    
    % Scatter SIC for first day
    sic_date = datetime(chosen_year, 1, 1);
    s = scatter(lons,lats,2,sic_mat(:,1),'filled');
    
    % Set title and x,y axis labels
    title(string(sic_date));
    xlabel('Longitude');
    ylabel('Latitude');
    
    % Preallocate movie matrix
    pos = get(gcf, 'Position');
    width = pos(3);
    height = pos(4);
    mov = zeros(height, width, 1, days_in_year, 'uint8');
    
    % Get frame for the initial day
    f = getframe(gcf);
    [mov(:,:,1,1), map] = rgb2ind(s.CData, 256, 'nodither');
    
    % Create the rest of the animation by updating color map
    for d=2:days_in_year
        
        % Updating the title
        sic_date = datetime(chosen_year, 1, d);
        title(string(sic_date));
        
        % Update the color map of the scatter object
        s.CData = sic_mat(:,d);
        
        % Get frame of GIF from figure
        f = getframe(gcf);
        
        % Create the rest of the frames
        mov(:,:,1,d) = rgb2ind(f.cdata, map, 'nodither');
        
    end
    
    % Create animated GIF
    imwrite(mov, map, save_fname, 'DelayTime', gif_delay, 'LoopCount', inf)
    
    clf;
end

