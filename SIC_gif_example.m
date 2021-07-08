tic
main_exec_sic_gif(2008);
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

    lons = coords(:,1);
    lats = coords(:,2);

    fig = figure('visible','off');

    days_in_year = size(sic_mat,2);
    im = cell(1,days_in_year);
    for d=1:days_in_year
        sic_date = datetime(chosen_year, 1, d);

        scatter(lons, lats, 20, sic_mat(:,d));
        title(string(sic_date));
        xlabel('Longitude');
        ylabel('Latitude');

        drawnow;
        frame = getframe(fig);
        im{d} = frame2im(frame);
    end
    clf;

    for idx = 1:days_in_year
        [A,map]=rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,save_fname,'gif',...
                    'LoopCount',inf,'DelayTime',gif_delay);
        else
            imwrite(A,map,save_fname,'gif',...
                    'WriteMode','append','DelayTime',gif_delay);
        end
    end
end

