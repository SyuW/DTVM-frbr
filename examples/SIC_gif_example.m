dir_sic=dir(['~/sic_data_2007/*2007*.dat']);
sic_mat = [];
date_vec = [];

fig = figure;
for ifile=1:length(dir_sic)
    fname=dir_sic(ifile).name;

    fdate=fname(15:22);
    fyear=fdate(1:4);
    fmonth=fdate(5:6);
    fday=fdate(7:8);
    fdate2=([ fyear '/' fmonth '/' fday]);

    %t=datetime(fdate2,'InputFormat','yyyy-MM-dd');
    %doy_tmp=day(t,'dayofyear');

    tmpdata=load(['~/sic_data_2007/' fname]);
    lons=tmpdata(:,1);
    lats=tmpdata(:,2);
    sic_day=tmpdata(:,3);

    scatter(lons, lats, 20, sic_day);
    title(fdate2);
    xlabel('Longitude');
    ylabel('Latitude');
    %c = colorbar('southoutside');
    %c.Label.String = 'Sea ice concentration';
    drawnow;
    frame = getframe(fig);
    im{ifile} = frame2im(frame);
end
close;

gif_file_name = 'Hudson_SIC_in_2007.gif';
% gif_delay = 0.1 is pretty smooth
gif_delay = 0.1;

for idx = 1:length(dir_sic)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,gif_file_name,'gif','LoopCount',Inf,'DelayTime',gif_delay);
    else
        imwrite(A,map,gif_file_name,'gif','WriteMode','append','DelayTime',gif_delay);
    end
end


