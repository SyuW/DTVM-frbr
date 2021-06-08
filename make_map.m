% Hudson: [50 70] [-100 -50]
% NWT: [68 78] [-140 -110]

shown_points = {[75 -119] [75 -132] [71 -119] [71 -132]};
plot_map_and_graphs([68 78], [-140 -110], shown_points);

function [] = plot_map_and_graphs(region_lat_bounds, region_long_bounds, shown_points)

    axesm('mercator','MapLatLimit',region_lat_bounds,'MapLonLimit',region_long_bounds);
    setm(gca,'mlabellocation',6,'plabellocation',3,'MLabelParallel','south');

    lon_min = region_long_bounds(1);
    lon_max = region_long_bounds(2);
    lat_min = region_lat_bounds(1);
    lat_max = region_lat_bounds(2);

    mlabel on;
    plabel on;
    framem on;
    axis off;

    %disp([region_lat_bounds;region_long_bounds])

    provinces = shaperead('canada_shape_files/PROVINCE.SHP','UseGeoCoords',true,'BoundingBox',[lon_min,lat_min;lon_max,lat_max]);

    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];

    setm(gca,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);

    hold on;
    for k = 1:length(shown_points)
        pt = shown_points{k};
        lat = pt(1);
        lon = pt(2);
        geoshow(lat, lon, 'DisplayType','Point','Marker','x','Color','red');
        txt = strcat('(',num2str(lat),',',num2str(lon),')');
        text(lat, lon, txt);
    end

    saveas(gca, '~/scratch/MO_outputs/region_map.png');
    close;

end
