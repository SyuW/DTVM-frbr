plot_map_and_graphs()

function [] = plot_map_and_graphs()

    % Create remaining subplots
    subplot(3,2,2);
    subplot(3,2,4);
    subplot(3,2,[5 6]);

    % Create the map
    subplot(3,2,[1 3]);

    axesm('mercator','MapLatLimit',[50 70],'MapLonLimit',[-100 -50]);
    setm(gca,'mlabellocation',10,'plabellocation',5,'MLabelParallel','south');

    mlabel on;
    plabel on;
    framem on;
    axis off;

    provinces = shaperead('canada_shape_files/PROVINCE.SHP','UseGeoCoords',true,'BoundingBox',[-100,50;-50,70]);

    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];

    setm(gca,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);

    saveas(gca, '~/scratch/MO_outputs/map.png');

end
