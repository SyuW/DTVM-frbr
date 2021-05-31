function [] = plot_map_and_graphs(points)

    figure;

    axesm('mercator','MapLatLimit',[50 70],'MapLonLimit',[-100 -50]);
    setm(gca,'mlabellocation',10,'plabellocation',5);

    mlabel on;
    plabel on;
    framem on;
    axis off;

    provinces = shaperead('canada_shape_files/PROVINCE.SHP','UseGeoCoords',true,'BoundingBox',[-100,50;-50,70]);

    LandColor = [0.298 0.6 0];
    OceanColor = [0.8 1 1];

    setm(gca,'FFaceColor',OceanColor);
    geoshow(provinces,'FaceColor',LandColor);
    saveas(gca, 'map.png');

end
