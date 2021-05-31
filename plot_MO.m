load lons.mat
load lats.mat

figure(1)
load MO_2003.mat
scatter(lons,lats,20,MO,'filled');colorbar
title('2003');

figure(2)
load MO_2004.mat
scatter(lons,lats,20,MO,'filled');colorbar
set(gca,'clim',[90,220]);set(gca,'clim',[90,220]);
title('2004');

figure(3)
load MO_2005.mat
scatter(lons,lats,20,MO,'filled');colorbar
set(gca,'clim',[90,220]);set(gca,'clim',[90,220]);
title('2005');

figure(4)
load MO_2006.mat
scatter(lons,lats,20,MO,'filled');colorbar
set(gca,'clim',[90,220]);
title('2006');

figure(5)
load MO_2007.mat
scatter(lons,lats,20,MO,'filled');colorbar
set(gca,'clim',[90,220]);
title('2007');

figure(6)
load MO_2008.mat
scatter(lons,lats,20,MO,'filled');colorbar
set(gca,'clim',[90,220]);
title('2008');

figure(7)
load MO_2009.mat
scatter(lons,lats,20,MO,'filled');colorbar
set(gca,'clim',[90,220]);
title('2009');
