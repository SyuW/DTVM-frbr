yearstr =int2str(2009)

for imonth=1:12

monthstr=int2str(imonth)

if(imonth<10)
  dir_sic=dir(['~/Hudson_Strait/ESACCI/' yearstr '/0' monthstr '/*.nc']);
end

if(imonth>=10)
  dir_sic=dir(['~/Hudson_Strait/ESACCI/' yearstr '/' monthstr '/*.nc']);
end

% Sams debugging
dir_sic

for ifile=1:length(dir_sic)

fname=dir_sic(ifile).name
fdate=fname(45:52);

if(imonth<10)
  ncid = netcdf.open([ yearstr '/0' monthstr '/' fname ]);
end
if(imonth>=10)
  ncid = netcdf.open([ yearstr '/' monthstr '/' fname ]);
end

%ncid = netcdf.open(['ESACCI-SEAICE-L4-SICONC-AMSR_25.0kmEASE2-NH-' fdate '-fv2.1.nc']);
%info = ncinfo('ESACCI-SEAICE-L4-SICONC-AMSR_25.0kmEASE2-NH-20140117-fv2.1.nc');

lat_esa=netcdf.getVar(ncid,5);
lon_esa=netcdf.getVar(ncid,6);
ice_conc=netcdf.getVar(ncid,7);

raw_ice_conc=netcdf.getVar(ncid,8);
total_error=netcdf.getVar(ncid,9);
smearing_error=netcdf.getVar(ncid,10);
algorithm_error=netcdf.getVar(ncid,11);
status_flag=netcdf.getVar(ncid,12);

latmin=45;
latmax=90;
lonmin=-180;
lonmax=180.0;

ikeep=0;
for i=1:size(lat_esa,1)
for j=1:size(lat_esa,2)
    if(lat_esa(i,j)<latmax & lat_esa(i,j)>latmin)
    if(lon_esa(i,j)<lonmax & lon_esa(i,j)>lonmin)
    if(ice_conc(i,j)>=0 & status_flag(i,j) ~= 1)
      ikeep=ikeep+1;
      latk(ikeep)=lat_esa(i,j);
      lonk(ikeep)=lon_esa(i,j);
      ice_conck(ikeep)=double(ice_conc(i,j))/10000;
      raw_ice_conck(ikeep)=double(raw_ice_conc(i,j))/10000;
      total_errork(ikeep)=total_error(i,j);
      smearing_errork(ikeep)=double(smearing_error(i,j))/10000;
      algorithm_errork(ikeep)=algorithm_error(i,j);
      status_flagk(ikeep)=status_flag(i,j);
      if(raw_ice_conck(ikeep)>-3)
         retrieved_ice_conck(ikeep)=raw_ice_conck(ikeep);
      end
      if(raw_ice_conck(ikeep)<=-3)
         retrieved_ice_conck(ikeep)=ice_conck(ikeep);
      end
    end
    end 
    end
end
end

%bin 0 - land
%bin 1 - lake
%bin 2 - weather filter
%bin 3 - land spill-over

fid2=fopen(['ice_esa_sicci_' fdate '.dat'],'w');

for i=1:length(lonk)
    fprintf(fid2,'%14.6e',lonk(i));
    fprintf(fid2,'%14.6e',latk(i));
    fprintf(fid2,'%14.6e',ice_conck(i));
    fprintf(fid2,'%14.6e',raw_ice_conck(i));
    fprintf(fid2,'%14.6e',retrieved_ice_conck(i));
    fprintf(fid2,'\n');
end

fclose(fid2)

end

end
