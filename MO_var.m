
dir_sic=dir(['~/sic_data_2007/'])
sic_mat = [];
date_vec = [];

% Read data in from directory
for ifile=3:length(dir_sic)
    fname=dir_sic(ifile).name;
    fdate=fname(15:22);
    fyear=fdate(1:4);
    fmonth=fdate(5:6);
    fday=fdate(7:8);
    fdate2=([ fyear '-' fmonth '-' fday]);
    t=datetime(fdate2,'InputFormat','yyyy-MM-dd');
    doy_tmp=day(t,'dayofyear');
    date_vec=[date_vec doy_tmp];
    tmpdata=load(['~/sic_data_2007/' fname]);
    sic_day=tmpdata(:,3);
    sic_mat=[sic_mat sic_day];
end

window = 5;

keyboard;

time_series_doy=date_vec;
days = date_vec;
sic_std_mat=[];
sic_mean_mat=[];

 % creating a time series of daily variance for every day in the series
 sic_std= nan(size(days));
 sic_mean= nan(size(days));
 for k = 1:length(days);
     logical = ((days(k)- time_series_doy)>-1);
     if(sum(logical)>2 & sum(logical<=5))
     daily_sic =sic_mat(:,logical);
     sic_std = std(daily_sic');                 
     sic_mean = mean(daily_sic');               
     sic_std_mat=[sic_std_mat sic_std'];
     sic_mean_mat=[sic_mean_mat sic_mean'];
     end 
 end
 
% creating an array of possible melt onset dates based on different thresholds
rd = linspace(0, 1); % 500 used thresholds
MO_index(1:size(sic_mean_mat,1),1:length(rd))=NaN;
for factor = 1:length(rd);
    acc_range = rd(factor);
    for j = 1:size(sic_mean_mat,1)
    for k = 1:size(sic_mean_mat,2)
        if sic_mean_mat(j,k) < acc_range && isnan(MO_index(j,factor));
           MO_index(j,factor)=days(k);
        end
    end
    end
end

%keyboard;

iqr_lim=20;
for j=1:size(sic_mean_mat,1)
  %MO_index_nonan = MO_index(~isnan(MO_index(j,:)));
  MO_index_nonan = MO_index(~isnan(MO_index(j,:)));
  MO_index_nonan = MO_index(j,:);
  filter = MO_index_nonan > 61 & MO_index_nonan  <300; % 61 corresponds to AHRA algorithm
  MO_filter = MO_index_nonan(filter);
  inter_q(j) = iqr(MO_filter);
  MO(j) = quantile(MO_filter, 0.25);
  %hist_ratio = sum(MO_index(j,:) > 61 & MO_index(j,:)  <200)/sum(MO_index(j,:)  <200);
  %if inter_q > iqr_lim || hist_ratio <0.5
  %    MO(j) = nan;
  %end
  if(mod(j,50)==0)
     keyboard;
  end
end

function [] = create_figure()
end