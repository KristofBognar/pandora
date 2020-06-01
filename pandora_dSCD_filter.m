function [ind_bad_ci, ind_bad_o4, ind_bad_rms] = pandora_dSCD_filter(data,cols_in,rms_lim,smooth_window,N_min)
% [smooth_ci, smooth_o4, bad_rms] = pandora_dSCD_filter(data,rm_rms,smooth_window,N_min)
% filter pandora MAX-DOAS dSCDs based on rms and CI, O4 smoothness
%
% This finction only returns indices, it does not modify the data
%
% INPUT: 
%   data: Raw QDOAS data saved in matlab format
%   cols_in: cell, names of analysis windows (or names of RMS columns) in the file
%   rms_lim: RMS cutoff for all elevation angles
%            Single value, or same size as cols_in
%   smooth_window: length of smoothing window to use in CI, O4
%                  smoothing. Actual input is fraction of datapoints per day,
%                  this is calculated as the smooth_window divided by the length of
%                  the measurement period for the given day
%   N_min: minimum number of datapoints for smoothing. Days with fewer
%          points (or duration shorter than smooth_window) are flagges as bad data. 
%           
% OUTPUT:logical indices of bad ci, O4 and RMS values, same lengths as the input data
% 
%@Kristof Bognar, May 2020

 
if length(rms_lim)==1
    rms_lim=repmat(rms_lim,1,length(cols_in));
elseif length(rms_lim)~=length(cols_in)
    error('rms_lim must be either single value, or same size as cols_in')
end

disp('RMS limit used:')
disp(rms_lim)

% filter high RMS -- remove bad datapoints so they're not included in the
% smoothness check (or the retrievals)
ind_bad_rms=false(size(data,1),length(cols_in));
for i=1:length(cols_in)
    
    % find RMS column
    if strcmp(cols_in{i}(end-2:end),'RMS')
        tmp=cols_in{i};
    else
        tmp=[cols_in{i} 'RMS'];
    end
    
    % get column index
    col_ind=find(strcmp(data.Properties.VariableNames,tmp));
    if isempty(col_ind), error('RMS column not found'); end
    
    % filter RMS
    ind_bad_rms(:,i)=data{:,col_ind}>rms_lim(i);
    
end

% calculate color index
% wavelengths (but not this exact pair) from Wagner et al. 2014, 2016
ci=data.Fluxes330./data.Fluxes440;

% break up measurements by local day (can ignore daylight savings here)
try
    day_local=day(data.DateTime-hours(5),'dayofyear');
catch
    error('Check if files were processed the same way as in merge_pandora_dscds.m (DateTime column missing)')    
end
N_days=unique(day_local);

% array for smoothed data
smooth_ci=NaN(size(ci));
smooth_o4=NaN(size(ci));

% loop over all days
nn=0;
for i=1:length(N_days)
    
    % display progress info
    disp_str=['Filtering day ', num2str(i), '/', num2str(length(N_days))];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,nn));
    fprintf(disp_str);
    nn=numel(disp_str);    
    
    % loop over all elevation angles
    for ea=[1,2,3,5,8,10,15,20,30,40,50]
        
        % indices of given elev angle on current day 
        ind=(data.Elevviewingangle==ea & day_local==N_days(i));

        % skip given elev if all data is bad
        if sum(ind)==0, continue, end
        
        % only proceed if enough measurements are present
        tmp=data.DateTime(ind);
        day_duration=hours(tmp(end)-tmp(1));

        if sum(ind)>N_min && day_duration>smooth_window

            % calculate fraction of data to use in smoothing (to match smooth_window)
            frac=smooth_window/day_duration;

            % smoth data
            smooth_ci(ind)=smooth(data.Fractionalday(ind),ci(ind),frac,'rloess');
            smooth_o4(ind)=smooth(data.Fractionalday(ind),data.NO2_VisSlColo4(ind),frac,'rloess');
            
        end
    end
end

% get indices of outliers (threshold from Gielen et al. 2014, Zhao et al. 2019)
% use negative of good values to make sure NaN smooth data registers as
% bad data (NaN conparison is always false)
ind_bad_ci=~(abs((ci-smooth_ci)./smooth_ci)<0.1);
ind_bad_o4=~(abs((data.NO2_VisSlColo4-smooth_o4)./smooth_o4)<=0.2);

fprintf('\n')

end

%%% old testing and plotting bits
%
% load('/home/kristof/work/PANDORA/profiling/QDOAS_output/spec_v4/Pandora_103_vis_2018.mat')
%
% for ea=[1,2,3,5,8,10,15,20,30,40,50]
%     
%     ind=data.Elevviewingangle==ea;
%     figure
%     dscatter(data.SZA(ind),data.NO2_VisRMS(ind))
%     title(['elev: ' num2str(ea) ' deg'])
%     ylim([0,0.006])
%     grid on
%     
%     tmp=sum(data.NO2_VisRMS(ind)>rms_lim)/sum(ind);
%     disp(['Filtered ' num2str(tmp*100) '% at ea = ' num2str(ea)])
%     
% end
%
% 
% 
% figure
% ind=data.Elevviewingangle==30;
% ax1=subplot(211);
% plot(data.Fractionalday(ind),ci(ind),'b.'), hold on
% plot(data.Fractionalday(ind & ind_bad_ci),ci(ind & ind_bad_ci),'rx'), hold on
% plot(data.Fractionalday(ind & ind_bad_o4),ci(ind & ind_bad_o4),'ro'), hold on
% plot(data.Fractionalday(ind),smooth_ci(ind),'r-')
% 
% ax2=subplot(212);
% plot(data.Fractionalday(ind),data.NO2_VisSlColo4(ind),'b.'), hold on
% plot(data.Fractionalday(ind & ind_bad_o4),data.NO2_VisSlColo4(ind & ind_bad_o4),'ro'), hold on
% plot(data.Fractionalday(ind & ind_bad_ci),data.NO2_VisSlColo4(ind & ind_bad_ci),'rx'), hold on
% plot(data.Fractionalday(ind),smooth_o4(ind),'r-')
% 
% linkaxes([ax1,ax2],'x')
