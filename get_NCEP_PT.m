function [ P_out, T_out, alt_grid_out ] = get_NCEP_PT( dates_in_datetime, alt_grid )
%GET_NCEP_PT get NCEP pressure and temperature profiles on set altitude grid
%   
% NCEP data is in pressure coordinates, interpolate to altitude grid
% Code assumes fixed file format for NCEP data
% One function call handles one year only
%
% INPUT 
%       dates_in_datetime: datetime array of measurement times (individual
%                 measurements, mean time of sequence, etc; anything works)
%       alt_grid (optional): specify altitude levels for P,T interpolation.
%                 Default is the (approximate) mean altitude of each
%                 pressure level in the NCEP files
%                 UNITS must be meters!!
%
% OUTPUT: P,T output is daily, and values are repeated for input times on
%         the same day, so the dimensions match
%       P_out: Pressure in hPa, interpolated linearly
%       T_out: Temperature in K, interpolated linearly
%       alt_grid_out: altitude in meters corresponding to P_out and T_out
%                     (returns default if alt_grid is not passed in)
%
%@Kristof Bognar, 2020

%% setup
if nargin==1,
   
    % default altitude grid, based on mean height of each pressure level in
    % 2018 Toronto NCEP data
    alt_grid=[0.15,1.5,3,5.5,7,9,10.5,12,14,16,18.5,20.5]*1000; % in m
    
end

% return altitude grid
alt_grid_out=alt_grid;

% find NCEP files
ncep_loc='/home/kristof/aurora/reanalysis/ncep/';

% location ID (eur, egb, mxc, tor)
loc_str='tor';

% pressure grid in hPa
p_grid=[0.4,1,2,5,10,30,50,70,100,150,200,250,300,400,500,700,850,1000];

% otput arrays
P_out=NaN(length(dates_in_datetime),length(alt_grid));
T_out=NaN(length(dates_in_datetime),length(alt_grid));

%% convert times to days

% year of interest (data should have one year only)
yr=unique(year(dates_in_datetime));
if length(yr)>1, error('Call GET_NCEP_PT one year at a time'), end

% get days in yyyymmdd format
dates_in=str2num(datestr(dateshift(dates_in_datetime, 'start', 'day'),'yyyymmdd'));
days_in=unique(dates_in);

%% read files

% select required file
fname_alt=[ncep_loc 'HgtNMC_' loc_str '_' num2str(yr) '.dat'];
fname_T=[ncep_loc 'TempNMC_' loc_str '_' num2str(yr) '.dat'];

% read altitude data, and separate date column
alt_data=dlmread(fname_alt,'',22,0);
days_ncep=alt_data(:,1);
alt_data(:,1)=[];

% read temperature data, and separate date column
T_data=dlmread(fname_T,'',22,0);
days_T=T_data(:,1);
T_data(:,1)=[];

% check if files are the same length
if ~isequal(days_ncep,days_T), error(['NCEP files incomplete for ' num2str(yr)]), end

%% interpolate data

% loop over days in the file (1 line per day)
for i=1:length(days_ncep)
    
    % if current date is needed, find P and T
    if ismember(days_ncep(i),days_in)
       
        % find which dates correspond to current day 
        ind=(dates_in==days_ncep(i));
        
%         % get pressure using barometric law
%         not great, deviates from model P -- use linear interpolation instead
%         tmp_fit=fit(alt_data(i,:)',p_grid','exp1'); % pressure at any altitude level
%         tmp_P=tmp_fit(alt_grid)'; % pressure on selected altitude grid
        
        % get pressure
        tmp_P=interp1(alt_data(i,:),p_grid,alt_grid,'linear','extrap');
        P_out(ind,:)=repmat(tmp_P,sum(ind),1);
        
        % get temperature
        tmp_T=interp1(alt_data(i,:),T_data(i,:),alt_grid,'linear','extrap');
        T_out(ind,:)=repmat(tmp_T,sum(ind),1);
        
    else
        continue
    end
    
end



end



