function [ P_out, T_out, alt_grid_out ] = get_ERA5_PT(dates_in_datetime, loc_in, alt_grid)
%GET_ERA5_PT get ERA5 pressure and temperature profiles on set altitude grid
%   
% ERA5 data is in pressure coordinates, interpolate to altitude grid
% Code assumes fixed file format for data
% One function call handles one year only
%
% Files might contain more than one profile per day, those are then
% averaged, and the mean daily profile is returned
%
% INPUT 
%       ERA5 data files: Using files with yearly data for each location.
%       The code finds the closest coordinate (only a few points per file).
%       The vertical grid is all pressure levels up to (and including)
%       10hPa, and UTC times of 11:00, 14:00, 17:00, 20:00, and 23:00 are
%       included (to get a decent daily average)
%       dates_in_datetime: datetime array of measurement times (individual
%                 measurements, mean time of sequence, etc; anything works)
%       alt_grid (optional): specify altitude levels for P,T interpolation.
%                 Default is the (approximate) mean geopotential height of each
%                 pressure level in the 2018 ERA5 files
%                 UNITS must be meters!!
%
% OUTPUT: P,T output is daily average (DAYTIME ONLY: mean of 11:00,14:00,17:00,20:00,23:00 UTC)
%         values are repeated for input times on the same day, so the dimensions match
%       P_out: Pressure in hPa, interpolated linearly
%       T_out: Temperature in K, interpolated linearly
%       alt_grid_out: altitude in meters corresponding to P_out and T_out
%                     (returns default if alt_grid is not passed in)
%
%@Kristof Bognar, 2020


%% setup

% find ERA5 files
era5_loc='/home/kristof/work/models/ERA5/';

if nargin<3,
   
    % default altitude grid, based on mean height of each pressure level in
    % 2018-2019 Downsview ERA5 data
    load([era5_loc 'alt_grid.mat']) % in m
        
end

if strcmpi(loc_in,'downsview')
    
    % ECCC Downsview   
    location.lat=43.7810; 
    location.lon=-79.4680;
    loc_str='Toronto';
    
elseif strcmpi(loc_in,'uoft')
    error('Set coords')
end

% year of interest (data should have one year only)
yr=unique(year(dates_in_datetime));
if length(yr)>1, error('Call GET_ERA5_PT one year at a time'), end

% current year
yr_now=year(datetime(now,'convertfrom','datenum'));

% return altitude grid
alt_grid_out=alt_grid;

% otput arrays
P_out=NaN(length(dates_in_datetime),length(alt_grid));
T_out=NaN(length(dates_in_datetime),length(alt_grid));

% get days in yyyymmdd format
% datenum_in=str2num(datestr(dateshift(dates_in_datetime, 'start', 'day'),'yyyymmdd'));
datenum_in=floor(datenum(dates_in_datetime));
days_in=unique(datenum_in);

%% read file

tmp=dir([era5_loc 'ERA5_' loc_str '*' num2str(yr) '*']);
fname={tmp.name};
fname=fname{1};

lon=double(ncread([era5_loc fname],'longitude'));
lat=double(ncread([era5_loc fname],'latitude'));
p_grid=double(ncread([era5_loc fname],'level')); % units: hPa
times=double(ncread([era5_loc fname],'time'));
z_data=double(ncread([era5_loc fname],'z'));
z_data=z_data./9.80665; % calculate geopotential height (meters)
t_data=double(ncread([era5_loc fname],'t')); % units: K

% convert time to datetime
% file time is hours since 1900-01-01 00:00:00.0
% fines should contain one profile per day, at 17:00 utc
pivot_time = datenum('1900-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');

% convert to matlab datenum
era5_datenum=floor(times/24 + pivot_time);
% era5_datetime=datetime(times/24 + pivot_time,'convertfrom','datenum');
% era5_dates=str2num(datestr(dateshift(era5_datetime, 'start', 'day'),'yyyymmdd'));

%% find nearest point

% coordinates
[~,ind_lon]=min(abs(lon-location.lon));
[~,ind_lat]=min(abs(lat-location.lat));

% check for ERA5T (near real time) preliminary data
% expver variable exist in file, 1 for ERA5, 5 for ERA5T
try % expver exists: recent data
    
    tmp=double(ncread([era5_loc fname],'expver'));
    %if yr~=yr_now, warning('Near real time data found in file'); end
    
    % reduce array sizes
    % both datasets are the same size, one dataset is all NaN where the
    % other has data
    z_data1=squeeze(z_data(ind_lon,ind_lat,:,1,:));
    t_data1=squeeze(t_data(ind_lon,ind_lat,:,1,:));
   
    z_data2=squeeze(z_data(ind_lon,ind_lat,:,2,:));
    t_data2=squeeze(t_data(ind_lon,ind_lat,:,2,:));

    % replace NaN values with data from second array
    z_data1(:,isnan(z_data1(1,:)))=z_data2(:,isnan(z_data1(1,:)));
    t_data1(:,isnan(t_data1(1,:)))=t_data2(:,isnan(t_data1(1,:)));
    
    z_data=z_data1;
    t_data=t_data1;
    
catch % expver doesn't exist, consolidated data

    % reduce array sizes
    z_data=squeeze(z_data(ind_lon,ind_lat,:,:));
    t_data=squeeze(t_data(ind_lon,ind_lat,:,:));

end

%% interpolate data
for i=1:length(days_in)

    % current chunk of data
    ind=(datenum_in==days_in(i));
    
    ind_era5=find(era5_datenum==days_in(i));
        
    if length(ind_era5)>1 % average all times for a given day
    
        % get pressure
        tmp_P=interp1(mean(z_data(:,ind_era5),2),p_grid,alt_grid,'linear','extrap');
        P_out(ind,:)=repmat(tmp_P,sum(ind),1);

        % get temperature
        tmp_T=interp1(mean(z_data(:,ind_era5),2),mean(t_data(:,ind_era5),2),alt_grid,'linear','extrap');
        T_out(ind,:)=repmat(tmp_T,sum(ind),1);
    
    elseif length(ind_era5)==1 % only one time per day
        
        % get pressure
        tmp_P=interp1(z_data(:,ind_era5),p_grid,alt_grid,'linear','extrap');
        P_out(ind,:)=repmat(tmp_P,sum(ind),1);

        % get temperature
        tmp_T=interp1(z_data(:,ind_era5),t_data(:,ind_era5),alt_grid,'linear','extrap');
        T_out(ind,:)=repmat(tmp_T,sum(ind),1);
        
    else
        error(['Date ' num2str(days_in(i)) ' not found in ERA5 file'])
    end
    
end

end

