function [ P_out, T_out ] = get_insitu_PT( dates_in_datetime, loc_in, out_type)
%GET_INSITU_PT get mean pressure and temperature from in situ measurements
% 
% INPUT 
%       met data files: from ECCC, see import_met function
%                       time is local; convert to UTC
%       dates_in_datetime: datetime array of measurement times (individual
%                 measurements, mean time of sequence, etc; anything works)
%       out_type: 'mean' for daily mean values (DAYTIME ONLY: mean of 11:00 to 23:00 UTC),
%                        values are repeated for input times on the same day
%                 'interp' for linearly interpolated values for each input time
%
% OUTPUT
%       P_out: Pressure in hPa
%       T_out: Temperature in K
%
%@Kristof Bognar, 2020

    %% setup
    met_dir='/home/kristof/work/PANDORA/in_situ/';
    
    if strcmpi(loc_in,'downsview')
    
        % ECCC Downsview   
        location.lat=43.7810; 
        location.lon=-79.4680;
        loc_str='Toronto';
        
        % load data
        try
            load([met_dir 'Toronto_North_34021_met_data_2017-2019.mat'])
        catch
            met_data=import_met(met_dir);
        end

    elseif strcmpi(loc_in,'uoft')
        error('Set coords')
    end
    
    % output variables
    P_out=NaN(size(dates_in_datetime));
    T_out=NaN(size(dates_in_datetime));

    
    if strcmpi(out_type,'mean')
        %% get daily averages (DAYTIME ONLY: mean of 11:00 to 23:00 UTC)

        % remove morning/evening hours that are not used in the average
        % times are already UTC
        met_data(met_data.DateTime.Hour<11,:)=[];
        met_data(met_data.DateTime.Hour>23,:)=[]; % 23:00 is the last entry in hourly data

        % unique days in input datetime array
        time_steps=day(dates_in_datetime,'dayofyear')+dates_in_datetime.Year*1000;    
        % unique dates in met data
        time_met=day(met_data.DateTime,'dayofyear')+met_data.DateTime.Year*1000;

        % loop over days/weeks/months
        for i=unique(time_steps)'

            ind=(time_steps==i);

            if ~isempty(ind)
                tmp=nanmean(met_data.PRES(time_met==i));
                P_out(ind)=repmat(tmp,sum(ind),1);

                tmp=nanmean(met_data.TEMP(time_met==i))+273.15;
                T_out(ind)=repmat(tmp,sum(ind),1);
            end

        end
        
    elseif strcmpi(out_type,'interp')
        %% interpolate P, T values to individual input dates
        
        P_out=interp1(met_data.DateTime,met_data.PRES,dates_in_datetime);
        T_out=interp1(met_data.DateTime,met_data.TEMP,dates_in_datetime);
        
        T_out=T_out+273.15;
        
    end
end

function met_data=import_met(met_dir)
    % Auto-generated by MATLAB on 2020/05/05 17:46:12

    %% Import the data
    
    [~, ~, raw] = xlsread([met_dir 'Toronto_North_34021_met_data_2017-2019.xlsx'],'Sheet1');
    raw = raw(6:end,:);
    raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
    cellVectors = raw(:,1);
    raw = raw(:,[2,3,4]);

    %% Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
    raw(R) = {NaN}; % Replace non-numeric cells

    %% Create output variable
    data = reshape([raw{:}],size(raw));

    %% Create table
    met_data = table;

    %% Allocate imported array to column variable names
    met_data.DateTime = cellVectors(:,1);

    % convert to datetime, need UTC
    % Original data is just EST! No daylight savings, UTC offset is 5h all year round   
    met_data.DateTime=datetime(met_data.DateTime,'InputFormat','dd/MM/yyyy HH:mm');
    met_data.DateTime=met_data.DateTime+hours(5);
    
    met_data.RH = data(:,1);
    met_data.TEMP = data(:,2);
    met_data.PRES = data(:,3);
    
    %% remove NaNs
    met_data=met_data(~any(ismissing(met_data),2),:);
    
    %% save as .mat file
    save([met_dir 'Toronto_North_34021_met_data_2017-2019.mat'],'met_data')
    
end