 function [pan_maxdoas_processed,scan_times]=...
            reformat_pandora_for_retrievals( p_num, uvvis, yr_in, savedir, version )
% reformat rqw QDOAS output into format required for retrieval input file
% creation
%
% code differs drastically for fifferent data versions
%
% version<5 (2 and 4 really): QDOAS was not saving failed spectra
%           references were not saved
%           Necessary to insert 90deg measurements in between scans
%
% version>=5: QDOAS can save references as well (!)
%           (Kristof: realized this after 5 years...)
%           Inserting 90deg redundant, some formatting and checks still
%           necessary
%



%%%%%%%%%%%%%%% original description for version <5 %%%%%%%%%%%%%%%%%
%
%Insert 90 deg measurements into QDOAS output files for Pandora data
%
% Reads Pandora MAX-DOAS table and breaks it up by day/week/month
%
% Pandora scans go from 90 down to lowest elev, then back up to 90 again
% Lowest elev (usually 1deg) is not repeated
%
% 90 deg placeholders are inseted using interpolated properties from two 30
% (or 50) deg measurements at the start and end of the scans
%       - first and last line are extrapolated
%
% Pandora measurements have two scanning modes:
%           30,     15,              2, 1, and up (tagged as shortscan==1)
%   50, 40, 30, 20, 15, 10, 8, 5, 3, 2, 1, and up (tagged as longscan==1)
%   90 deg dummies in between short and long scans are tagged as belonging
%   to both, to make separating scan types easier
%
%%   
% INPUT:    p_num: pandora instrument number
%           uvvis: 'uv' for UV data, 'vis' for visible data
%           year: year to process
%           savedir (optional): if provided, pan_maxdoas_processed and scan_times are saved here
%
% OUTPUT:   pan_maxdoas_processed: full dSCD table with dummy 90 deg lines inserted
%               (includes scans across file boundaries that were not written out)
%           scan_times: scan lengths (90 to next 90)
%
%
% Kristof Bognar, August 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~any([2,4,5]==version), error('Incorrect version number'); end

disp('Make sure QDOAS files are read by merge_pandora_dscds.m,')
disp(' ')

%% setup
% load QDOAS data (has to be saved by merge_pandora_dscds.m)
filename=['Pandora_' num2str(p_num) '_' uvvis '_' num2str(yr_in) '.mat'];

filedir=['/home/kristof/work/PANDORA/profiling/QDOAS_output/spec_v' num2str(version) '/'];
disp(['Using QDOAS_output/spec_v' num2str(version) ' as source'])
disp(' ')

try
    load([filedir filename]);
    table_in=data;
    clearvars data
catch
    error('Read QDOAS files using merge_pandora_dscds,m')
end


% find fit windows in the file (each window has one RMS output field)
tmp=table_in.Properties.VariableNames;
cols_in_file=tmp(find_in_cell(tmp,'RMS'));

% convert to window name
for i=1:length(cols_in_file)
    cols_in_file{i}=cols_in_file{i}(1:end-3);
end

disp('Analysis windows found in file:')
disp(cols_in_file)


% find first retrieval column (measurement info columns come first; assume
% that first column is the RMS)
% general info columns might change as people add columns with extra info
switch uvvis
    case 'vis'
        retr_start_col=find(strcmp(table_in.Properties.VariableNames,[cols_in_file{1} 'RMS']));
    case 'uv'
        error('specify format')
end

% version<5: make sure format of first few columns is the same (if not, change section
% that adds dummy  90deg lines in middle of file)
%%% Expected columns:
% Spec No
% Year
% Fractional day
% Fractional time
% SZA
% Solar Azimuth Angle
% Elev. viewing angle
% Azim. viewing angle
% Date (DD/MM/YYYY)
% Time (hh:mm:ss)
% Tint
% Total Experiment Time (sec)

% expected possible elevation angles in file
elevs_expected=[1,2,3,5,8,10,15,20,30,40,50,90];


% gap tolerance
% max accepted time difference between consecutive measurements -- if time
% is greater, it's considered a gap
%
% there are frequent 5-15 min gaps when direct sun measurements are made,
% 90deg dummies are inserted in the middle of some of these short gaps,
% slightliy far from the next scan (should be OK for HEIPRO, doesn't matter
% for MAPA)
if version<5
    gap_tolerance=10/(24*60); % 10 minutes, converted to days
%     disp(['Using ' num2str(gap_tolerance*(24*60)) ' min as cutoff for data gaps'])
%     disp(' ')
end

% round angles to nearest integer (done by read_maxdoas_table_pandora.m)
% table_in.Elevviewingangle=round(table_in.Elevviewingangle);
% table_in(table_in.NO2_VisSlColno2>1e20,:)=[];

%% remove scans with uncommon azimuth angles
% necessary for both versions

% consider only elevs <=50, anithing above is either used as zenith reference
% or ignored for analysis
% az angles for senith reference measurements don't matter much
az_vals=table_in.Azimviewingangle(table_in.Elevviewingangle<=50);

az_list=unique(az_vals);
az_count=histc(az_vals,unique(az_vals));

if length(az_list)>1
    
    disp('Multiple azimuth angles present (angle + n.o. repetitions):') 
    disp([az_list, az_count])

    disp('Enter cutoff for number of measurements with given AZ angle')
    disp('Angles that repeat fewer times than given are deleted (default:25)')
    
    tmp=input('','s');

    if ~isempty(tmp)
        az_cutoff=str2double(tmp);
    else
        az_cutoff=25;
    end

    for i=1:length(az_list)
        if az_count(i)<az_cutoff
            table_in(table_in.Azimviewingangle==az_list(i) & table_in.Elevviewingangle<=50,:)=[]; 
        end
    end

    % check results
    az_vals=table_in.Azimviewingangle(table_in.Elevviewingangle<=50);
    az_list=unique(az_vals);
    az_count=histc(az_vals,unique(az_vals));
    disp('Remaining azimuth angles:') 
    disp([az_list, az_count])

end

%% remove non-zenith elev angles
table_in(table_in.Elevviewingangle<80 & table_in.Elevviewingangle>50,:)=[];


%% format date/time fields properly, and get rid of missing times
table_in.Timehhmmss.Format='HH:mm:ss';
table_in.DateDDMMYYYY.Format='dd/MM/yyyy';

table_in.Timehhmmss.Year=table_in.DateDDMMYYYY.Year;
table_in.Timehhmmss.Month=table_in.DateDDMMYYYY.Month;
table_in.Timehhmmss.Day=table_in.DateDDMMYYYY.Day;

% get rid of missing times
ind_nat=isnat(table_in.Timehhmmss);
ind_ok=~isnat(table_in.Timehhmmss);

if sum(ind_nat)~=0

    table_in.Timehhmmss(ind_nat)=...
             ft_to_date(table_in.Fractionalday(ind_nat)-1,table_in.Year(ind_nat));

    if sum(isnat(table_in.Timehhmmss)), error('Could not fill all NaTs'); end
end

scan_times=table;


%% add smoothness filter
% do it before any formatting so bad data (e.g. days with high RMS) can be removed
% filter values will be messed up for 90deg dummies, but those lines are
% discarded for MAPA anyway (and extra columns don't impact HEIPRO processing)

% set flter parameters
switch uvvis
    case 'vis'
        rms_lim=0.003; % rms limit (could specify different values for each column)
    case 'uv'
        error('specify RMS limits')
end

smooth_window=3; % smoothing window in hours, also min interval of daily data to do smoothing
N_min=10; % min number of data points to do smoothing

% run smoothing function 
[ind_bad_ci,ind_bad_o4,ind_bad_rms]=pandora_dSCD_filter(table_in,cols_in_file,...
                                                        rms_lim,smooth_window,N_min);

% add cloud flag (possible values: 0, 1, 2, 3)
table_in.cloud_flag=zeros(size(table_in,1),1); % good data: flag is 0
table_in.cloud_flag(ind_bad_ci & ~ind_bad_o4)=1; % only ci smoothness off
table_in.cloud_flag(ind_bad_o4 & ~ind_bad_ci)=2; % only O4 smoothness off
table_in.cloud_flag(ind_bad_ci & ind_bad_o4)=3; % both ci and O4 smoothness off

% Mask data with high RMS
for i=1:length(cols_in_file)
    % find data for given analysis window 
    col_ind=find_in_cell(table_in.Properties.VariableNames,cols_in_file{i});
    table_in{ind_bad_rms(:,i),col_ind}=NaN;
end

% add color index column for later
table_in.ci=table_in.Fluxes330./table_in.Fluxes440;

% remove datetime column (was needed for pandora_dSCD_filter, will be recalculated at the end)
table_in.DateTime=[];

%% format data
if version<5
    %%%%
    %%%% insert 90 deg dummy lines into entire dataset (likely one year)
    %%%%

    %% first 90
    % just extrapolate from first 2 measurements, no need for serious checks

    % replicate first line, and modify time so it passes as 90 deg dummy
    table_in=[table_in(1,:); table_in];

    table_in.SpecNo(1)=table_in.SpecNo(1)-1;
    table_in.Fractionalday(1)=table_in.Fractionalday(1)-diff(table_in.Fractionalday(2:3));
    table_in.Fractionaltime(1)=table_in.Fractionaltime(1)-diff(table_in.Fractionaltime(2:3));
    table_in.SZA(1)=table_in.SZA(1)-diff(table_in.SZA(2:3));
    table_in.SolarAzimuthAngle(1)=table_in.SolarAzimuthAngle(1)-diff(table_in.SolarAzimuthAngle(2:3));

    table_in.Elevviewingangle(1)=90;

    table_in.Timehhmmss(1)=table_in.Timehhmmss(1)-diff(table_in.Timehhmmss(2:3));

    table_in(1,retr_start_col:end)=num2cell(zeros(size(table_in(1,retr_start_col:end))));

    % check if time slipped back to prev. day
    if hour(table_in.Timehhmmss(1))==23 && hour(table_in.Timehhmmss(2))==0

        table_in.Timehhmmss(1).Hour=0;
        table_in.Timehhmmss(1).Minute=0;
        table_in.Timehhmmss(1).Second=2;
        table_in.Timehhmmss(1).Day=table_in.DateDDMMYYYY(1).Day;
        table_in.Timehhmmss(1).Month=table_in.DateDDMMYYYY(1).Month;
        table_in.Fractionaltime(1)=0.0005556;

    end

    %% mid-file 90

    % find places where spectrum number skips 1 (or more), that's where 90 deg
    % placeholders should go
    ind90=find(table_in.SpecNo(1:end-1)-table_in.SpecNo(2:end)~=-1);
    ind90=ind90+1;

    % handle data gaps (nighttime, other modes, data breaks)
    % need 90deg dummy on both ends, instead of the middle
    ind90_gap=find(table_in.Fractionalday(1:end-1)-table_in.Fractionalday(2:end) < -gap_tolerance);
    ind90_gap=ind90_gap+1;

    % make sure last gap (if any) has at least two measurements after it
    % (not too important here, left ower from daily processing of PGBS data)
    if ~isempty(ind90_gap) && (size(table_in,1)-ind90_gap(end)<2) 

        % remove corresponding 90 index, so code doesn't attempt to add 90deg dummy
        ind90(ind90==ind90_gap(end))=[];
        table_in(ind90_gap(end):end,:)=[];
        ind90_gap(end)=[];

        % check the same for the now truncated file
        search=1;
        while search
            if ~isempty(ind90_gap) && (size(table_in,1)-ind90_gap(end)<2) 
                ind90(ind90==ind90_gap(end))=[];
                table_in(ind90_gap(end):end,:)=[];
                ind90_gap(end)=[];
            else
                search=0;
            end
        end    
    end

    % make sure there's enough measurements at the end of the dataset to interpolate dummy times
    if length(ind90)>1 && ind90(end)==size(table_in,1) 
        % ind90 indicates start of scan
        % only one measurement in last scan: delete measurement and index
        ind90(end)=[];
        table_in(end,:)=[]; 
    end


    % insert 90 deg dummies
    n=0;
    for j=1:length(ind90)
        % display progress info
        disp_str=['Processing data from: ' datestr(table_in.Timehhmmss(ind90(j)))];
        % stuff to delete last line and reprint updated message
        fprintf(repmat('\b',1,n));
        fprintf(disp_str);
        n=numel(disp_str);    


        % check if there's a gap
        if ismember(ind90(j),ind90_gap)

            % if yes, insert two dummies, on each end of the gap
            table_in=[table_in(1:ind90(j)-1,:);...
                       table_in(ind90(j)-1,:);...
                       table_in(ind90(j),:);...
                       table_in(ind90(j):end,:)];

            gap_start=ind90(j);
            gap_end=ind90(j)+1;

            %%% add dummy at start of the gap
            table_in.SpecNo(gap_start)=table_in.SpecNo(gap_start)+1;
            table_in.Fractionalday(gap_start)=table_in.Fractionalday(gap_start)+diff(table_in.Fractionalday(gap_start-2:gap_start-1));
            table_in.Fractionaltime(gap_start)=table_in.Fractionaltime(gap_start)+diff(table_in.Fractionaltime(gap_start-2:gap_start-1));
            table_in.SZA(gap_start)=table_in.SZA(gap_start)+diff(table_in.SZA(gap_start-2:gap_start-1));
            table_in.SolarAzimuthAngle(gap_start)=table_in.SolarAzimuthAngle(gap_start)+diff(table_in.SolarAzimuthAngle(gap_start-2:gap_start-1));

            table_in.Elevviewingangle(gap_start)=90;

            table_in.Timehhmmss(gap_start)=table_in.Timehhmmss(gap_start)+diff(table_in.Timehhmmss(gap_start-2:gap_start-1));

            % check if time shifted to next day, adjust date (relevant for PEARL data only)
            if hour(table_in.Timehhmmss(gap_start))==0 && hour(table_in.Timehhmmss(gap_start-1))==23
                table_in.DateDDMMYYYY(gap_start)=table_in.DateDDMMYYYY(gap_start)+days(1);
                table_in.Fractionaltime(gap_start)=table_in.Fractionaltime(gap_start)-24;
            end

            table_in(gap_start,retr_start_col:end)=num2cell(zeros(size(table_in(1,retr_start_col:end))));

            %%% add dummy at the end of the gap
            table_in.SpecNo(gap_end)=table_in.SpecNo(gap_end)-1;
            table_in.Fractionalday(gap_end)=table_in.Fractionalday(gap_end)-diff(table_in.Fractionalday(gap_end+1:gap_end+2));
            table_in.Fractionaltime(gap_end)=table_in.Fractionaltime(gap_end)-diff(table_in.Fractionaltime(gap_end+1:gap_end+2));
            table_in.SZA(gap_end)=table_in.SZA(gap_end)-diff(table_in.SZA(gap_end+1:gap_end+2));
            table_in.SolarAzimuthAngle(gap_end)=table_in.SolarAzimuthAngle(gap_end)-diff(table_in.SolarAzimuthAngle(gap_end+1:gap_end+2));

            table_in.Elevviewingangle(gap_end)=90;

            table_in.Timehhmmss(gap_end)=table_in.Timehhmmss(gap_end)-diff(table_in.Timehhmmss(gap_end+1:gap_end+2));

            % check if time shifted to prev. day, adjust date (relevant for PEARL data only)
            if hour(table_in.Timehhmmss(gap_end))==23 && hour(table_in.Timehhmmss(gap_end+1))==0
                table_in.DateDDMMYYYY(gap_end)=table_in.DateDDMMYYYY(gap_end)-days(1);
                table_in.Fractionaltime(gap_end)=table_in.Fractionaltime(gap_end)+24;
            end

            table_in(gap_end,retr_start_col:end)=num2cell(zeros(size(table_in(1,retr_start_col:end))));

            % account for added rows
            ind90=ind90+2;
            ind90_gap=ind90_gap+2;

        else % if not a gap, interpolate between two scans

            % check if we're in the middle of a scan -- nothing should be inserted there
            % scans end with 30, 40, or 50deg meas.: one scan might be
            % incomplete, but if neither end of the gap is 30 or above, then
            % it's likely a missed elevation angle
            if table_in.Elevviewingangle(ind90(j)-1)<30 && table_in.Elevviewingangle(ind90(j))<30
                continue
            end

            % insert 90deg dummy
            table_in=[table_in(1:ind90(j)-1,:);...
                       table_in(ind90(j),:);...
                       table_in(ind90(j):end,:)];

            % average spectra details, replace necessary entries below
            table_in(ind90(j),1:8)=num2cell(mean(table_in{ind90(j)-1:2:ind90(j)+1,1:8}));

            % fix spectrum number and solar azimuth
            table_in.SpecNo(ind90(j))=round(table_in.SpecNo(ind90(j)));        
            if table_in.SpecNo(ind90(j)) < table_in.SpecNo(ind90(j)-1),
                table_in.SpecNo(ind90(j)) = table_in.SpecNo(ind90(j)-1)+1;
            end

            % calculate date and time
            table_in.Timehhmmss(ind90(j))=mean([table_in.Timehhmmss(ind90(j)-1),table_in.Timehhmmss(ind90(j)+1)]);

            table_in.DateDDMMYYYY(ind90(j))=table_in.Timehhmmss(ind90(j));
            table_in.DateDDMMYYYY(ind90(j)).Hour=0;
            table_in.DateDDMMYYYY(ind90(j)).Minute=0;
            table_in.DateDDMMYYYY(ind90(j)).Second=0;

            table_in.Fractionaltime(ind90(j))=(table_in.Fractionalday(ind90(j))-floor(table_in.Fractionalday(ind90(j))))*24;

            % recalculate solar azimuth when measuring across midnight (PEARL only)
            if table_in.SolarAzimuthAngle(ind90(j)-1)>0 && table_in.SolarAzimuthAngle(ind90(j)+1)<0
                table_in.SolarAzimuthAngle(ind90(j))=...
                    mean([table_in.SolarAzimuthAngle(ind90(j)-1),360+table_in.SolarAzimuthAngle(ind90(j)+1)]);
                if table_in.SolarAzimuthAngle(ind90(j))>180
                    table_in.SolarAzimuthAngle(ind90(j))=table_in.SolarAzimuthAngle(ind90(j))-360;
                end
            end

            % add 90 for elev viewing angle
            table_in.Elevviewingangle(ind90(j))=90;
            % set dSCD data to 0
            table_in(ind90(j),retr_start_col:end)=num2cell(zeros(size(table_in(ind90(j),retr_start_col:end))));

            % account for added row
            ind90=ind90+1;
            ind90_gap=ind90_gap+1;

        end
    end

    fprintf('\n');

    %% last 90

    % replicate last line, and modify time so it passes as 90 deg dummy
    table_in=[table_in; table_in(end,:)];

    table_in.SpecNo(end)=table_in.SpecNo(end)+1;
    table_in.Fractionalday(end)=table_in.Fractionalday(end)+diff(table_in.Fractionalday(end-2:end-1));
    table_in.Fractionaltime(end)=table_in.Fractionaltime(end)+diff(table_in.Fractionaltime(end-2:end-1));
    table_in.SZA(end)=table_in.SZA(end)+diff(table_in.SZA(end-2:end-1));
    table_in.SolarAzimuthAngle(end)=table_in.SolarAzimuthAngle(end)+diff(table_in.SolarAzimuthAngle(end-2:end-1));

    table_in.Elevviewingangle(end)=90;

    table_in.Timehhmmss(end)=table_in.Timehhmmss(end)+diff(table_in.Timehhmmss(end-2:end-1));

    table_in(end,retr_start_col:end)=num2cell(zeros(size(table_in(1,retr_start_col:end))));


    % check if time slipped to next day
    if hour(table_in.Timehhmmss(end)) == 0

        table_in.Timehhmmss(end).Hour=23;
        table_in.Timehhmmss(end).Minute=59;
        table_in.Timehhmmss(end).Second=58;
        table_in.Timehhmmss(end).Day=table_in.DateDDMMYYYY(end).Day;
        table_in.Timehhmmss(end).Month=table_in.DateDDMMYYYY(end).Month;
        table_in.Fractionaltime(end)=23.9994444;

    end
    
else
    
    %%%%
    %%%% edit 90deg lines to match expected format
    %%%%
    
    % replace all reference elevs with 90
    table_in.Elevviewingangle(table_in.Elevviewingangle>80)=90;
    
    % if more than 2 consecutive 90 lines: leave only first and last
    % QDOAS outputs all 90deg measurements at beginning and end of day
    ind=find(table_in.Elevviewingangle==90);
    tmp=find(ind(3:end)-ind(1:end-2)==2);
    tmp=tmp+1;
    
    remove=ind(tmp);
    
    table_in(remove,:)=[];
    
    % leave only 1 90deg line at start/end of file
    if isequal(table_in.Elevviewingangle(1:2),[90;90]), table_in(1,:)=[]; end
    if isequal(table_in.Elevviewingangle(end-1:end),[90;90]), table_in(end,:)=[]; end
    
    % replace 90deg line content with zeroes (leave fluxes intact)
    ind=(table_in.Elevviewingangle==90);    
    for i=1:length(cols_in_file)
        % find data for given analysis window 
        col_ind=find_in_cell(table_in.Properties.VariableNames,cols_in_file{i});
        table_in{ind,col_ind}=0;
    end    
    
    % replace az angles for zenith measurements with az of nearest scan
    table_in.Azimviewingangle(ind)=interp1(table_in.Fractionalday(~ind),...
                                           table_in.Azimviewingangle(~ind),...
                                           table_in.Fractionalday(ind),'nearest','extrap');
    
end

%% Recalculate datetime column

table_in.DateTime=table_in.DateDDMMYYYY+timeofday(table_in.Timehhmmss);
table_in.DateTime.Format='dd/MM/uuuu HH:mm:ss';

%% get scan lengths

% all inserted 90deg dummies
ind=find(table_in.Elevviewingangle==90);

% time diffs between 90 deg lines
scan_times=table_in.Timehhmmss(ind(2:end))-table_in.Timehhmmss(ind(1:end-1));

% remove gaps (consecutive 90 deg lines, no scans in between)
remove=find(ind(2:end)-ind(1:end-1)==1);
scan_times(remove)=[];
ind(remove)=[];

%% deal with exceptions

% any scan left that's longer than 20 min is probably a scan with no zenith
% measurements (based on 2018 data from 103 -- using 10min gaps for 90 deg
% dummies, 90 to 90 duration is 17 min)

% some 'zenith' measurements have elev angles of <80, not used as reference in QDOAS

tmp=find(scan_times>minutes(20));

if ~isempty(tmp)
    disp('')
    
    % save indices of bas scans
    bad_ranges=[];
    for i=1:length(tmp)
        bad_ranges(i,:)=[ind(tmp(i)),ind(tmp(i)+1)];
    end
    
    disp('Start/stop indices of bad scans:')
    disp(bad_ranges)
    
    % check for off-angle zenith measurements -- redundant, since they are
    % removed at the start
%     elevs=unique(table_in.Elevviewingangle)';
%     bad_elevs=setdiff(elevs,elevs_expected);
%     
%     if ~isempty(bad_elevs)
%         disp('Indices of bad elevation angles:')
% 
%         for i=1:length(bad_elevs)
%             disp(['elev: ' num2str(bad_elevs(i))])
%             disp(find(table_in.Elevviewingangle==bad_elevs(i)))
%             % check if all bad elev angles are within the bad scans?
%         end
%     end
    
    % delete bad scans
    bad_inds=[];
    for i=1:size(bad_ranges,1)
        % make sure a 90deg measurement is left at each end
        if table_in.Elevviewingangle(bad_ranges(i,1)-1)~=90
            bad_ranges(i,1)=bad_ranges(i,1)+1;
        end
        if table_in.Elevviewingangle(bad_ranges(i,2)+1)~=90
            bad_ranges(i,2)=bad_ranges(i,2)-1;
        end
        
        bad_inds=[bad_inds,bad_ranges(i,1):bad_ranges(i,2)];
    end
    
    table_in(bad_inds,:)=[];    
    
    disp('')
    disp('Deleted bad scans')
    disp('')
    
    % check for bad elevs again
%     elevs=unique(table_in.Elevviewingangle)';
%     bad_elevs=setdiff(elevs,elevs_expected);
%     
%     if ~isempty(bad_elevs)
%         disp('Bad elevation angles left:')
%         for i=1:length(bad_elevs)
%             disp(['elev: ' num2str(bad_elevs(i))])
%             disp(find(table_in.Elevviewingangle==bad_elevs(i)))
%         end
%     end
        
else % no scans >gap tolerance, check for bad elevs anyway
    
%     elevs=unique(table_in.Elevviewingangle)';
%     bad_elevs=setdiff(elevs,elevs_expected);
%     
%     if ~isempty(bad_elevs)
%         disp('Bad elevation angles (all scans <20 min):')
%         for i=1:length(bad_elevs)
%             disp(['elev: ' num2str(bad_elevs(i))])
%             disp(find(table_in.Elevviewingangle==bad_elevs(i)))
%             % check if all bad elev angles are within the bad scans?
%         end
%     end
    
end


%% separate long and short scans
% short scans:       30,   15,         2,1
% long scans:  50,40,30,20,15,10,8,5,3,2,1

% add temporary column to indicate scan length 
% long: idexed as 1
% short: indexed as -1
% 90 deg: indexed as 0 (included in both)
table_in.longscan=zeros(size(table_in.Elevviewingangle));
table_in.shortscan=zeros(size(table_in.Elevviewingangle));
    
% redo zenith rows index
ind90=find(table_in.Elevviewingangle==90);

% tag each measurement
longscan=[];
shortscan=[];

for i=1:length(ind90)-1

    % elevs in between 90 deg measurements
    curr_ind=ind90(i)+1:ind90(i+1)-1;
    tmp=table_in.Elevviewingangle(curr_ind);

    if any(ismember(tmp,[50,40,20,10,8,5,3])) % elevs that only appear in long scans
        longscan=[longscan;curr_ind'];
    else
        shortscan=[shortscan;curr_ind'];
    end
end

% include 90 deg before/after each scan
longscan=union(longscan,longscan+1);
longscan=union(longscan,longscan-1);
shortscan=union(shortscan,shortscan+1);
shortscan=union(shortscan,shortscan-1);

% save indices in data table
table_in.longscan(longscan)=1;
table_in.shortscan(shortscan)=1;

pan_maxdoas_processed=table_in;

%% save results

if ~strcmp(savedir(end),'/'), savedir=[savedir, '/']; end
if ~exist(savedir,'dir'), mkdir(savedir); end

save([savedir 'p' num2str(p_num) '_' uvvis '_' num2str(yr_in) '_maxdoas_processed.mat'],'pan_maxdoas_processed')

end

