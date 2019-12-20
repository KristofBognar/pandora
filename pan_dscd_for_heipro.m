% function [table_out,scan_times,start_times,end_times]=...
%             pan_dscd_for_heipro( table_in, savedir, out_type, split_scans )
%INSERT_90 Insert 90 deg measurements into QDOAS output files for Pandora data
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
%%   
% INPUT:    table_in: MAX-DOAS dSCD file (straight from QDOAS), imported as a
%               matlab table
%           savedir: directory where output files will be saved
%           out_type: 'day', 'week', or 'month' for daily, weekly or monthly files
%           split_scans: 0 to write all scans into one file
%                        1 to separate short and long scans into different files
%                        2 to do 1 and 2 at the same time
%                        
%
% OUTPUT:   daily/weekly/monthly dSCD files with dummy 90 deg lines inserted
%           table_out: full dSCD table with dummy 90 deg lines inserted
%               (includes scans across file boundaries that were not written out)
%           scan_times: scan lengths (90 to next 90)
%
%
% Kristof Bognar, August 2019

%%%%%% Make sure QDOAS files are read by read_maxdoas_table_pandora.m,
%%%%%% there are some inportent filters/format changes in there

%% setup
% find first retrieval column (measurement info columns come first; assume
% that first retrieval window is O3)
% general info columns might change as people add columns with extra info
retr_start_col=find(strcmp(table_in.Properties.VariableNames,'O3_VisRMS'));

% make sure format of first few columns is the same (if not, change section
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


% remove datetime column (recalculated at the end)
try
    table_in.DateTime=[];
catch
    error('Check if files were processed the same way as in read_maxdoas_table_pandora.m (DateTime column missing)')
end
% gap tolerance
% max accepted time difference between consecutive measurements -- if time
% is greater, it's considered a gap
gap_tolerance=30/(24*60); % 30 minutes, converted to days

% round angles to nearest integer (done by read_maxdoas_table_pandora.m)
% table_in.Elevviewingangle=round(table_in.Elevviewingangle);
% table_in(table_in.NO2_VisSlColno2>1e20,:)=[];


%% format date/time fields properly
table_in.Timehhmmss.Format='HH:mm:ss';
table_in.DateDDMMYYYY.Format='dd/MM/yyyy';

table_in.Timehhmmss.Year=table_in.DateDDMMYYYY.Year;
table_in.Timehhmmss.Month=table_in.DateDDMMYYYY.Month;
table_in.Timehhmmss.Day=table_in.DateDDMMYYYY.Day;

%% get rid of missing times
ind_nat=isnat(table_in.Timehhmmss);
ind_ok=~isnat(table_in.Timehhmmss);

if sum(ind_nat)~=0

    table_in.Timehhmmss(ind_nat)=...
             ft_to_date(table_in.Fractionalday(ind_nat)-1,table_in.Year(ind_nat));

    if sum(isnat(table_in.Timehhmmss)), error('Could not fill all NaTs'); end
end

scan_times=table;

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

% make sure there's enough measurements at the end of the day to
% interpolate dummy times
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
        if table_in.Elevviewingangle(ind90(j))<30 && table_in.Elevviewingangle(ind90(j))<30
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
% any gap left that's longer than 30 min is a scan with no zenith measurements
% some 'zenith' measurements have elev angles of <80, not used as reference in QDOAS
tmp=find(scan_times>minutes(30));

if ~isempty(tmp)
    disp('')
    
    % save indices of bas scans
    bad_ranges=[];
    for i=1:length(tmp)
        bad_ranges(i,:)=[ind(tmp(i)),ind(tmp(i)+1)];
    end
    
    % check off angle zenith measurements
    
    elevs=unique(table_in.Elevviewingangle)';
    bad_elevs=setdiff(elevs,elevs_expected);
    
    disp('Start/stop indices of bad scans:')
    disp(bad_ranges)
    disp('Indices of bad elevation angles:')
    
    for i=1:length(bad_elevs)
        disp(['elev: ' num2str(bad_elevs(i))])
        disp(find(table_in.Elevviewingangle==bad_elevs(i)))
        % check if all bad elev angles are within the bad scans?
    end
    
    % delete bad scans
    for i=1:size(bad_ranges,1)
        % make sure a 90deg measurement is left at each end
        if table_in.Elevviewingangle(bad_ranges(i,1)-1)~=90
            bad_ranges(i,1)=bad_ranges(i,1)+1;
        end
        if table_in.Elevviewingangle(bad_ranges(i,2)+1)~=90
            bad_ranges(i,2)=bad_ranges(i,2)-1;
        end
        
        table_in(bad_ranges(i,1):bad_ranges(i,2),:)=[];
    end
    
    disp('')
    disp('Deleted bad scans')
    disp('')
    
    % check for bad elevs again
    elevs=unique(table_in.Elevviewingangle)';
    bad_elevs=setdiff(elevs,elevs_expected);
    
    if ~isempty(bad_elevs)
        disp('Bad elevation angles left:')
        for i=1:length(bad_elevs)
            disp(['elev: ' num2str(bad_elevs(i))])
            disp(find(table_in.Elevviewingangle==bad_elevs(i)))
        end
    end
    
    % get scan lengths again

    % all inserted 90deg dummies
    ind=find(table_in.Elevviewingangle==90);

    % time diffs between 90 deg lines
    scan_times=table_in.Timehhmmss(ind(2:end))-table_in.Timehhmmss(ind(1:end-1));

    % remove gaps (consecutive 90 deg lines, no scans in between)
    remove=find(ind(2:end)-ind(1:end-1)==1);
    scan_times(remove)=[];
    ind(remove)=[];
    
else % no scans >30 min, check for bad elevs anyway
    
    elevs=unique(table_in.Elevviewingangle)';
    bad_elevs=setdiff(elevs,elevs_expected);
    
    if ~isempty(bad_elevs)
        disp('Bad elevation angles (all scans <30 min):')
        for i=1:length(bad_elevs)
            disp(['elev: ' num2str(bad_elevs(i))])
            disp(find(table_in.Elevviewingangle==bad_elevs(i)))
            % check if all bad elev angles are within the bad scans?
        end
    end
    
end


table_out=table_in;


%% separate long and short scans, if required
% short scans:       30,   15,         2,1
% long scans:  50,40,30,20,15,10,8,5,3,2,1

% add temporary column to indicate scan length 
% long: idexed as 1
% short: indexed as -1
% 90 deg: indexed as 0 (included in both)
table_in.longscan=zeros(size(table_out.Elevviewingangle));
table_in.shortscan=zeros(size(table_out.Elevviewingangle));

if split_scans
    
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

end

%% write files

if ~exist(savedir,'dir'), mkdir(savedir); end

switch out_type
    case 'day' % get DOY (fractionalday in table is not necessarily correct...)
        tmp=day(table_in.DateTime,'dayofyear');
    case 'week' % get week number (Matlab starts the week on Sunday...)
        tmp=week(table_in.DateTime);
    case 'month' % get month number
        tmp=month(table_in.DateTime);
end
     
% loop over days/weeks/months
for i=unique(tmp)'

    % select data from given day
    to_write=table_in(tmp==i,:);

    % discard partial scans before first and after last 90deg dummy
    % there are a very large number of scans even for each day, no point
    % trying to save some of the scans that might be almost complete
    ind=find(to_write.Elevviewingangle==90);
    % use 90 immediately before/after scans (discard double 90s)
    if ind(2)-ind(1)==1, ind(1)=[]; end
    if ind(end)-ind(end-1)==1, ind(end)=[]; end

    to_write=to_write(ind(1):ind(end),:);

    % write files
    if split_scans
        
        %%% long scans
        to_write_tmp=to_write(to_write.longscan==1,:);
        
        yyyy=num2str(year(to_write_tmp.DateDDMMYYYY(1)));
        mm=num2str(month(to_write_tmp.DateDDMMYYYY(1)));
        dd=num2str(day(to_write_tmp.DateDDMMYYYY(1)));

        if length(mm)==1, mm=['0' mm]; end
        if length(dd)==1, dd=['0' dd]; end

        switch out_type
            case 'day' 
                fname=[savedir 'DSCD_' yyyy '_' mm '_' dd '_long.dat'];
            case 'week' 
                fname=[savedir 'DSCD_' yyyy '_week_' num2str(i) '_long.dat'];
            case 'month'
                fname=[savedir 'DSCD_' yyyy '_' mm '_long.dat'];
        end

        to_write_tmp.longscan=[];
        to_write_tmp.shortscan=[];
        writetable(to_write_tmp,fname,'Delimiter',',');
        
        %%% short scans
        to_write_tmp=to_write(to_write.shortscan==1,:);
        
        yyyy=num2str(year(to_write_tmp.DateDDMMYYYY(1)));
        mm=num2str(month(to_write_tmp.DateDDMMYYYY(1)));
        dd=num2str(day(to_write_tmp.DateDDMMYYYY(1)));

        if length(mm)==1, mm=['0' mm]; end
        if length(dd)==1, dd=['0' dd]; end

        switch out_type
            case 'day' 
                fname=[savedir 'DSCD_' yyyy '_' mm '_' dd '_short.dat'];
            case 'week' 
                fname=[savedir 'DSCD_' yyyy '_week_' num2str(i) '_short.dat'];
            case 'month'
                fname=[savedir 'DSCD_' yyyy '_' mm '_short.dat'];
        end

        to_write_tmp.longscan=[];
        to_write_tmp.shortscan=[];
        writetable(to_write_tmp,fname,'Delimiter',',');
        
    end
    
    if split_scans~=1
        % write all scans in one file
        
        yyyy=num2str(year(to_write.DateDDMMYYYY(1)));
        mm=num2str(month(to_write.DateDDMMYYYY(1)));
        dd=num2str(day(to_write.DateDDMMYYYY(1)));

        if length(mm)==1, mm=['0' mm]; end
        if length(dd)==1, dd=['0' dd]; end

        switch out_type
            case 'day' 
                fname=[savedir 'DSCD_' yyyy '_' mm '_' dd '.dat'];
            case 'week' 
                fname=[savedir 'DSCD_' yyyy '_week_' num2str(i) '.dat'];
            case 'month'
                fname=[savedir 'DSCD_' yyyy '_' mm '.dat'];
        end

        to_write.longscan=[];
        to_write.shortscan=[];
        writetable(to_write,fname,'Delimiter',',');

    end
        
end
        
% % % end