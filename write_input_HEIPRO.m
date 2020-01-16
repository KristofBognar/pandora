function write_input_HEIPRO( table_in, savedir, out_type, split_scans )
% Write Pandora dSCDs to input files used by HEIPRO
%
%%   
% INPUT:    table_in: MAX-DOAS data processed by reformat_pandora_for_retrievals.m
%                     Data is filtered, gaps are dealt with, and 90 deg
%                     dummies are inserted by that function
%           savedir: directory where output files will be saved
%           out_type: 'day', 'week', or 'month' for daily, weekly or monthly files
%           split_scans: 0 to write all scans into one file
%                        1 to separate short and long scans into different files
%                        2 to do 1 and 2 at the same time
%
% OUTPUT:   daily/weekly/monthly dSCD files with dummy 90 deg lines inserted
%
%
% Kristof Bognar, August 2019

%% write files

if ~strcmp(savedir(end),'/'), savedir=[savedir, '/']; end
if ~exist(savedir,'dir'), mkdir(savedir); end

switch out_type
    case 'day' % get DOY (fractionalday in table is not necessarily correct...)
        tmp=day(table_in.DateTime,'dayofyear');
    case 'week' % get week number (Matlab starts the week on Sunday...)
        tmp=week(table_in.DateTime);
        % redefine Sunday as end of the week 
        tmp2=day(table_in.DateTime,'dayofweek'); % returns 1 for Sunday
        tmp(tmp2==1)=tmp(tmp2==1)-1;
    case 'month' % get month number
        tmp=month(table_in.DateTime);
end
     
% loop over days/weeks/months
for i=unique(tmp)'

    % select data from given period
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
        
end