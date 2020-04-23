function write_input_HEIPRO( p_num, uvvis, yr_in, file_dir, out_type, split_scans )
% Write Pandora dSCDs to input files used by HEIPRO
%
%%   
% INPUT:    Code loads MAX-DOAS data processed by reformat_pandora_for_retrievals.m
%           Data is filtered, gaps are dealt with, and 90 deg dummies are
%           inserted by that function.
%           p_num: Pandora instrument number. Add location details as
%                  new instruments are included
%           uvvis: 'uv or 'vis', according to wavelength range used for
%                  processing (used for output file name only)
%           yr_in: year of data to process
%           file_dir: directory where processed data is found
%           out_type: 'day', 'week', or 'month' for daily, weekly or monthly files
%           split_scans: 0 to write all scans into one file
%                        1 to separate short and long scans into different files
%                        2 to do 1 and 2 at the same time
%
% OUTPUT:   daily/weekly/monthly dSCD files with dummy 90 deg lines inserted
%           Data is saved in a new folder in file_dir
%
%
% Kristof Bognar, August 2019

%% load data
% files saved by reformat_pandora_for_retrievals.m

error('match mapa input file naming')

if ~strcmp(file_dir(end),'/'), file_dir=[file_dir, '/']; end

fname=[file_dir 'p' num2str(p_num) '_' uvvis '_' num2str(yr_in) '_maxdoas_processed.mat'];

try
    load(fname);
catch
    error('Check data folder, or save data using reformat_pandora_for_retrievals.m')
end

table_in=pan_maxdoas_processed;


%% write files

% create output folder
% data is saved separately for each instument and wavelength range
savedir=[file_dir 'HEIPRO_input_p' num2str(p_num) '_' uvvis '/'];

% create folder if necessary
if ~exist(savedir,'dir'), mkdir(savedir), end

switch out_type
    case 'day' % get DOY (fractionalday in table is not necessarily correct...)
        time_steps=day(table_in.DateTime,'dayofyear');
    case 'week' % get week number (Matlab starts the week on Sunday...)
        time_steps=week(table_in.DateTime);
        % redefine Sunday as end of the week 
        tmp2=day(table_in.DateTime,'dayofweek'); % returns 1 for Sunday
        time_steps(tmp2==1)=time_steps(tmp2==1)-1;
    case 'month' % get month number
        time_steps=month(table_in.DateTime);
end
     
% loop over days/weeks/months
for i=unique(time_steps)'

    % select data from given period
    to_write=table_in(time_steps==i,:);

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
        
        % output file name
        fname=[savedir 'long_p' num2str(p_num) '_' uvvis '_' num2str(yr_in)...
               '_' out_type num2str(i) '.dat'];
        
        to_write_tmp.longscan=[];
        to_write_tmp.shortscan=[];
        writetable(to_write_tmp,fname,'Delimiter',',');
        
        
        %%% short scans
        to_write_tmp=to_write(to_write.shortscan==1,:);
        
        % output file name
        fname=[savedir 'short_p' num2str(p_num) '_' uvvis '_' num2str(yr_in)...
               '_' out_type num2str(i) '.dat'];
        
        to_write_tmp.longscan=[];
        to_write_tmp.shortscan=[];
        writetable(to_write_tmp,fname,'Delimiter',',');
        
    end
    
    if split_scans~=1
        % write all scans in one file
        
        % output file name
        fname=[savedir 'p' num2str(p_num) '_' uvvis '_' num2str(yr_in)...
               '_' out_type num2str(i) '.dat'];
        
        to_write.longscan=[];
        to_write.shortscan=[];
        writetable(to_write,fname,'Delimiter',',');

    end
        
end
        
end