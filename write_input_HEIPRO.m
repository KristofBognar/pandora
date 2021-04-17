function write_input_HEIPRO( p_num, uvvis, yr_in, file_dir, scan_type )
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
%           scan_type: 'long' for long scans only
%                      'short' for short scans only
%
% OUTPUT:   daily dSCD files with dummy 90 deg lines inserted
%           Data is saved in a new folder in file_dir
%
%
% Kristof Bognar, August 2019

% fix daily output, code to find dt and start/end times works for daily
% input only
out_type='day';

%% set input/output folder naming
if ~strcmp(file_dir(end),'/'), file_dir=[file_dir, '/']; end

% last folder in path
input_version=strsplit(file_dir,'/');
input_version=input_version{end-1};

% default location
if strcmp(input_version,'retrieval_input'), input_version=''; end

% output folder name, given input and a priori settings
output_version=['p' num2str(p_num) '_' input_version '/'];

disp(' ')
disp('Output folder based on input file location and settings in the code:')
disp(output_version)
disp('Continue? [y]/n')

tmp=input('','s');

if isempty(tmp) || strcmpi(tmp,'y')
    input_version=[input_version '/'];
    file_dir=file_dir(1:end-length(input_version));
else
    return
end
    

%% load data
% files saved by reformat_pandora_for_retrievals.m
fname=[file_dir input_version 'p' num2str(p_num) '_' uvvis '_' num2str(yr_in) ...
       '_maxdoas_processed.mat'];

try
    load(fname);
catch
    error('Check data folder, or save data using reformat_pandora_for_retrievals.m')
end

table_in=pan_maxdoas_processed;


%% write files

% create output folder
% data is saved separately for each instument and wavelength range
savedir_dscd=[file_dir 'HEIPRO_input/' output_version 'dSCD_files/'];
savedir_inp=[file_dir 'HEIPRO_input/' output_version 'input_files/'];

% create folder if necessary
if ~exist(savedir_dscd,'dir'), mkdir(savedir_dscd), end
if ~exist(savedir_inp,'dir'), mkdir(savedir_inp), end

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
steps_list=unique(time_steps);

daily_times_all=[];
for i=1:length(steps_list)

    % select data from given period
    to_write=table_in(time_steps==steps_list(i),:);

    % discard partial scans before first and after last 90deg dummy
    % there are a very large number of scans even for each day, no point
    % trying to save some of the scans that might be almost complete
    ind=find(to_write.Elevviewingangle==90);
    % use 90 immediately before/after scans (discard double 90s)
    if ind(2)-ind(1)==1, ind(1)=[]; end
    if ind(end)-ind(end-1)==1, ind(end)=[]; end

    to_write=to_write(ind(1):ind(end),:);

    % write files
    if strcmp(scan_type,'long')
        
        %%% long scans
        to_write_tmp=to_write(to_write.longscan==1,:);
        
        % get dt and start/end times
        [ daily_times, to_write_tmp ] = find_dt_HEIPRO(to_write_tmp,scan_type);
        
        % output file name
        dscd_fname=['long_p' num2str(p_num) '_' uvvis '_' num2str(yr_in)...
               '_' out_type num2str(steps_list(i)) '.dat'];
        fname=[savedir_dscd dscd_fname];
        
        % write dscd file
        to_write_tmp.longscan=[];
        to_write_tmp.shortscan=[];
        writetable(to_write_tmp,fname,'Delimiter',',');
        
        % create heipro control files
        date_in=to_write_tmp.DateTime(1);
        create_retrieval_files_pandora(output_version,uvvis,dscd_fname,date_in,daily_times,savedir_inp);
        
        % save daily start/end times and dt
        daily_times_all=[daily_times_all; [{datestr(date_in,'yyyymmdd')}, daily_times]];

    elseif strcmp(scan_type,'short')

        %%% long scans
        to_write_tmp=to_write(to_write.shortscan==1,:);
        
        % get dt and start/end times
        [ daily_times, to_write_tmp ] = find_dt_HEIPRO(to_write_tmp,scan_type);
        
        % output file name
        dscd_fname=['short_p' num2str(p_num) '_' uvvis '_' num2str(yr_in)...
               '_' out_type num2str(steps_list(i)) '.dat'];
        fname=[savedir_dscd dscd_fname];
        
        % write dscd file
        to_write_tmp.longscan=[];
        to_write_tmp.shortscan=[];
        writetable(to_write_tmp,fname,'Delimiter',',');
        
        % create heipro control files
        date_in=to_write_tmp.DateTime(1);
        create_retrieval_files_pandora([output_version(1:end-1) '_short'],...
                                       uvvis,dscd_fname,date_in,daily_times,savedir_inp);
        
        % save daily start/end times and dt
        daily_times_all=[daily_times_all; [{datestr(date_in,'yyyymmdd')}, daily_times]];
        
    end
            
end
      
save([file_dir 'HEIPRO_input/' output_version scan_type '_p' num2str(p_num) '_' uvvis '_' ...
      num2str(yr_in) '_retr_times.mat'],'daily_times_all')

end

