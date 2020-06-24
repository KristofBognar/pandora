function create_retrieval_files_pandora(out_folder_name,dscd_fname,date_in,daily_times,fpath)
% Create daily retrieval settings file for pandora aerosol and NO2 retrievals
%
% % % % dSCD files created by reformat_pandora_for_retrievals.m and write_input_HEIPRO.m
% % % %
% % % % Enter start-end times and intervals manually, so files are checked over before the
% % % % retrievals; there are too many exceptions and issues for automating the process
% % % %   
% % % % Aerosol a priori is fixed
% % % % BrO a priori determined from BrO dSCDs and T inversion
% % % %       Surface conc.: 5 when 5deg dscd consistently above 1e14
% % % %       Scale height: from get_sonde_PT.m (check T profile plots as well)
% % % % 
% % % % retrieval performed on Windows side; separate code (python) modifies
% % % % yealy control file with sonde filenames, and copies appropriate retrieval
% % % % file to xx_retrieval.inp
%
% @Kristof Bognar, June 2020
%

if out_folder_name(end)=='/', out_folder_name=out_folder_name(1:end-1); end

%% instrument altitude (km above surface, not sea level)

p_num=str2double(out_folder_name(2:4));

if any(p_num==[103,104])
 
    % ECCC Downsview: 15m
    instr_alt_str='0.015';
    loc_str='Downsview';
        
else
    error(['Add location details for Pandora ' num2str(instr_num)]);
end

%% do both aerosol and tg input
for aer=[1,0]
    %% read template
    if aer
        fid=fopen(['/media/kristof/Windows7_OS/SCIATRAN2/AEROSOL_RETRIEVAL_v-1-2/',...
                       'IDL_execute/aerosol_retrieval_template_pandora_vis.inp'],'r');
    else
        fid=fopen(['/media/kristof/Windows7_OS/SCIATRAN2/TRACEGAS_RETRIEVAL_v-1-2/',...
                       'IDL_execute/tracegas_retrieval_template_pandora_vis.inp'],'r');
    end

    i = 1;
    infile={};
    tline = fgetl(fid);
    infile{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        infile{i} = tline;
    end
    fclose(fid);

    outfile=infile;

    %% modify relevant lines

    %%% look for field description for dates/times, since placeholder and field
    %%% header repeat

    % start date
    tmp=find_in_file(outfile,'Start day of considered time period');
    outfile{tmp+2}=datestr(date_in,'dd/mm/yyyy');

    % stop date
    tmp=find_in_file(outfile,'Stop day of considered time period');
    outfile{tmp+2}=datestr(date_in,'dd/mm/yyyy');

    % start time
    tmp=find_in_file(outfile,'Start time for the considered day');
    outfile{tmp+2}=daily_times{1};

    % stop time
    tmp=find_in_file(outfile,'Stop time for the considered day');
    outfile{tmp+2}=daily_times{2};

    % time interval
    tmp=find_in_file(outfile,'Interval length for the temporal retrieval sequences');
    outfile{tmp+2}=daily_times{3};

    %%% output folder
    tmp=find_in_file(outfile,'Folder for the output (will be generated if necessary)');
    outfile{tmp+2}=out_folder_name;

    %%% instrument position
    tmp=find_in_file(outfile,'# Height position of instrument in [km]');
    outfile{tmp+2}=instr_alt_str;
    
    %%% look for placeholder for dSCD and sonde files

    % dSCD input file
    tmp=find_in_file(outfile,'dscd_xxxx.dat');
    outfile{tmp}=dscd_fname;

    % radiosonde input file (indices specific to path!!)
    tmp=find_in_file(outfile,'ERA5_xxxxxx.dat');
%     outfile{tmp}(end-9:end-4)=datestr(date_in,'yymmdd');
    outfile{tmp}=['C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\Complement_Data\Radiosonde\' ...
                  'ERA5_' loc_str '_' datestr(date_in,'yymmdd') '.dat'];

    % directory with aerosol results for tg retrieval
    if ~aer
        tmp=find_in_file(outfile,'Path to the aerosol retrieval results');
        outfile{tmp+3}=['C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\'...
                        out_folder_name '\'];
    end


    %% write file

%         fpath=['/home/kristof/Drive/data_for_profile_retrieval/input_files/pandora/' ...
%                out_folder_name '/'];

    if ~exist(fpath,'dir'), mkdir(fpath); end

    if aer
        fname=['aerosol_retrieval_' out_folder_name '_' datestr(date_in,'yymmdd') '.inp'];
    else
        fname=['tracegas_retrieval_' out_folder_name '_' datestr(date_in,'yymmdd') '.inp'];
    end
    
    file_tmp = fopen([fpath fname],'w');
    fprintf(file_tmp,'%s\n',outfile{:});
    fclose(file_tmp);        

    
end

end

