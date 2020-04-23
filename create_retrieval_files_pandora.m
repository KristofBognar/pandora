function create_retrieval_files_pandora(pnum,uvvis)
% Create daily retrieval settings file for pandora aerosol and NO2 retrievals
%
% dSCD files created by reformat_pandora_for_retrievals.m and write_input_HEIPRO.m
%
% Enter start-end times and intervals manually, so files are checked over before the
% retrievals; there are too many exceptions and issues for automating the process
%   
% Aerosol a priori is fixed
% BrO a priori determined from BrO dSCDs and T inversion
%       Surface conc.: 5 when 5deg dscd consistently above 1e14
%       Scale height: from get_sonde_PT.m (check T profile plots as well)
% 
% retrieval performed on Windows side; separate code (python) modifies
% yealy control file with sonde filenames, and copies appropriate retrieval
% file to xx_retrieval.inp
%
% @Kristof Bognar, sometime in 2018
%
%% setup

% aer=0; % 1 for aerosol, 0 for BrO
% year=2018;

%%%%%%%%%%%%% do it as a function of pnum and uvvis
error('to do: out folder automation, date/time selection, dscd file naming')

out_folder_name=['p' num2str(pnum) '_' uvvis];


% start/stop times and BrO a priori are selected manually, and saved
%%%%%%%%%%%%%%%%[dates,daily_times,apriori_BrO] = variable_init(year);

%% read template
for aer=[1,0]
    
    if aer
        fid=fopen(['/media/kristof/Windows7_OS/SCIATRAN2/AEROSOL_RETRIEVAL_v-1-2/',...
                       'IDL_execute/aerosol_retrieval_template_pandora_vis.inp'],'r');
    else
        fid=fopen(['/media/kristof/Windows7_OS/SCIATRAN2/TRACEGAS_RETRIEVAL_v-1-2/',...
                       'IDL_execute/tracegas_retrieval_template_pandora vis.inp'],'r');
    end

    i = 1;
    tline = fgetl(fid);
    infile{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        infile{i} = tline;
    end
    fclose(fid);

    %% loop over desired dates
    for i=1:length(dates)

        outfile=infile;

        %% modify relevant lines

        %%% look for field description for dates/times, since placeholder and field
        %%% header repeat

        % start date
        tmp=find_in_file(outfile,'Start day of considered time period');
        outfile{tmp+2}=datestr(dates(i),'dd/mm/yyyy');

        % stop date
        tmp=find_in_file(outfile,'Stop day of considered time period');
        outfile{tmp+2}=datestr(dates(i),'dd/mm/yyyy');

        % start time
        tmp=find_in_file(outfile,'Start time for the considered day');
        outfile{tmp+2}=daily_times{i,1};

        % stop time
        tmp=find_in_file(outfile,'Stop time for the considered day');
        outfile{tmp+2}=daily_times{i,2};

        % time interval
        tmp=find_in_file(outfile,'Interval length for the temporal retrieval sequences');
        outfile{tmp+2}=daily_times{i,3};

        %%% output folder
        tmp=find_in_file(outfile,'Folder for the output (will be generated if necessary)');
        outfile{tmp+2}=out_folder_name;

        %%% look for placeholder for dSCD and sonde files

        % dSCD input file
        tmp=find_in_file(outfile,'dscdxxxx.dat');
    % % %     outfile{tmp}=['DSCD_' datestr(dates(i),'yyyy_mm_dd') '.dat'];

        % radiosonde input file (indices specific to path!!)
        tmp=find_in_file(outfile,'radiosonde_xxxxxx.dat');
    %     outfile{tmp}(end-9:end-4)=datestr(dates(i),'yymmdd');

        % directory with aerosol results for tg retrieval
        if ~aer
            tmp=find_in_file(outfile,'Path to the aerosol retrieval results');
            outfile{tmp+3}=['C:\SCIATRAN2\AEROSOL_RETRIEVAL_v-1-2\Campaign\'...
                            out_folder_name '\'];
        end


        %% write file

        fpath=['/home/kristof/Drive/data_for_profile_retrieval/input_files/pandora/' ...
               out_folder_name '/'];

        if ~exist(fpath,'dir'), mkdir(fpath); end

        if aer
            fname=['aerosol_retrieval_' out_folder_name '_' datestr(dates(i),'yymmdd') '.inp'];
        else
            fname=['tracegas_retrieval_' out_folder_name '_' datestr(dates(i),'yymmdd') '.inp'];
        end

        file_tmp = fopen([fpath fname],'w');
        fprintf(file_tmp,'%s\n',outfile{:});
        fclose(file_tmp);        


    end
    
end

end



function [dates,daily_times,apriori_BrO] = variable_init(year)

    %%% start and end times + time interval (match manually to times in dSCD file)
    % (start,stop,interval, mmdd, comments), and
    %%% BrO a priori
    % determined by visual insection of dSCDs and T inversion (Zhao et al., 2016)
    % surf VMR (x1e-6), scake height (km), mmdd
        
    if year==2019
        
        % day range 
        doy_range=[64:151]; 

        % remove missing days
        missing=[112,118];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);
        
        % times and apriori info        
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2019.mat')
       
    end
    
    if year==2018
        %% dates

        % day range 
        doy_range=[64:151]; 
%         doy_range=[64:151,213:216,218:221,233,234]; 

        % remove missing days
        doy_range(doy_range==133)=[]; 
        
        % dates
        dates=ft_to_date(doy_range-1,year);

        % times and apriori info        
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2018.mat')
                 
    elseif year==2017
        %% dates

        % day range 
        doy_range=[66:94]; 

        % remove missing days
%         doy_range(doy_range==133)=[]; 

        % dates
        dates=ft_to_date(doy_range-1,year);
        
        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2017.mat')       
   
    elseif year==2016
        %% dates

        % day range 
        doy_range=[66:132]; 

        % remove missing days
        missing=[73:75,84:91];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);

        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2016.mat')       

    elseif year==2015
        %% dates

        % day range 
        doy_range=[62:151]; 

        % remove missing days
        missing=[63,87,98,99,106,123,145];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);

        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2015.mat')       
      
    elseif year==2013
        %% dates

        % day range 
        doy_range=[69:113]; 

        % remove missing days
        doy_range(doy_range==93)=[]; 

        % dates
        dates=ft_to_date(doy_range-1,year);
                
        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2013.mat')       

    elseif year==2011
        %% dates

        % day range 
        doy_range=[70:100]; 

        % remove missing days
        missing=[83,84,96,98];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);
                
        % times and apriori info
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2011.mat')       
        
    elseif year==2010
        %% dates

        % day range 
        %%% there are more measurements; only process part of the dataset
        doy_range=[100:140]; 

        % remove days that are not needed
        missing=[107:122,130];
        doy_range=setdiff(doy_range, missing);

        % dates
        dates=ft_to_date(doy_range-1,year);
                
        % times and apriori info
        % used 1h profile time window for most days, since scans are very
        % long for some reason (there are lots of gaps too)
        load('/home/kristof/work/profile_retrievals/profile_results/profile_details/prof_info_2010.mat')       
        

    end
end




