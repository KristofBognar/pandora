function write_input_MAPA( table_in, savedir, out_type, split_scans )
% Write Pandora dSCDs to input files used by MAPA
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
% OUTPUT:   daily/weekly/monthly dSCD files in MAPA's netCDF format
%
%
% Kristof Bognar, January 2020

% sample input file
% '/home/kristof/work/MAPA/Release_0.991/sample_files/input/sample/sample_input_file_v0991.nc'

%% reorganize input data

%%% debug
% load('pan_103_vis_2018_maxdoas_processed.mat')
% table_in=pan_maxdoas_processed;
%%%

% calculate relative azimuth angle
table_in.RAA=table_in.SolarAzimuthAngle-(table_in.Azimviewingangle-180);

% keep relevant columns only
columns_orig = table_in.Properties.VariableNames;
columns_new = {'DateTime','Elevviewingangle','Azimviewingangle',...
               'SZA','SolarAzimuthAngle','RAA',...     
               'NO2_VisSlColno2','NO2_VisSlErrno2','NO2_VisSlColo4','NO2_VisSlErro4'};
[~,tmp] = ismember(columns_new,columns_orig);
table_in = table_in(:,tmp);


%% elevation angles -- some of the code is hard coded for:
%           30,     15,              2, 1, and up (tagged as shortscan==1)
%   50, 40, 30, 20, 15, 10, 8, 5, 3, 2, 1, and up (tagged as longscan==1)

% all possible angles (must be sorted from low to high!)
elevs_all=[1,2,3,5,8,10,15,20,30,40,50,90];
elevs_saved=[1,2,3,5,8,10,15,20,30,40,50];

% assuming down-up scans: what is the lowest elevation
middle_elev=1; % 1deg measurement

% full down and up scan
elevs_long=[fliplr(elevs_all),elevs_all(2:end)]; 
elevs_short=[90,30,15,2,1,2,15,30,90];

% indices for rearranging short and long arrays to scan# by elevs format
lines_empty=array2table(NaN(length(elevs_all)*2,length(columns_new)),...
                        'Variablenames',columns_new);
lines_empty.DateTime=-1*ones(size(lines_empty.DateTime));
                    
lines_empty_len=size(lines_empty,1);

%% break up data by date

switch out_type
    case 'day' % get DOY (fractionalday in table is not necessarily correct...)
        time_steps=day(table_in.DateTime,'dayofyear');
    case 'week' % get week number (Matlab starts the week on Sunday...)
        time_steps=week(table_in.DateTime);
        % put sundays as end of the week 
        tmp2=day(table_in.DateTime,'dayofweek'); % returns 1 for sunday
        time_steps(tmp2==1)=time_steps(tmp2==1)-1;
    case 'month' % get month number
        time_steps=month(table_in.DateTime);
end

% reformat datetime (need integers, but keep it as double for processing)
table_in.DateTime=str2num(datestr(table_in.DateTime,'yyyymmddHHMMSS'));


%% loop over days/weeks/months
for i=unique(time_steps)'

    % select data from given period
    to_write=table_in(time_steps==i,:);
    
    % discard double 90deg lines
    ind_90=find(to_write.Elevviewingangle==90);
    diffs=[0;diff(ind_90)];

    to_write(ind_90(diffs==1),:)=[];
    
    % 90deg lines now indicate finished down-up scans
    ind_90=find(to_write.Elevviewingangle==90);
    
    % loop over each down-up scan
    for j=1:length(ind_90)
        
        % first scan might be truncated
        if j==1 
            if ind_90(j)~=1
                % data before 1st 90deg: partial scan from prev. UTC day -- remove
                to_write(1:ind_90(j)-1,:)=[];
                % recalc 90deg indices
                ind_90=find(to_write.Elevviewingangle==90);
            end
            
            continue
        end
                    
        % double check 90deg indices
        if to_write.Elevviewingangle(ind_90(j))~=90 || ...
           to_write.Elevviewingangle(ind_90(j-1))~=90 
            disp([int64(to_write.DateTime(ind_90(j-1))),int64(to_write.DateTime(ind_90(j)))])
            error('90 deg indices are off')
        end
        
        % indices of current down-up scan
        seq_curr=ind_90(j-1):ind_90(j);
        
        % current elevs
        elevs_curr=to_write.Elevviewingangle(seq_curr);
        
        %% check for messed up scans first
        
        % number of repetitions for each angle
        tmp=hist(elevs_curr,unique(elevs_curr));
        
        if any(tmp>2) ||... % too many measurements, likely a jumble of truncated scans
           length(seq_curr)<6 % only part of a short scan, not enough for retrievals
            
            % remove current segment, except for last 90deg line (seq_curr(end))
            to_write(seq_curr(1:end-1),:)=[];
            
            % account for removed lines
            ind_90=ind_90-(length(seq_curr)-1);
            
            % move on to next chunk
            continue
        end

        %% reorganize good scans
        
        if isequal(elevs_curr',elevs_short) 
            % if it's a complete short scan, assign accordingly

            % reorder current down-up sequence as up-up, with 1deg
            % measurement duplicated
            seq_ordered=[ind_90(j-1)+4:-1:ind_90(j-1),ind_90(j-1)+4:ind_90(j)];

            % copy empty table segment that fits 2 long scans
            tmp=lines_empty;

            % fill in entries from short scan, leaving the rest as NaN
            tmp([1,2,7,9,12,13,14,19,21,24],:)=to_write(seq_ordered,:);

            % replace original down-up short scan with new expanded segment
            to_write=[to_write(1:seq_curr(1)-1,:);...
                      tmp;...
                      to_write(seq_curr(end)+1:end,:)];

        elseif isequal(elevs_curr',elevs_long)
            % if it's a complete short scan, assign accordingly

            % reorder current down-up sequence as up-up, with 1deg
            % measurement duplicated
            seq_ordered=[ind_90(j-1)+11:-1:ind_90(j-1),ind_90(j-1)+11:ind_90(j)];
            
            to_write=[to_write(1:seq_curr(1)-1,:);...
                      to_write(seq_ordered,:);...
                      to_write(seq_curr(end)+1:end,:)];

        elseif length(elevs_curr)<length(elevs_long)
            % some elev angles missing
            
            % copy empty table segment that fits 2 long scans
            lines_to_fill=lines_empty;
            
            % fill with available angles:
            % loop through ONE set of angles!
            middle_ind=[];
            for k=1:length(elevs_all)
                
                % occurences of selected angle in the data
                elev_ind=find(elevs_curr==elevs_all(k));
                
                if elevs_all(k)==middle_elev
                    % deal with 1deg measurement
                    
                    if length(elev_ind)==1
                        % if elev angle is 1, duplicate
                        lines_to_fill(k,:)=to_write(seq_curr(elev_ind),:);
                        lines_to_fill(k+length(elevs_all),:)=to_write(seq_curr(elev_ind),:);
                        % save as middle of the scan for later
                        middle_ind=elev_ind; 
                    else
                        % no 1deg measurement -- middle of the scan is
                        % where elevation angles start to increase
                        to_sort=[diff(elevs_curr)',0];
                        [~,tmp]=sort(to_sort);
                        middle_ind=find(tmp==length(to_sort))-0.5; % -0.5 to avoid equalities
                    end
                    
                else
                    % all other elevs

                    if length(elev_ind)==2
                        % if angle appears twice, assign both to the desired
                        % location in the extended array
                        lines_to_fill(k,:)=to_write(seq_curr(elev_ind(1)),:);
                        lines_to_fill(k+length(elevs_all),:)=to_write(seq_curr(elev_ind(2)),:);
                    elseif length(elev_ind)==1
                        % if angle appears once:
                        % need to find if first or second scan
                        if elev_ind<middle_ind
                            lines_to_fill(k,:)=to_write(seq_curr(elev_ind),:);
                        else
                            lines_to_fill(k+length(elevs_all),:)=to_write(seq_curr(elev_ind),:);
                        end
                    end
                
                end
                    
                    
            end
            
            % replace original down-up short scan with new expanded segment
            to_write=[to_write(1:seq_curr(1)-1,:);...
                      lines_to_fill;...
                      to_write(seq_curr(end)+1:end,:)];

        else
            error('No scans should be here...')
        end
              
        % account for the number of extra rows added
        ind_90=ind_90+(lines_empty_len-length(seq_curr));
        
        % last line is not 90 (scan crosses UTC midnight)
        if j==length(ind_90) && ind_90(j)~=size(to_write,1)
            
            % double check 90 indices are correct
            tmp=find(to_write.Elevviewingangle==90);
            if any(tmp>ind_90(j)), error('fix this'), end
            
            % remove measurements after last 90
            to_write(ind_90(j)+1:end,:)=[];
        end
            
    end
    
    %% reformat to shape required by MAPA
    
    % remove 90deg dummies (code above removes half of them)
    to_write(to_write.Elevviewingangle==90,:)=[];
    
    % number of elevation angles
    n_elevs=length(elevs_saved); % remove 90deg from elev list
    
    % number of scans
    n_seq=size(to_write,1)/n_elevs; % dataset is now filled up so each scan is the same length 

    %%% convert to array
    to_write=table2array(to_write);
    
    % reshape as elevs by scans array, with each page coresponding to one parameter
    to_write=reshape(to_write,n_elevs,n_seq,length(columns_new));
    
    % transpose each page in 3D array to get required scans by elevs array
    to_write=permute(to_write,[2 1 3]);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Write data to file
    
    n_pt=339;

    f_out=['test_w' num2str(i) '.nc'];
    if exist(f_out,'file'), delete(f_out); end

    %% create file
    ncid = netcdf.create(f_out,'netcdf4');

    % gobal attributes
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,varid,'data_name','Pandora');
    netcdf.putAtt(ncid,varid,'data_version','1.0');
    netcdf.putAtt(ncid,varid,'data_contact','xiaoyi.zhao@canada.ca, kbognar@physics.utoronto.ca');
    netcdf.putAtt(ncid,varid,'data_institute','Environment and Climate Change Canada');

    % set dimensions
    dim_seq = netcdf.defDim(ncid,'dim_sequences',n_seq);
    dim_elevs = netcdf.defDim(ncid,'dim_elevations',n_elevs);
    dim_pt = netcdf.defDim(ncid,'dim_profile',n_pt);
    dim_scalar = netcdf.defDim(ncid,'dim_scalar',1);

    %% data groups and variables
    % aerosols
    tmp = netcdf.defGrp(ncid,'aerosol');
    grp_aer=netcdf.defGrp(tmp,'o4_425to490nm');

    var_o4 = netcdf.defVar(grp_aer,'value','NC_DOUBLE',[dim_seq,dim_elevs]);
    set_attr(grp_aer,var_o4,NaN,'o4_slant_column','molec**2/cm**5')

    var_o4err = netcdf.defVar(grp_aer,'error','NC_DOUBLE',[dim_seq,dim_elevs]);
    set_attr(grp_aer,var_o4err,NaN,'o4_slant_column_error','molec**2/cm**5')


    % tracegas
    tmp = netcdf.defGrp(ncid,'tracegas');
    grp_tg=netcdf.defGrp(tmp,'no2_425to490nm');

    var_no2 = netcdf.defVar(grp_tg,'value','NC_DOUBLE',[dim_seq,dim_elevs]);
    set_attr(grp_tg,var_no2,NaN,'no2_slant_column','molec/cm**2')

    var_no2err = netcdf.defVar(grp_tg,'error','NC_DOUBLE',[dim_seq,dim_elevs]);
    set_attr(grp_tg,var_no2err,NaN,'no2_slant_column_error','molec/cm**2')


    % measurement info
    grp_meas = netcdf.defGrp(ncid,'measurement_info');

    var_ind = netcdf.defVar(grp_meas,'scan_idx','NC_INT',dim_seq);
    set_attr(grp_meas,var_ind,-1,'index','')

    var_raa = netcdf.defVar(grp_meas,'raa','NC_DOUBLE',[dim_seq,dim_elevs]);
    set_attr(grp_meas,var_raa,NaN,'relative_azimuth_angle','degree')

    var_time = netcdf.defVar(grp_meas,'time','NC_INT64',[dim_seq,dim_elevs]);
    set_attr(grp_meas,var_time,-1,'measurement_time_(UTC)','yyyymmddHHMMSSZ')

    var_sza = netcdf.defVar(grp_meas,'sza','NC_DOUBLE',[dim_seq,dim_elevs]);
    set_attr(grp_meas,var_sza,NaN,'solar_zenith_angle','degree')

    var_elev = netcdf.defVar(grp_meas,'elevation_angle','NC_DOUBLE',dim_elevs);
    set_attr(grp_meas,var_elev,NaN,'elevation_angle','degree')


    % a priori temperature, pressure profiles
    grp_pt = netcdf.defGrp(ncid,'atmosphere');

    var_h = netcdf.defVar(grp_pt,'height','NC_DOUBLE',[dim_seq,dim_pt]);
    set_attr(grp_pt,var_h,NaN,'Height','km')

    var_t = netcdf.defVar(grp_pt,'temperature','NC_DOUBLE',[dim_seq,dim_pt]);
    set_attr(grp_pt,var_t,NaN,'Temperature','K')

    var_p = netcdf.defVar(grp_pt,'pressure','NC_DOUBLE',[dim_seq,dim_pt]);
    set_attr(grp_pt,var_p,NaN,'Pressure','hPa')


    % instrument location
    grp_loc = netcdf.defGrp(ncid,'instrument_location');

    station_alt = netcdf.defVar(grp_loc,'altitude_of_station','NC_DOUBLE',dim_scalar);
    set_attr(grp_loc,station_alt,NaN,'Altitude of the station above sea level','m')

    instr_lon = netcdf.defVar(grp_loc,'longitude','NC_DOUBLE',dim_scalar);
    set_attr(grp_loc,instr_lon,NaN,'Longitude of the instrument (positive East)','degree_east')

    instr_lat = netcdf.defVar(grp_loc,'latitude','NC_DOUBLE',dim_scalar);
    set_attr(grp_loc,instr_lat,NaN,'Latitude of the instrument (positive North)','degree_north')

    instr_alt = netcdf.defVar(grp_loc,'altitude_of_instrument','NC_DOUBLE',dim_scalar);
    set_attr(grp_loc,instr_alt,NaN,'Altitude of the instrument above sea level','m')

    %% write variables
    netcdf.endDef(ncid); %enter data mode

    %%% aerosols
    tmp=find_in_cell(columns_new,'NO2_VisSlColo4');
    netcdf.putVar(grp_aer,var_o4,to_write(:,:,tmp));

    tmp=find_in_cell(columns_new,'NO2_VisSlErro4');
    netcdf.putVar(grp_aer,var_o4err,to_write(:,:,tmp));
    
    %%% tracegas
    tmp=find_in_cell(columns_new,'NO2_VisSlColno2');
    netcdf.putVar(grp_tg,var_no2,to_write(:,:,tmp));
    
    tmp=find_in_cell(columns_new,'NO2_VisSlErrno2');
    netcdf.putVar(grp_tg,var_no2err,to_write(:,:,tmp));
    
    %%% measurement info
    % scan number
    netcdf.putVar(grp_meas,var_ind,[1:n_seq]');

    % relative azimuth angle
    tmp=find_in_cell(columns_new,'RAA');
    netcdf.putVar(grp_meas,var_raa,to_write(:,:,tmp));

    % time (format??)
    tmp=find_in_cell(columns_new,'DateTime');
    netcdf.putVar(grp_meas,var_time,int64(to_write(:,:,tmp)));

    % solar zenigh angle
    tmp=find_in_cell(columns_new,'SZA');
    netcdf.putVar(grp_meas,var_sza,to_write(:,:,tmp));

    % elevation angle
    netcdf.putVar(grp_meas,var_elev,elevs_saved);

    %%% instrument location
    % station altitude
    netcdf.putVar(grp_loc,station_alt,NaN);
    % instrument altitude
    netcdf.putVar(grp_loc,instr_alt,NaN);
    % station latitude
    netcdf.putVar(grp_loc,instr_lat,NaN);
    % station longitude
    netcdf.putVar(grp_loc,instr_lon,NaN);
    



    
    netcdf.close(ncid);
end

end

function set_attr(ncid,varid,fillval,longname,units)

    netcdf.defVarFill(ncid,varid,false,fillval);
    netcdf.putAtt(ncid,varid,'long_name',longname);
    netcdf.putAtt(ncid,varid,'units',units);

end

