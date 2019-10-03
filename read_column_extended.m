function [ data ] = read_column_extended( f_name, ncol )
% %[ data ] = read_column_extended( f_name )
% 
% Read in spectra from QDOAS column extended format
% If file contains wavelength calibration and spectrum counts, set ncol=2
%
% Output is a table with fields of
% {'DateTime','view_el','view_az','t_meas','t_exp','SZA','SAA','spec'}
%
% Code is grossly inefficient since it reads file line by line
%
% Should be able to read files with different headers as long as the key
% lines (info that goes in the table) are unchanged
%
%@Kristof Bognar, July 2019


%f_name=['/home/kristof/work/PANDORA/profiling_test/',...
%        'OFFAXIS_OPEN_Pandora103s1_Downsview_20180214_L1_smca1c2p1-5.spe'];
%/home/kristof/atmosp_servers/net/aurora/ground/toronto/pandora/Pandora104/spe_BlickP_L1_v2

%% setup the variables

if nargin==1, ncol=1; end

% detector size -- number of lines per spectrum
det_size=NaN;

% number of spectra
num_spec=NaN;

% number of header lines per spectrum
num_head=0;


% varnames={'Date','Time','View_el','View_az','Type','t_meas','t_exp','SZA','SAA'};
% data=cell2table(cell(0,length(varnames)), 'VariableNames', varnames);

DateTime_arr=[];
view_el_arr=[];
view_az_arr=[];
% type_arr=[];
t_meas_arr=[];
t_exp_arr=[];
SZA_arr=[];
SAA_arr=[];

spec=[];
wl=[];
tmp_wl=[];
skip_line=0;

% temporary variables
date_tmp=[];
time_tmp=[];
view_el_tmp=[];
view_az_tmp=[];
t_meas_tmp=[];
t_exp_tmp=[];
SZA_tmp=[];
SAA_tmp=[];


%% read file

fid=fopen(f_name,'r');

% loop over file line by line 
find_header=1;
    

%% read file header and first spectrum header
while find_header

    % get next line
    line = fgets(fid);

    % find detector size
    if find_in_cell(cellstr(line),'# Size of the detector',true),
        % break up by equal sign
        ind=cell2mat(strfind(cellstr(line),'='));
        % convert stuff following equal sign to a number (text might follow)
        det_size=str2double(line(ind+1:end));
        if isempty(det_size) || isnan(det_size)
            det_size=str2double(line(ind+1:ind+6));
        end
    end

    % find number of spectra
    if find_in_cell(cellstr(line),'# Total number of records',true),
        % break up by equal sign
        ind=cell2mat(strfind(cellstr(line),'='));
        % convert stuff following equal sign to a number (text might follow)
        num_spec=str2double(line(ind+1:end));
        if isempty(num_spec) || isnan(num_spec)
            num_spec=str2double(line(ind+1:ind+6));
        end
    end

    % figure out how many lines in spectrum headers (and read first one)
    if isletter(line(1))

        num_head=num_head+1;

        [date_tmp,time_tmp,view_el_tmp,view_az_tmp,t_meas_tmp,t_exp_tmp,...
                  SZA_tmp,SAA_tmp]=assign_header(line,...
         date_tmp,time_tmp,view_el_tmp,view_az_tmp,t_meas_tmp,t_exp_tmp,SZA_tmp,SAA_tmp);
              
    elseif ~strcmp(line(1),'#')

        %turn off header bloc and save datetime
        find_header=0;
        DateTime_arr=datetime([date_tmp,time_tmp],'inputformat','dd/MM/yyyy HH:mm:ss');
          
        view_el_arr=view_el_tmp;
        view_az_arr=view_az_tmp;
        t_meas_arr=t_meas_tmp;
        t_exp_arr=t_exp_tmp;
        SZA_arr=SZA_tmp;
        SAA_arr=SAA_tmp;
%         type_arr=type_tpm;
  
        % read first spectrum
        if ncol==2
            % split by tab, first number is wavelength, second is counts
            % first entry might be a tab
            tmp=strsplit(line,'\t');
            tmp_wl=str2double(tmp{end-1});
            tmp_sp=str2double(tmp{end});
        else
            tmp_sp=str2double(line);
        end
    
        % test if there are empty lines between spectrum counts
        test_line=fgets(fid);
        
        if length(test_line)<=2 %% might not be 2 for all files...
            
            skip_line=1;

            for i=2:det_size
                
                line=fgets(fid);
                tmp=fgets(fid);
                
                if ncol==2
                    % split by tab, first number is wavelength, second is counts
                    % first entry might be a tab
                    tmp=strsplit(line,'\t');
                    tmp_wl(i)=str2double(tmp{end-1});
                    tmp_sp(i)=str2double(tmp{end});
                else
                    tmp_sp(i)=str2double(line);
                end
                
            end
            
            
        else
            
            tmp=strsplit(test_line,'\t');
            tmp_wl(2)=str2double(tmp{end-1});
            tmp_sp(2)=str2double(tmp{end});
            
            for i=3:det_size
                
                line = fgets(fid);
                
                if ncol==2
                    % split by tab, first number is wavelength, second is counts
                    % first entry might be a tab
                    tmp=strsplit(line,'\t');
                    tmp_wl(i)=str2double(tmp{end-1});
                    tmp_sp(i)=str2double(tmp{end});
                else
                    tmp_sp(i)=str2double(line);
                end
            end

            
        end
        
        spec=tmp_sp;
        wl=tmp_wl;
        
    end

end


%% read spectra and headers

for i=2:num_spec

    % spectrum header
    for j=1:num_head

        % get next line
        line = fgets(fid);

        [date_tmp,time_tmp,view_el_tmp,view_az_tmp,t_meas_tmp,t_exp_tmp,...
                  SZA_tmp,SAA_tmp]=assign_header(line,...
         date_tmp,time_tmp,view_el_tmp,view_az_tmp,t_meas_tmp,t_exp_tmp,SZA_tmp,SAA_tmp);

    end

    DateTime_arr=[DateTime_arr;...
          datetime([date_tmp,time_tmp],'inputformat','dd/MM/yyyy HH:mm:ss')];

    view_el_arr=[view_el_arr;view_el_tmp];
    view_az_arr=[view_az_arr;view_az_tmp];
    t_meas_arr=[t_meas_arr;t_meas_tmp];
    t_exp_arr=[t_exp_arr;t_exp_tmp];
    SZA_arr=[SZA_arr;SZA_tmp];
    SAA_arr=[SAA_arr;SAA_tmp];
%         type_arr=[type_arr,type_tmp];
    
    
    % spectrum 
    for j=1:det_size

        line=fgets(fid);
        if skip_line, tmp=fgets(fid); end
        
        if ncol==2
            % split by tab, first number is wavelength, second is counts
            % first entry might be a tab
            tmp=strsplit(line,'\t');
            tmp_wl(j)=str2double(tmp{end-1});
            tmp_sp(j)=str2double(tmp{end});
        else
            tmp_sp(j)=str2double(line);
        end

    end

    spec=[spec;tmp_sp];
    wl=[wl;tmp_wl];
    
end

% test if file really ends
line = fgets(fid);
if ischar(line), error('File longer than expected'), end

    
fclose(fid);


%% addign results

data=table;

data.DateTime=DateTime_arr;
data.view_el=view_el_arr;
data.view_az=view_az_arr;
data.t_meas=t_meas_arr;
data.t_exp=t_exp_arr;
data.SZA=SZA_arr;
data.SAA=SAA_arr;

data.spec=spec;

end

function [date_tmp,time_tmp,view_el_tmp,view_az_tmp,t_meas_tmp,t_exp_tmp,...
          SZA_tmp,SAA_tmp]=assign_header(line,...
          date_tmp,time_tmp,view_el_tmp,view_az_tmp,t_meas_tmp,t_exp_tmp,SZA_tmp,SAA_tmp)


        % read date
        if find_in_cell(cellstr(line),'Date',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            date_tmp=line(ind+1:end);
        end

        % read time
        if find_in_cell(cellstr(line),'UTC Time',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            time_tmp=line(ind+1:end);
        end

        % read viewing elev
        if find_in_cell(cellstr(line),'Viewing Elevation',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            view_el_tmp=str2double(line(ind+1:end));
        end
        
        % read viewing az
        if find_in_cell(cellstr(line),'Viewing Azimuth',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            view_az_tmp=str2double(line(ind+1:end));
        end
        
        % read total meas time
        if find_in_cell(cellstr(line),'Total Measurement Time',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            t_meas_tmp=str2double(line(ind+1:end));
        end
        
        % read exposure time
        if find_in_cell(cellstr(line),'Exposure Time',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            t_exp_tmp=str2double(line(ind+1:end));
        end
        
        % read SZA
        if find_in_cell(cellstr(line),'Solar Zenith Angle',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            SZA_tmp=str2double(line(ind+1:end));
        end
        
        % read SAA
        if find_in_cell(cellstr(line),'Solar Azimuth Angle',true)
            ind=cell2mat(strfind(cellstr(line),'='));
            ind2=cell2mat(strfind(cellstr(line),'('));
            SAA_tmp=str2double(line(ind(1)+1:ind2(2)-1));
        end
%             % read meas type
%             if find_in_cell(cellstr(line),'Measurement Type',true)
%                 ind=cell2mat(strfind(cellstr(line),'='));
%                 type_tmp=line(ind+1:end);
%             end

end








% 
% end
% 
