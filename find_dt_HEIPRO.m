function [ times_out, table_out ] = find_dt_HEIPRO(table_in,scan_type)
%[ dt_out, t_lim, table_out ] = find_dt_HEIPRO(table_in)
%   find optimal retrieval window for daily dSCD files

%% setup

% start/end times of daily interval (minutes rounded down/up)
ft_start=table_in.DateTime(1).Hour+(table_in.DateTime(1).Minute)/60;
ft_end=table_in.DateTime(end).Hour+(table_in.DateTime(end).Minute+1)/60;

% times of 1deg and 90deg measurements
loc1=table_in.Fractionaltime(table_in.Elevviewingangle==1);
loc90=table_in.Fractionaltime(table_in.Elevviewingangle==90);

% use existence of 1deg meas as marker of complete down-up scan 
% remove scans that don't have 1deg measurement
if length(loc90)~=length(loc1)*2

    for i=2:length(loc90)
        cond=table_in.Elevviewingangle(table_in.Fractionaltime<loc90(i) & ...
                                     table_in.Fractionaltime>loc90(i-1));
        if ~isempty(cond) && ~any(cond == 1)
            table_in(table_in.Fractionaltime<=loc90(i) & ...
                     table_in.Fractionaltime>=loc90(i-1), :)=[];
        end
    end
    
end

% return filtered table
table_out=table_in;

%% get time window stats

% list of time intervals (in minutes)
if strcmp(scan_type,'long')
    dt_list=20:1:120;
elseif strcmp(scan_type,'short')
    dt_list=3:1:60;
end

% get time window stats
[out_arr,t_lim_arr]=scan_timing(dt_list,loc1,loc90,ft_start,ft_end,scan_type);

%% select optimal time interval

% number of scans
n_scans=length(loc1);

% do we have a dt that aligns scans perfectly? (need first match only)
ind=find(out_arr(:,1)==n_scans,1);

if ind % yes -- found the answer

    bad_dt=0;
    
    ind_final=ind;
    
else % no -- need to find best compromise

    bad_dt=1;
    
    % find highest number of complete down-up scans
    tmp=max(out_arr(:,1)); % gives first match only
    ind=find(out_arr(:,1)==tmp); % find all dt that give similar max. n.o. good windows

    if length(ind)==1 % clear second place
        
        ind_final=ind;

    else % multiple dt with same n.o. complete scans
        
        % decide by number of usable half-scans
        tmp=max(out_arr(ind,2)); % gives first match only
        ind2=find(out_arr(ind,2)==tmp,1); % take first match -- no practical difference if n.o. 
                                      % complete and partial scans is the same
        
        ind_final=ind(ind2);
        
    end
    
end

%% assign final values

times_out={};

% start and end times
tmp=t_lim_arr{ind_final};

times_out{1}=datestr(datetime(tmp(1)/24,'convertfrom','datenum'),'HH:MM:SS');
times_out{2}=datestr(datetime(tmp(end)/24,'convertfrom','datenum'),'HH:MM:SS');

% time interval
times_out{3}=num2str(dt_list(ind_final),'%.1f');


%% check results
if bad_dt

    disp([datestr(table_in.DateTime(1),'mmm dd') ': half (extra) scans: ' ...
         num2str(out_arr(ind_final,2)) '(' num2str(out_arr(ind_final,3)) ')/' num2str(n_scans)])
    
%     figure(99)
% 
%     plot(table_in.Fractionaltime,table_in.Elevviewingangle,'kx-'), hold on
%     for x=tmp
%         plot([x,x],[0,92],'r--','linewidth',1.5)
%     end
%     
%     title(datestr(table_in.DateTime(1),'mmm dd'))
%     ylim([-20,110])
% %     uiwait;
%     waitforbuttonpress;
%     clf;

else
    
    disp([datestr(table_in.DateTime(1), 'mmm dd') ': good dt'])
    
end

end

function [ out_arr, t_lim_arr ] = scan_timing(dt,loc1,loc90,ft_start,ft_end,scan_type)
    %%% based on input dt and start/end times, outputs number of time
    %%% windows that contain one full down-up scan
    
    % 23:59:59, to limit end time to same day
    ft_end_max=23+59/60+59/3600;

    % initialize output array
    out_arr=zeros(length(dt),3);
    t_lim_arr={};
    
    % loop over each dt value
    for i=1:length(dt)

        % time bins to consider
        t_lim=ft_start:dt(i)/60:ft_end;
        
        % append end time: add another step to make sure last scan is
        % included (make sure we stay on the same day though)
        if t_lim(end)<ft_end, t_lim=[t_lim, min(t_lim(end)+dt(i)/60,ft_end_max)]; end
        
        % if dt larger than time interval, stop
        if length(t_lim)==1, break, end
        
        % save time limits array
        t_lim_arr{i}=t_lim;
        
        if strcmp(scan_type,'long')
            
            % simple approach -- already used in retrievals
            % method for 'short' might be slightly better for long scans as
            % well, but likely not much of a difference
            
            % number of 1deg and 90deg lines in each time bin
            count1=histcounts(loc1,t_lim);
            count90=histcounts(loc90,t_lim);

            %%% number of windows with:
            % correct mix of elev angles
            out_arr(i,1)=sum( count1==1 & count90==2 );
            % half a scan -- still OK
            out_arr(i,2)=sum( count1==1 & count90==1 );
            % more than one scan
            out_arr(i,3)=sum( count90>2 );
            
        elseif strcmp(scan_type,'short')
            
            % halves of two scans often fit in one window -- code above
            % counts those as good scans
            % shouldn't me much of an issue for long scans, those are more
            % spaced out (and gaps are longer)
            for j=2:length(t_lim)
                
                ind90=find(loc90>=t_lim(j-1) & loc90<t_lim(j));
                ind1=find(loc1>=t_lim(j-1) & loc1<t_lim(j));
                
                %%% number of windows with:
                % correct mix of elev angles
                if length(ind90)==2 && length(ind1)==1 
                    if loc1(ind1)>loc90(ind90(1)) && loc1(ind1)<loc90(ind90(2))
                        % actual complete scan
                        out_arr(i,1)=out_arr(i,1)+1;
                    else
                        % 90deg from 2 different scans, count as half scan
                        out_arr(i,2)=out_arr(i,2)+1;
                    end
                % half a scan    
                elseif length(ind90)==1 && length(ind1)==1
                    out_arr(i,2)=out_arr(i,2)+1;
                % more than one scan
                elseif length(ind90)>2
                    out_arr(i,3)=out_arr(i,3)+1;
                end
                
            end
        end
    end
end


