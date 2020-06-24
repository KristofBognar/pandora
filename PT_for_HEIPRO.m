function PT_for_HEIPRO( yr_in, loc_str )
%PT_FOR_HEIPRO Summary of this function goes here
%   Detailed explanation goes here

savedir='/home/kristof/work/PANDORA/profiling/retrieval_input/HEIPRO_input/ERA5_files/';

t1 = datetime(yr_in,1,1,12,0,0);
t2 = datetime(yr_in,12,31,12,0,0);
daily_times = t1:t2;

% retrieve daily ERA5 profiles
[ hPa_all, K_all, asl_km_orig ] = get_ERA5_PT(daily_times,loc_str);

% convert alt grid to km (P is lready in hPa and T is in K)
asl_km_orig=asl_km_orig/1000;

for i=1:length(daily_times)
    % interpolate to 100m grid
    grid=[0.1:0.1:30];
    hPa=interp1(asl_km_orig,hPa_all(i,:),grid,'linear','extrap')';
    K=interp1(asl_km_orig,K_all(i,:),grid,'linear','extrap')';
    asl_km=grid';

    % create table from data
    radiosonde=table(asl_km,hPa,K);

    % create filename and write file
    savename=[savedir 'ERA5_' loc_str '_' datestr(daily_times(i),'yymmdd') '.dat'];
    writetable(radiosonde,savename,'Delimiter',',');
    
end

end

