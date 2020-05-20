function [ P_out, T_out, alt_grid_out ] = construct_PT_apriori( ...
            dates_in_datetime, loc_in, ap_model )
%CONSTRUCT_PT_APRIORI generate pressure, temperaure a priori profiles using
%surface measurements and reanalysis profiles

error('use reanalysis profiles only')

%% get reanalysis profiles
if strcmp(ap_model,'NCEP')

    % retrieve NCEP profiles (daily, all profiles for given day will have
    % identical P, T profiles)
    [ P_prof, T_prof, alt_grid ] = get_NCEP_PT(dates_in_datetime);

elseif strcmp(ap_model,'ERA5')

    % retrieve ERA5 profiles (daily, all profiles for given day will have
    % identical P, T profiles)
    [ P_prof, T_prof, alt_grid ] = get_ERA5_PT(dates_in_datetime,loc_in);

end

    
%% get surface data
[ P_surf, T_surf ] = get_insitu_PT(dates_in_datetime,loc_in);
    
%% combine data
P_out = [ P_surf, P_prof ];
T_out = [ T_surf, T_prof ];
alt_grid_out = [ 10, alt_grid];



plot(P_out,alt_grid_out,'b-'), hold on
%%% leaves bump at merging point
% tmp_fit=fit(alt_grid_out(1:6)',P_out(1:6)','exp1'); % pressure at any altitude level
% P_out(1:6)=tmp_fit(alt_grid_out(1:6)'); % pressure on selected altitude grid

%%% keep surface P fixed: fit surface and one hihger point only
tmp_fit=fit(alt_grid_out(1:9:10)',P_out(1:9:10)','exp1'); % pressure at any altitude level
P_out(1:10)=tmp_fit(alt_grid_out(1:10)'); % pressure on selected altitude grid
plot(P_out,alt_grid_out,'r--')

end

