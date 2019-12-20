function merge_pandora_dscds( p_num, year, lambda )
%reads in and merges Pandora dSCDs for given instrument and year (lambda is
%optional, default is 'vis')


if nargin==2, lambda='vis'; end

%% find files

cd('/home/kristof/work/PANDORA/profiling_test/QDOAS_output')

% make list of all files
tmp = dir(['Pandora_' num2str(p_num) '_' lambda '_' num2str(year) '_*']); 
f_list = {tmp.name}; % cell array of file names

% remove backup files
tmp=find_in_cell(f_list,'~');
if ~isempty(tmp), f_list(tmp)=[]; end

%% read files

data=[];

for i=1:length(f_list)
    disp(['Reading ' f_list{i}])
    tmp=read_maxdoas_table_pandora('./',f_list{i});
    data=[data;tmp];
end

%% save yearly dSCD file

    save(['Pandora_' num2str(p_num) '_' lambda '_' num2str(year) '.mat'],'data');
    

end

