% to plot maxdoas results for Pandora

%% load and filter data

load('/home/kristof/work/PANDORA/profiling_test/QDOAS_output/Pandora_104_vis_2018.mat')

% data.Elevviewingangle=round(data.Elevviewingangle);
% 
% data(data.NO2_VisSlColno2>1e20,:)=[];

%% control variables
x_date=1;
plot_o4=1;

msize=16;

day1 = 64; 
day2 = 152;


%% set up the plot

% elevation angles in the file (might change from year to year)
elevs=unique(data.Elevviewingangle)';
elevs_expected=[1,2,3,5,8,10,15,20,30,40,50];

if ~isempty(setdiff(elevs,elevs_expected))
    
    warning('Unexpected elevation angles')
    elevs=intersect(elevs,elevs_expected);
    
end

% colors for plotting
colors={'r.','g.','b.','y.','m.','c.','k.','ro','go','bo','yo','mo','co','ko'};

% empty legend cell
legends=cell(1,length(elevs));

% initialize figure
figure(1)
% figure('Position', [100, 100, 1410, 1050])
set(gcf, 'Position', [100, 100, 900, 600]);

if plot_o4, ax_tg=subplot(2,1,1); end
box on
hold on;
grid on
% grid minor
% ax.GridAlpha=0.4;
% grid minor

% select plotting axis
if x_date % datetime
    x_arr=data.DateTime;
else % DOY (-1 for fractional date)
    x_arr=data.Fractionalday-1;
end

%% loop over elevations
for i=1:length(elevs)
   
    % find indices with current elevation
    ind=data.Elevviewingangle==elevs(i);
         
    % plot
    if ~isempty( strfind(colors{i},'o') ), msize_tmp=4; else msize_tmp=msize; end
    plot(x_arr(ind),data.NO2_VisSlColno2(ind),colors{i},'markersize', msize_tmp);     
        
    % create legend entry
    legends{i}=[num2str(elevs(i)) '\circ'];
    
end 

% ylim([-0.5e14 8e14]);
% xlabel('Fractional day')
ylabel('NO_2 dSCD (mol/cm^2)')
legend(legends,'location','northeast','Orientation','horizontal')

%%%%%%%%%%%%%
if plot_o4
    
    ax_o4=subplot(2,1,2);
    box on
    hold on;
    grid on
    % grid minor
    % ax.GridAlpha=0.4;
    % grid minor

    for i=1:length(elevs)

        % find indices with current elevation
        ind=data.Elevviewingangle==elevs(i);

        % plot
        if ~isempty( strfind(colors{i},'o') ), msize_tmp=4; else msize_tmp=msize; end
        plot(x_arr(ind),data.NO2_VisSlColo4(ind),colors{i},'markersize', msize_tmp);     

        % create legend entry
        legends{i}=[num2str(elevs(i)) '\circ'];

    end 
    
%     ylim([-1e3 10e3]);
    % xlabel('Fractional day')
    ylabel('O_4 DSCD (mol/cm^2)')
    legend(legends,'location','northeast','Orientation','horizontal')
    
end

try linkaxes([ax_tg,ax_o3],'x'); end
try linkaxes([ax_tg,ax_o4],'x'); end


% set font on plots
set(findall(gcf,'-property','FontSize'),'FontSize',12)
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman') 
% 
% f=gcf; 
% figpos=getpixelposition(f); 
% resolution=get(0,'ScreenPixelsPerInch'); 
% set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); 
% path='/home/kristof/work/summer_school/poster/'; 
% name='BrO_march_19-22'; 
% print(f,fullfile(path,name),'-dpng','-r300','-opengl') %save file

