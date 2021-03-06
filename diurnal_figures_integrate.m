%% Diurnal pattern of SIF and GPP averaged for a certain period of time
% For example, for a month, or ten days

% Xi Yang, geoxiyang@gmail.com
% History:
% Aug.10, 2014: v1.0 Read data, calculate the average

%% 1. Read the data

% clc, clear all
datapath = '/Volumes/XiYangResearch/Projects/9.Fluorescence/11.Matlab_data/';
load([datapath,'SIF760daily.mat'],'halfhourly_result','raw_final_result');
load([datapath,'HF_2013_GPP.mat']); 
load([datapath,'hf_barn_2013_env.mat'],'apar')

%% 2. Set the integration window

% types = {'monthly' 'tenday' 'halfmonthly'};
% selected_type = cellstr(types{1});                                           %Choose the integration type
% 
% switch selected_type
%     case 'monthly'
%         n_period = ceil(halfhourly_result(end,1)/30) - ceil(halfhourly_result(1,1)/30) + 1;
%     case 'tenday'
%         
%     case 'halfmonthly'
% end

        n_period = ceil(halfhourly_result(end,1)/30) - ceil(halfhourly_result(1,1)/30) + 1;


%% 3. Calculate mean daily pattern for each period of time
%     and make figures

font_size = 24;


for ii = 1:n_period %1:n_period
        
        month_str = {'June','July','August','September','October'};
        opath1     = ['/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/' month_str{ii} '.png'];
    
        lb = datenum(2013,ceil(halfhourly_result(1,1)/30)+ii-1,1) - datenum(2013,1,1) + 1;
        if ii == n_period
            ub = ceil(halfhourly_result(end,1));
        else
            ub = datenum(2013,ceil(halfhourly_result(1,1)/30)+ii,1)   - datenum(2013,1,1) + 1;
        end
        hour      = 24.0 * (halfhourly_result(1:48,1) - halfhourly_result(1,1));
        temp_sif  = halfhourly_result(halfhourly_result(:,1)>= lb & halfhourly_result(:,1)<ub,2);
        temp_gpp  = gpp_raw(doy>= lb & doy <ub);
        temp_apar = apar(doy>= lb & doy <ub);
        temp_fday = doy(doy>= lb);
        n_days    = ub - temp_fday(1);
        
        % Reshape sif so that data from each day occupy one column
        sif_cube = reshape(temp_sif,48,n_days);
        gpp_cube = reshape(temp_gpp,48,n_days);
        apar_cube= reshape(temp_apar,48,n_days);
        sif_cube(sif_cube <= 0.0 | sif_cube >= 4.0) = NaN;
        apar_cube(apar_cube<=0.0) = NaN;
        gpp_cube(gpp_cube<=0.0) = NaN;
        % Calculate half-hourly average
        sif_mean = nanmean(sif_cube,2);
        gpp_mean = nanmean(gpp_cube,2);
        apar_mean= nanmean(apar_cube,2);
        sif_sd   = nanstd(sif_cube,0,2);
        gpp_sd   = nanstd(gpp_cube,0,2);
        apar_sd  = nanstd(apar_cube,0,2);
        
        sif_mean(isnan(sif_mean)) = 0.00;
        sif_sd(isnan(sif_sd))     = 0.00;
        apar_mean(isnan(apar_mean)) = 0.00;
        apar_sd(isnan(apar_sd))     = 0.00;        
        
        
%         sify_cube = sif_cube./apar_cube;
%         lue_cube  = gpp_cube./apar_cube;
%        
%         
%         sify_mean = nanmean(sify_cube,2);
%         lue_mean  = nanmean(lue_cube,2);
%         lue_mean = lue_mean(20:30,1);
%         sify_mean = sify_mean(20:30,1);
%         lue_mean(lue_mean > 0.1) = NaN;
%         sify_mean(sify_mean > 0.1) = NaN;
%         
%         corr(sify_mean(~isnan(sify_mean)),lue_mean(~isnan(sify_mean)))^2
%         
%         plot(1:1:48,lue_mean,'ro',1:1:48,sify_mean*10,'ko')
%         ylim([0 0.1])
%        Make graph
        
        % Plot SIF
        figure('units','normalized','position',[0 0 1 1])
        h1     = gca;
        set(h1,'Position',h1.Position - [0.05 0 0 0],...
               'YColor','r',...
               'XAxisLocation','bottom',... 
               'YAxisLocation','left',...
               'FontSize',font_size,...
               'FontName','Whitney',...
               'ylim',[0,3],...
               'xlim',[0,24],...
               'NextPlot','add');
        [hb1,p1] = boundedline(hour, sif_mean, sif_sd, 'ro-', h1,'alpha','transparency', 0.1);
        ylabel(h1, 'SIF(mw/m^{2}/sr/nm)','FontName','Whitney');
        l1 = xlabel('Hours','FontSize',font_size,'FontName','Whitney');
        h1_pos = h1.Position;
        set(hb1,'MarkerSize',16);
        
        pos = [0.5,2.9];
        text(pos(1)+23,pos(2)-0.2,month_str{ii},'FontSize',font_size+16.0,'FontName','Whitney','HorizontalAlignment','right');
        text(pos(1),pos(2)-0.2,['R^{2}=',num2str(corr(gpp_mean,sif_mean)^2,'% 5.2f')],'FontSize',font_size+16.0,'FontName','Whitney','Color','b');

        %         % Plot APAR
        
        h3     = axes('Position',get(h1,'Position')+[0 0 (1/24)*h1_pos(1,3) 0],...
                       'XAxisLocation','top',...
                       'XTick',[],...
                       'XTickLabel',[],...
                       'YAxisLocation','right',...                       
                       'Color','none',...
                       'XColor',get(gcf,'Color'),'YColor','m',...
                       'FontSize',font_size,...
                       'FontName','Whitney-Book',...
                       'ylim',[0,2000],...
                       'xlim',[0,25]    );
        h3.YTick = [0,500,1000,1500,2000];
        ylabel(h3, 'GPP or APAR(umol/m^{2}/second)','FontName','Whitney','Color','k');
        [hb3,p3] = boundedline(hour, apar_mean, apar_sd, 'mo-', h3,'alpha','transparency', 0.1);
        text(pos(1),pos(2)+1500,['R^{2}=',num2str(corr(apar_mean,sif_mean)^2,'% 5.2f')],'FontSize',font_size+16.0,'FontName','Whitney','Color','m');
        set(hb3,'MarkerSize',16);
        
        % Plot GPP
        h2     = axes( 'Position',h1.Position,...
                       'XAxisLocation','top',...
                       'XTick',[],...
                       'XTickLabel',[],...
                       'YAxisLocation','right',...
                       'Color','none',...
                       'XColor','k','YColor','b',...
                       'FontSize',font_size,...
                       'FontName','Whitney-Book',...
                       'ylim',[-5,40],...
                       'xlim',[0,24]    );
        h2.YTick = [0,10,20,30,40];        
        %ylabel(h2, 'GPP(umol/m^{2}/second)','FontName','Whitney');
        [hb2,p2] = boundedline(hour, gpp_mean, gpp_sd, 'bo-', h2,'alpha','transparency', 0.1);
        set(hb2,'MarkerSize',16);

        a = get(gcf,'OuterPosition');
        set(gcf,'OuterPosition',[a(1),a(2),a(3)+0.2,a(4)+0.2*a(4)/a(3)]);

        
        set(gcf,'paperPositionMode','manual','PaperPosition',[0,0,14,8]) % make the print as big as the figure
        print(gcf, '-dpng','-r300', opath1);
        close(gcf);
        
        
        
        
end

