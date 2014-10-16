%% Diurnal pattern of SIF and GPP averaged for a certain period of time
% For example, for a month, or ten days

% Xi Yang, geoxiyang@gmail.com
% History:
% Aug.10, 2014: v1.0 Read data, calculate the average

%% 1. Read the data

clc, clear all

load('SIF760daily.mat','halfhourly_result','raw_final_result');
load('HF_2013_GPP.mat'); 

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
for ii = 1:n_period
        
        month_str = {'June','July','August','September','October'};
        opath1     = ['/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/' month_str{ii} '.png'];
        opath2     = ['/Volumes/XiYangResearch/Projects/9.Fluorescence/4.JPG/' month_str{ii} '.eps'];
    
        lb = datenum(2013,ceil(halfhourly_result(1,1)/30)+ii-1,1) - datenum(2013,1,1) + 1;
        if ii == n_period
            ub = ceil(halfhourly_result(end,1));
        else
            ub = datenum(2013,ceil(halfhourly_result(1,1)/30)+ii,1)   - datenum(2013,1,1) + 1;
        end
        hour      = 24.0 * (halfhourly_result(1:48,1) - halfhourly_result(1,1));
        temp_sif  = halfhourly_result(halfhourly_result(:,1)>= lb & halfhourly_result(:,1)<ub,2);
        temp_gpp  = gpp_raw(doy>= lb & doy <ub);
        temp_fday = doy(doy>= lb);
        n_days    = ub - temp_fday(1);
        
        % Reshape sif so that data from each day occupy one column
        sif_cube = reshape(temp_sif,48,n_days);
        gpp_cube = reshape(temp_gpp,48,n_days);
        sif_cube(sif_cube <= 0.0 | sif_cube >= 4.0) = NaN;
        % Calculate half-hourly average
        sif_mean = nanmean(sif_cube,2);
        gpp_mean = mean(gpp_cube,2);
        sif_sd   = nanstd(sif_cube,0,2);
        gpp_sd   = std(gpp_cube,0,2);
        
        sif_mean(isnan(sif_mean)) = 0.00;
        sif_sd(isnan(sif_sd))     = 0.00;
        
        % Corrleation
        month_str(ii)
        corr(gpp_mean,sif_mean)^2
        
        % Make graph

        figure
        h1     = gca;
        set(h1,'YColor','r',...
               'YAxisLocation','left',...
               'FontSize',30,...
               'FontName','Whitney',...
               'ylim',[0,5],...
               'xlim',[0,24]);
        boundedline(hour, sif_mean, sif_sd, 'ro-', h1,'transparency', 0.1);
        ylabel(h1, 'SIF(mw/m^{2}/sr/nm)','FontName','Whitney');

        h2     = axes('Position',get(h1,'Position'),...
                       'XAxisLocation','bottom',...
                       'XTick',[],...
                       'XTickLabel',[],...
                       'YAxisLocation','right',...
                       'Color','none',...
                       'XColor','k','YColor','b',...
                       'FontSize',30,...
                       'FontName','Whitney-Book',...
                       'ylim',[-5,40],...
                       'xlim',[0,24]    );
        ylabel(h2, 'GPP(umol/m^{2}/second)','FontName','Whitney');
        l1     = xlabel('Hours','FontSize',30,'FontName','Whitney');
        set(l1,'Position',get(l1,'Position')-[0,1.0,0]);
        
        
        boundedline(hour, gpp_mean, gpp_sd, 'bo-', h2,'transparency', 0.1);
        
        xlim_fig = xlim;
        ylim_fig = ylim;
% 
%         wdth = xlim_fig(2)-xlim_fig(1);
%         ht = ylim_fig(2)-ylim_fig(1);
%         pos = [xlim_fig(1)+0.02*wdth ylim_fig(2)-0.15*ht]
        pos = [0.5,40];
        text(pos(1),pos(2),month_str{ii},'FontSize',30,'FontName','Whitney');
        
        a = get(gcf,'OuterPosition');
        set(gcf,'OuterPosition',[a(1),a(2)+5,a(3)+400,a(4)+300]);
        
        
%         a = get(gcf,'PaperPosition');
%         set(gcf,'paperPositionMode','manual');
%         set(gcf,'PaperPosition',[a(1),a(2)+5,a(3)+25,a(4)+25]);
        
        
        set(gcf,'paperPositionMode','auto') % make the print as big as the figure
        print(gcf, '-dpng','-r300', opath1);
        print(gcf, '-depsc', '-r300', opath2);
        close(gcf);
        
        
        
        
end
