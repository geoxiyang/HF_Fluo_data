% Extrapolate 15min data

% First we extrapolate the 15min measurements to the time when SIF
% measurments were made

[data_ems,txt,raw] = xlsread('F:\01.ResearchProject\9.Fluorescence\10.GPP_Proc\hf2013gpp_v2\gapfill_gpp.xlsx');
load('SIF760_result.mat','final_result_time');

gpp = data_ems(:,144);
doy = data_ems(:,149);

sub0 = doy>0 & gpp>0;

gpp = gpp(sub0);
doy = doy(sub0);

gppfill = NaN(15681,1);

% par15doy = data_15min(:,8);
% par15 = data_15min(:,6);
% vpd15 = data_15min(:,11);
% airt15 = data_15min(:,9);

% par15fill = zeros(15681,1);
% vpd15fill = zeros(15681,1);
% airt15fill = zeros(15681,1);

for i=1:15681
   if isnan(final_result_time(i,1))
%         par15fill(i,1) = NaN;
%         vpd15fill(i,1) = NaN; 
%         airt15fill(i,1) = NaN; 
          gppfill(i,1) = NaN;
        continue
   end
    
   [subs,dists] = knnsearch(doy,final_result_time(i,1),'K',2);
   
   weight1 = 1/abs(doy(subs(1))-final_result_time(i,1));
   weight2 = 1/abs(doy(subs(2))-final_result_time(i,1));
   
   gppfill(i,1) = (weight1*gpp(subs(1))+weight2*gpp(subs(2)))/(weight1+weight2);

%    par15fill(i,1) = (weight1*par15(subs(1))+weight2*par15(subs(2)))/(weight1+weight2);
%    vpd15fill(i,1) = (weight1*vpd15(subs(1))+weight2*vpd15(subs(2)))/(weight1+weight2);
%    airt15fill(i,1) = (weight1*airt15(subs(1))+weight2*airt15(subs(2)))/(weight1+weight2);
    
end

save('SIF760_result.mat','gppfill','-append');

clearvars