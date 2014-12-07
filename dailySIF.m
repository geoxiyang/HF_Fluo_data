% Calculate the daily sum of SIF
% Calculate SIF and VPD at a specific time: 0930,1400

clear all

% load('SIF760_result.mat','final_result_time');
% load('SIF760_2014_result.mat','raw_final_result');
load('hf_2013_svd.mat','final_svd');
load('hf_barn_2013_env.mat','cloud_ratio_daily');

% DOY 170-299: 130 ; 2014: 127-280
daily_raw_result = zeros(130,3);  %154,3
daily_raw_result(:,1) = 170:1:299; %170:1:299;127:1:280

halfhourly_result = zeros(130*48,2);  %154*48,2
halfhourly_result(:,1) = 170:(1/48):(300-1/48);  %127:(1/48):(281-1/48);

% fivemin_result = zeros(130*48*6,2);
% fivemin_result(:,1) = 170:(1/288):(300-1/288);

SIF_0930 = zeros(130,4);    %130
SIF_0930(:,1) = 170:1:299;  %170:1:299 ; 127:1:288

SIF_1330 = zeros(130,4);
SIF_1330(:,1) = 170:1:299;

SIF_1400 = zeros(130,4);
SIF_1400(:,1) = 170:1:299;

SIF_1430 = zeros(130,4);
SIF_1430(:,1) = 170:1:299;

SIF_mean = zeros(130,4);
SIF_mean(:,1) = 170:1:299;

SIF_max = zeros(130,4);
SIF_max(:,1) = 170:1:299;

% VPD_0930 = zeros(130,4);
% SIF_0930(:,1) = 170:1:299;
% 
% VPD_1400 = zeros(130,4);
% VPD_1400(:,1) = 170:1:299;

for uni_i = 1:130 %1:130; 1:154
   
   % Define the time
   lb = uni_i-1.0+170.;  %127
   ub = uni_i+170;
   % step 1: select good days
%    sub_temp = raw_final_result(:,1) >= lb & raw_final_result(:,1) <= ub & raw_final_result(:,3) >=0.90;
%    morning_sub = raw_final_result(:,1) >= lb & raw_final_result(:,1) <= (lb+0.5) & raw_final_result(:,3) >=0.90;
%    afternoon_sub = raw_final_result(:,1) >= (lb+0.5) & raw_final_result(:,1) <= ub & raw_final_result(:,3) >=0.90;
%    if (sum(sub_temp) == 0) | (sum(morning_sub) < 10) | (sum(afternoon_sub) < 10)
%       continue 
%    end
   sub_temp = final_svd(:,1) >= lb & final_svd(:,1) <= ub & final_svd(:,3) <=0.05;
   morning_sub = final_svd(:,1) >= lb & final_svd(:,1) <= (lb+0.5) & final_svd(:,3) <=0.05;
   afternoon_sub = final_svd(:,1) >= (lb+0.5) & final_svd(:,1) <= ub & final_svd(:,3) <=0.05;
   if (sum(sub_temp) == 0) | (sum(morning_sub) < 10) | (sum(afternoon_sub) < 10)
      continue 
   end   
   % step 2: linear-interpolation to every 30 minutes between 6am to 6pm
   
%    temp_time    = raw_final_result(sub_temp, 1);
%    temp_array   = raw_final_result(sub_temp, 2);
   temp_time    = final_svd(sub_temp, 1);
   temp_array   = final_svd(sub_temp, 2);
   %vpd_array    = vpd15fill(sub_temp,1);
   % negative values are assigned to zero
   temp_array(temp_array<0) = 0;
   
   sub_0930 = knnsearch(temp_time,9.5/24+lb,'K',1);
   SIF_0930(uni_i,2) = temp_array(sub_0930);   
   SIF_mean(uni_i,2) = mean(temp_array(temp_array>0));
   SIF_max(uni_i,2) = max(temp_array(temp_array>0));
   sub_1330 = knnsearch(temp_time,13.5/24+lb,'K',1);
   SIF_1330(uni_i,2) = temp_array(sub_1330);
   sub_1400 = knnsearch(temp_time,14/24+lb,'K',1);
   SIF_1400(uni_i,2) = temp_array(sub_1400);
   sub_1430 = knnsearch(temp_time,14.5/24+lb,'K',1);
   SIF_1430(uni_i,2) = temp_array(sub_1430);
   
   % For VPD
%    VPD_0930(uni_i,2) = vpd_array(sub_0930);
%    VPD_1400(uni_i,2) = vpd_array(sub_1400);
%    
%    
   if cloud_ratio_daily(uni_i) < 0.50 
       SIF_0930(uni_i,3) = temp_array(sub_0930);   
       SIF_mean(uni_i,3) = mean(temp_array(temp_array>0));
       SIF_max(uni_i,3)  = max(temp_array(temp_array>0));
       SIF_1330(uni_i,3) = temp_array(sub_1330);
       SIF_1400(uni_i,3) = temp_array(sub_1400);
       SIF_1430(uni_i,3) = temp_array(sub_1430);
%        VPD_0930(uni_i,3) = vpd_array(sub_0930);
%        VPD_1400(uni_i,3) = vpd_array(sub_1400);
   else
       SIF_0930(uni_i,4) = temp_array(sub_0930);   
       SIF_mean(uni_i,4) = mean(temp_array(temp_array>0));
       SIF_max(uni_i,4)  = max(temp_array(temp_array>0));
       SIF_1330(uni_i,4) = temp_array(sub_1330);
       SIF_1400(uni_i,4) = temp_array(sub_1400);
       SIF_1430(uni_i,4) = temp_array(sub_1430);
%        VPD_0930(uni_i,4) = vpd_array(sub_0930);
%        VPD_1400(uni_i,4) = vpd_array(sub_1400);    
   end
%           
    time_sub = find(halfhourly_result(((uni_i-1)*48+1):(uni_i*48),1)<temp_time(1) | halfhourly_result(((uni_i-1)*48+1):(uni_i*48),1)>temp_time(end));
    
   for jj = 1:48
     if (halfhourly_result((uni_i-1)*48+jj,1)<temp_time(1) | halfhourly_result((uni_i-1)*48+jj,1)>temp_time(end))
         halfhourly_result((uni_i-1)*48+jj,2) = 0.00;
     else
         halfhourly_result((uni_i-1)*48+jj,2) = mean(temp_array(temp_time >= halfhourly_result((uni_i-1)*48+jj,1) & temp_time < halfhourly_result((uni_i-1)*48+jj,1)+1/48));
%          [subs,dists] = knnsearch(temp_time,halfhourly_result((uni_i-1)*48+jj,1),'K',2);
%          weight1 = 1/abs(temp_time(subs(1))-halfhourly_result((uni_i-1)*48+jj,1));
%          weight2 = 1/abs(temp_time(subs(2))-halfhourly_result((uni_i-1)*48+jj,1));
%          if weight1 == Inf 
%              halfhourly_result((uni_i-1)*48+jj,2) = temp_array(subs(1));
%          elseif weight2 == Inf
%              halfhourly_result((uni_i-1)*48+jj,2) = temp_array(subs(2));
%          else
%              halfhourly_result((uni_i-1)*48+jj,2) = (weight1*temp_array(subs(1))+weight2*temp_array(subs(2)))/(weight1+weight2);
%          end
     end
   end
%    
%    for kk = 1:288
%      if (fivemin_result((uni_i-1)*288+kk,1)<temp_time(1) | fivemin_result((uni_i-1)*288+kk,1)>temp_time(end))
%          fivemin_result((uni_i-1)*288+kk,2) = 0.00;
%      else
%          [subs_2,dists_2] = knnsearch(temp_time,fivemin_result((uni_i-1)*288+kk,1),'K',2);
%          weight1 = 1/abs(temp_time(subs_2(1))-fivemin_result((uni_i-1)*288+kk,1));
%          weight2 = 1/abs(temp_time(subs_2(2))-fivemin_result((uni_i-1)*288+kk,1));
%          if weight1 == Inf 
%              fivemin_result((uni_i-1)*288+kk,2) = temp_array(subs(1));
%          elseif weight2 == Inf
%              fivemin_result((uni_i-1)*288+kk,2) = temp_array(subs(2));
%          else
%              fivemin_result((uni_i-1)*288+kk,2) = (weight1*temp_array(subs_2(1))+weight2*temp_array(subs_2(2)))/(weight1+weight2);
%          end         
%      end
%    end
%    
%    % step 3: integrate SIFflux to SIFdaily
%    daily_raw_result(uni_i,2) = sum(halfhourly_result(((uni_i-1)*48):(uni_i*48),2)*30*60);
%    daily_raw_result(uni_i,3) = sum(fivemin_result(((uni_i-1)*288):(uni_i*288),2)*5*60);
end


save('SIF760daily_2013_SVD.mat');


