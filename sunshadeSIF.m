%% calculate sun-shade leaf portion and SIFsun and SIFshade

%
clear all
clc

load('SIF760_result.mat','raw_ref_result','wl');
load('hf_barn_2013_env.mat','cloud_ratio','doy');
ref = raw_ref_result(:,2:end);
[xdim_ref,ydim_ref] = size(ref);


for uni_i = 6:6 %1:130
   
    
    
    
    
    
   % 1.Limit the data to the day of interest
   lb       = uni_i-1.0+170.;
   ub       = uni_i+170;
   sub_temp = raw_ref_result(:,1) >= lb & raw_ref_result(:,1) < ub;
   sub_temp1= doy >= lb & doy <= ub;
   
   ref1     = raw_ref_result(sub_temp,wl>=680 & wl < 690);
   [xdim,ydim] = size(ref1);
   diffuse  = cloud_ratio(sub_temp1);
   %refsg    = sgolayfilt(ref1,3,7);
   ref30    = zeros(48,ydim+1);
   ref30(:,1) = 0:0.5:23.5;
   time_ref   = raw_ref_result(sub_temp,1) - lb;
   
   for jj = 1:48
     if ref30(jj,1)<time_ref(1) || ref30(jj,1)>time_ref(end)
         ref30(jj,2:end) = 0.00;
     else
         spect_sub       = knnsearch(time_ref,jj*0.5,'K',1);
         ref30(jj,2:end) = ref1(spect_sub,:);
     end
    
    ydata = ref30(jj,2:end);
    xdata(1,:) = repmat(diffuse(jj),1,length(ydata));
    xdata(2,:) = (wl(wl>=680 & wl < 690))'; 

    x0 = [0.5,1,1,0.1,0.1,0.1,0.1];
    lb = [0,0,0,-inf,-inf,-inf,-inf];
    ub = [1,1,1,inf,inf,inf,inf];
    options = optimset('TolX',1e-10,'TolFun',1e-10,'Display','off');
    [x,SSresid] = lsqcurvefit(@linearmix,x0,xdata,ydata,lb,ub,options);
    
    SStotal = (length(ydata)-1) * var(ydata);
    rsq(jj) = 1 - SSresid/SStotal;
    rmse(jj) = (SSresid/(length(ydata)))^(0.5);
    shade(jj) = x(1);
     
     
     
   end
   
   
   
   
   
   
   
   
   

end