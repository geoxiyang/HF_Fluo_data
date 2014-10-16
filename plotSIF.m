% plot ChlF data


sub = rsq >0.9 & rmse <1 & f760 >0;
final_result_sub = [time(sub), f760(sub)];

final_result_time = sortrows(final_result_sub,1);
uni_doy = unique(int16(final_result_time(:,1)));
final_result_time((final_result_time(:,2) <0 | final_result_time(:,2) >5),2) = NaN;
doy_size = size(uni_doy);
f760_daily = NaN(doy_size);

for uni_i = 1:doy_size(1)
   
   lb = double(uni_doy(uni_i))+0.;
   ub = double(uni_doy(uni_i))+1.;
   
   f760_daily(uni_i) = sum(final_result_time(final_result_time(:,1) >= lb ...
       & final_result_time(:,1) <= ub, 2));
end

final_result_daily = [uni_doy,f760_daily];

%plot(uni_doy,f760_daily,'g.');

[AX,H1,H2] = plotyy(uni_doy,f760_daily,DOY_HF_2013,Precip_HF_2013);
