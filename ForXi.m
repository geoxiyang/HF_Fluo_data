%This is the flux partitioning program.  You'll need to define some
%variables before you run this program. There are two different methods to
%partition the flux.  The first uses nightime NEE and a Q10 function to
%separate GEE and ER, and the second uses daytime NEE to fit a light
%response curve and then separates the values afterwards.
%Variables that are needed to run the program are:
%NEEA1=Your NEE that you want to partition
%PARoA1= The PAR or radiation value you want to use
%AirTA1 = The temperature value that you want to use
%For this code you'll end up with 2 GEE values; one derived from the Q10
%response of nighttime ecosystem respiration and the other derived from the
%light response curve (see manuscript).  


%%%%%This is where you define your input variables (you'll also need to
%%%%%filter out low ustar or any other bad data (improper fetch) prior to
%%%%%using this code). You can modify the equations used to fit the data
%%%%%with the Q10 or lightcurve m-files, just remember that if you change
%%%%%them you'll need to change their formulation in this code when you
%%%%%call them up
%DOY = 
%NEEA1 =
%PARoA1 =
%AirTA1 = 


%%%%%This part will find ecosystem respiration with a Q10 value
NEEunfill1 = NEEA1; %This creates a variable that will be manipulated to filter out bad values 
BAD1 = NEEA1<0 | (PARoA1)>10; %This is the filter that takes out nightime values; may need to include ustar filter
NEEunfill1(BAD1) = nan; %This makes "bad" values nan
Daywind = 14; %no of days for your parameterization (THIS CAN BE CHANGED)
Inc = 1; %Time step increment for your parameterization (THIS CAN BE CHANGED, but 1 day creates new parameters for each day) 
binwid = (Daywind); %This sets up the temporal window for your parameterization
start = min(DOY);
last = max(DOY);
bincenter = start:1:(last+1); 
binnum=length(bincenter);
ERfit1 = NEEunfill1 * NaN;
DayQ1 = []; %This is where the values for your parameterization will be stored
Resp01 = []; %This is where the values for your parameterization will be stored; This is Basal Respiration
Q101 = []; %This is where the values for your parameterization will be stored; This is Q10
RQ1=[]; %This is where the values for your parameterization will be stored; This is R2 for regression

for i = 1:1:(binnum);
    focusday = (DOY >= (bincenter(i)-(binwid/2)) & (DOY <= (bincenter(i) + (binwid/2)))); %This sets up window
    temp = AirTA1(focusday & isfinite(NEEunfill1)& isfinite(AirTA1)); %Leaves only values for each window
    ER = NEEunfill1(focusday & isfinite(NEEunfill1)& isfinite(AirTA1));%Leaves only values for each window
    figure(1)
    plot(temp, ER,'.')
    hold on
if length(ER) > 5
        P = nlinfit(temp, ER,'Q10',[1 2]); %THIS IS WHERE YOU CALL THE Q10 FUNCTION AND FIT THE DATA TO IT
        ERfit1(focusday) = P(1).*P(2).^((AirTA1(focusday)-10)/10); %This creates new variable ERfit that is your fitted respriation values
plot(AirTA1(focusday), ERfit1(focusday),'-')
hold off
C = corrcoef(ER,ERfit1(focusday& isfinite(NEEunfill1)& isfinite(AirTA1))); %Find R2 of regression
rsq1 = C(1,2)^2; %Find R2 of regression
RQ1=[RQ1 rsq1];
Day =  nanmean(DOY(focusday));
DayQ1 =[DayQ1 Day]; %This is the day associated with the fit
Resp01 = [Resp01 P(1)]; %This is the basal respiration associated with the fit
Q101 = [Q101 P(2)]; %This is the QQ10 associated with the fit
    end
end

GEEQ10= ERfit1 - NEEA1; %This is your derived photosynthesis for the canopy based on nighttime parameterized Q10


%%%%This will partition the data with a daytime light response curve to get
%%%%respiration and photosynthesis; It is set up similarly than the
%%%%previous case and will use the same window and step size defined
%%%%before.
NEEunfill1 = NEEA1;
NEEunfill1(~BAD1) = nan; %!!!!!!!!!!!!!Be careful here if you use Ustar filter; This only works if ustar data is already filtered out
binwid = (Daywind);
start = min(DOY);
last = max(DOY);
bincenter = start:1:(last+1); 
binnum=length(bincenter);
NEEfit1 = NEEunfill1 * NaN;
GEE2001 = NEEunfill1 * NaN;
GEE1 = NEEunfill1 * NaN;
NEE2001 = NEEunfill1*NaN;
alpha11=[]; %Light use effieincy for light response curve
beta11=[]; %NEE max for light response curve
Resp11=[]; %Ecosystem Respiration for light response curve
DayL1=[]; %Day when light response curve is parameterized
RL1=[]; %R2 for light response curve
PAR1=PARoA1;
for i = 1:(binnum);
    focusday = (DOY >= (bincenter(i)-(binwid/2)) & (DOY <= (bincenter(i) + (binwid/2))));
    par_day = PAR1(focusday & isfinite(NEEunfill1)& isfinite(PAR1));
    nee_day = NEEunfill1(focusday & isfinite(NEEunfill1)& isfinite(PAR1));
    figure(1)
    plot(par_day, nee_day,'.')
    hold on
           
    if length(nee_day) > 10 %This will only fit a curve if there is >10 points in period
        P = nlinfit(par_day, nee_day,'lightcurve',[4 0.02 6]);
        NEEfit1(focusday) = P(1) - (P(2)*P(3).*PAR1(focusday))./(P(2).*PAR1(focusday)+P(3)); %Change
        GEE1(focusday) = P(1) - NEEA1(focusday); %THIS IS YOUR GEE DERIVED FROM LIGHT RESPONSE CURVE
plot(PAR1(focusday), NEEfit1(focusday),'-')
hold off
C = corrcoef(nee_day,NEEfit1(focusday & isfinite(NEEunfill1)& isfinite(PAR1))); %Find R2 of regression
rsq2 = C(1,2)^2;
RL1=[RL1 rsq2];
Day =  nanmean(DOY(focusday));
DayL1 =[DayL1 Day]; %This is the day associated with the fit
 Resp11 = [Resp11 P(1)]; %This is the ecosystem respiration associated with the fit
 alpha11 = [alpha11 P(2)];%This is the light use efficiency associated with the fit
 beta11 = [beta11 P(3)]; %This is the NEEmax associated with the fit

    end
end