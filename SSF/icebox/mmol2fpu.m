function fpu = mmol2fpu(con,kcat,km)
%MMOL2FPU(concentration,kcat,km) Convert the given amount of enzyme in mmol to Filter Paper Units
%   See for details: https://www.researchgate.net/profile/Lukasz_Wajda/post/Enzyme_activity_and_its_unit_conversion-can_anyone_help/attachment/59d62e49c49f478072e9ef17/AS%3A273573757816834%401442236474228/download/Measurement+of+cellulase+activities.pdf
%   Assumes the reaction product is glucose, or similarly weighted reducing
%   sugars.
%   Assumes units for kcat are in s^-1

%See https://www.nrel.gov/docs/gen/fy08/42628.pdf for instructions

%Step 1: find the dilution factor s.t. 0.5mL of enzyme solution yields 2mg
%glc from 50mg cellulose in 1h

%use formula for Michaelis-Menten with inhibition

%ApparentK_M = K_M(1 + [I]/K_i)
%use K_i = 90 for D-glucose inhibition of celliobiohyrdolase: http://brenda-enzymes.org/enzyme.php?ecno=3.2.1.4#Ki%20VALUE%20[mM]
%   alternative: should we use K_i = 1.6 for cellobiose inhibition instead?
ki = 90;
tsteps = 3600;
options = odeset('NonNegative',1);
tspan = 1:tsteps;
threshold = 0.01; %how close do we want to get to 2mg?


glctarget = 2/0.18; %umoles
glctarget = glctarget / 1000; %millimoles
glctargetupper = glctarget * (1 + threshold);
glctargetlower = glctarget * (1 - threshold);

%initcel = glctarget * 25; %the target is 4% conversion
initcel = 1;

%iterate dilutions until you get close on either end of 2mg
df = 1;
lower = 0;
upper = 0;
lowres = glctargetlower;
upres = glctargetupper;
while (lower == 0 || upper == 0)
    enzcon = con / (df); %is this doing it right? REVISIT THIS LINE
    [t,y] = ode45(@(t,x) dydt(x(1),x(2),enzcon,kcat,km,ki),tspan,[initcel,0],options);
    res = y(2,end)
    if res < lowres
        df = df * .9;
    elseif res > upres
            df = df * 1.1;
    elseif res > lowres & res < glctarget
        lower = df;
        lowres = res;
        df = df * .9;
    elseif res < upres & res > glctarget
        upper = df;
        upres = res;
        df = df * 1.1;
    end
end

%find res=2 on the line bt the two closest points using log of dilution
%factor on the y axis
ll = log10(lower);
lu = log10(upper);

coeffs = polyfit([lower, upper],[ll, lu], 1);
%y = mx+b so when x = 2, y = 2slope + int
logd = 2 * coeffs(2) + coeffs(1);
d = 10^logd;

%Step X: FPU = (0.37 / dilution factor) units/mL
fpu = 0.37 / d;

end

function dydt = dydt(sub,glc,enz,kcat,km,ki)
kma = km * (1 + (glc / ki));
v = kcat * enz * sub / (kma + sub);
ratio = 1;
dydt = [ -v/ratio; %substrate   %INCLUDE CONVERSION RATIO BT CELLULOSE AND GLC! 
    v]; %glc

end