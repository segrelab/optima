function mmol = fpu2mmol(fpu,kcat,km)
%FPU2MMOL Convert the given amount of enzyme in Filter Paper Units to mmol

%1 FPU releases 2mg of glucose from 50mg of filter paper in 1 hour
%So we have to find the concentration that releases 2FPU mg in 1h

ki = 90;
tsteps = 3600;
options = odeset('NonNegative',1);
tspan = 1:tsteps;
threshold = 0.01; %how close do we want to get to 2mg?

df = 1;
lower = 0;
upper = 0;
lowres = targetlower;
upres = targetupper;

enzcon = fpu;
while (lower == 0 || upper == 0)
    
    [t,y] = ode45(@(t,x) dydt(t,x,enzcon,kcat,km,ki),tspan,[initcel,0],options);

end

function dydt = dydt(t,x,enz,kcat,km,ki)
 %sub and glc are in mg, enz is in mmol
sub = x(1,end);
glc = x(2,end);
%kcat is in s^-1. Convert it to mg/s
kcat_mg = kcat * 180.16; %mg Glc / s
km_mg = km * 180.16;

glc_con = glc / 180.16;

kma = km * (1 + (glc_con / ki));
kma_mg = kma * 180.16;
v = kcat_mg * enz * sub / (kma_mg + sub);
ratio = 1;
dydt = [ -v * ratio; %substrate 
    v]; %glc

end