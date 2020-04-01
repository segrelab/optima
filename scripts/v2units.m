function res = v2units(delta,timestep)
%return the number of Enzyme Units based on the given rate of production 

%Enzyme Units are in umol/min
%delta is in millimoles
d = abs(delta * 1000); %umol/timestep
%timestep is in hours
ts = timestep * 60; %now timestep is in minutes

res = d / ts; %umol/min / min
end

