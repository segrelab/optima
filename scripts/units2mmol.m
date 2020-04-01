function res = units2mmol(units,kcat,km,sub)
%UNITS2MMOL convert the given value of enzyme units to mmol
%   Assume one Unit is the amount of enzyme that produces 1 umol product per minute 

% % 
% % enz_umol = 1000 * enz_mmol;
% % sub_umol = 1000 * sub_mmol;
% % km_umol = 1000 * km;
% % 
% % vcat_umol = enz_umol * kcat * sub_umol / (sub_umol + km_umol)
% % umolPerMin = 60 * vcat_umol

%Units = how many umol of product are released per minute
vcat_umol = units/60; %vcat in seconds
%vcat_umol = enz_umol * kcat * sub_umol / (sub_umol + km_umol)
%-> vcat_umol * (sub_umol + km_umol) / (kcat * sub_umol) = enz_umol
sub_umol = 1000 * sub;
km_umol = 1000 * km;
enz_umol = vcat_umol * (sub_umol + km_umol) / (kcat * sub_umol);
res = enz_umol / 1000;
end

