function res = mmol2units(mmol,kcat,km,sub)
%MMOL2UNITS Convert the given amount of enzyme to Enzyme Units
%Assume one Unit is the amount of enzyme that produces 1 umol product per minute 
%Equivalent to asking how many umol of product per minute does this amount
%produce? (assuming [substrate] >> [enzyme])
%kcat is in seconds^-1, km and sub are in millimoles

enz_umol = 1000 * mmol;
sub_umol = 1000 * sub;
km_umol = 1000 * km;

vcat_umol = enz_umol * kcat * sub_umol / (sub_umol + km_umol);
res = 60 * vcat_umol;

end

