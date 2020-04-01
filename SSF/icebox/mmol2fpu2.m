function fpu = mmol2fpu2(con,kcat,km)
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
tsteps = 60;
options = odeset('NonNegative',1);
tspan = 1:tsteps;
threshold = 0.01; %how close in % do we want to get to the target product?


%find how many Units (U) per volume of enzyme we have
%1U = amount of enzyme that will hydrolyze 1umol of glc in 1 minute

%C = A * m / V, where C = concentration, A is activity in Units per mass of
%material, m is mass of material used to make the resultant concentration,
%and V is the final volume