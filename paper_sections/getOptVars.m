function v = getOptVars()
%GETOPTVARS Return a data structure that contains variables for
%cellulase optimization runs

v.alpha = .2;
v.enzdecayperhour = 1/24;
v.deathrate = 0;
v.dilutionrate = 0; %fraction of media replaced per hour
v.mode = {'Proportional'};%{'Constitutive','Proportional'};
v.layouttype = {'1x1'}; %{'1x1','chemostat'};
v.diff_enz = 1e-6; %Default is 1e-6;
v.costfactor = 1;

v.richmedia = false; %provide a bunch of everything

%%STRATEGY VARIABLES
v.defaultbound = 10; %default exchange upper bound
%v.maxEnzymeRate = 1; %maximum enzyme:biomass ratio by weight for fixed production. 1:1 w/ alpha=1 should mean half of all amino acids go into enzyme
v.alpha = 0.001; %fixed enzyme expression rate relative to growth by mmol:g

%%ENZYME VARIABLES. See https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.4
%values used: wild-type T.ressei on 4-methylumbelliferyl cellobioside. https://www.brenda-enzymes.org/literature.php?e=3.2.1.4&r=737187
v.kcat_cel = 32; %mmol converted / mmol enzyme / second . Source DOI:10.1002/bit.260331112
v.km_cel = 0.7; %Units = mmoles/cm^3. Source PMID:3134347
v.n_enzymes = 2; %number of enzymes in the modified organism from the denHaan experiment

v.enz_weight =  1.5979532e-19 + 7.638479e-20; %grams. Combination of S. fibuligera BGL1 and T reesei EGI/Cel7b
    %for 3-enzyme complex, use weight = 8.5853183e-20
v.atp_per_peptide = 8; %used to calculate energy costs in getEnzymeStoich
v.gtp_per_peptide = 4; %see https://www.ncbi.nlm.nih.gov/books/NBK224633/
%v.enzByWeight = false; %Don't set this to true! Tiny fluxes get rounded
%off!

%%Glucose variables
v.km_glc = 0.01; %millimoles/cm^3.
v.vmax_glc = 1; %10;
v.glcpercellulose = 218.2;%how many units of glucose get released per unit cellulose? 1:1 assumes units are by weight
%Avicel degree of polymerization cite : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3217711/
%v.biomassperglc = 1/(1.447e9); %how many units of biomass are produced for every unit glc consumed?
%weight of E coli / weight of glc molecule = 433fg / (2.992e-22g) = 4.33e-13 / (2.992e-22) = 1.447e9
v.glcuptakerate = -1; %bound for the GLC import reaction

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0; %Diffusion constants
%v.diff_enz = 1e-6; %Default is 1e-6;
v.diff_glc = 1e-6;
%v.enzdecayrate = 0; %degradation of enzyme over time
%v.deathrate = 0;
%v.enzymeperglc = 0.05; %how many units of enzyme are produced for every unit glc consumed?
v.initcellulose = 1;
v.initglc = 1e-3; %needed to bootstrap production
v.initmedia = 1; %amount of all other media added
v.initialpop = 1e-3;
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.initergst = 0;
v.initzymst = 0;
v.initoleate = 0;
v.initpalmitoleate = 0;

%%COMETS EXECUTION VARIABLES
v.maxcycles = 1000;
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 1000000;
v.spaceWidth = 10^(1/3);
v.lograte = 1;