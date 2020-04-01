%build the "celtab" table with information for multiple simulations
function celtab = buildCelOptDataTable()

%% INPUT CONFIGURATION
%mode = {'p' 'c'}; %proportional or constitutive
mode = {'p'};
%alphas = [0.5:0.05:0.8 0.9:0.1:1.6];
%alphas = [0.01 0.015 0.02 0.025 0.03 0.05 0.075 0.1 0.15 0.2 0.3 0.5 0.75 1 1.5 2];
%alphas = [0.01 0.025 0.05 0.075 0.1 0.25 0.5 1 1.5 2];
alphas = 0.05;
%cel_amts = [0.1 0.2 0.35 0.5 0.75 1 1.5 2 2.5 3 4];
%cel_amts = [0.1 0.5 1 1.5 2];
%cel_amts = [0.5 0.75 1 1.25 1.5];
cel_amts = [1.58];

%deathrates = [0];
%deathrates = [0.001 0.005 0.01 0.05 0.1 0.2];%units here PER HOUR, not per timestep
%deathrates = 0.01:0.01:0.05;
%deathrates = [0.001 0.0025 0.005 0.0075];
deathrates = 0;
%Inputs to COMETS should be death rate per HOUR, since the code is 
  %deltaBiomass[i] -= cParams.getDeathRate() * biomass[i] * cParams.getTimeStep();
  %biomass[i] += deltaBiomass[i];
%v.deathrate will report per timestep, v.deathperhour will report this
%value

%enzdecayperhour = 0.1; %this row is in fraction/hour
%convert fraction/hour to fraction/second, which COMETS accepts
%enzdecaypersec = enzdecayperhour / 3600;
enzdecaypersec = 0;
enzdecayperhour = enzdecaypersec * 3600;
%enzhalflife = log(2)./enzdecaypersec;

%% STRATEGY VARIABLES
v.defaultbound = 10; %default exchange upper bound
v.maxEnzymeRate = 1; %maximum enzyme:biomass ratio by weight for fixed production. 1:1 w/ alpha=1 should mean half of all amino acids go into enzyme
v.alpha = 0.1; %fixed percentage of the maximum enzyme rate being expressed


%% ENZYME VARIABLES. See https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.4
%values used: wild-type T.ressei on 4-methylumbelliferyl cellobioside. https://www.brenda-enzymes.org/literature.php?e=3.2.1.4&r=737187
%v.kcat_cel = 0.01; %mmol converted / mmol enzyme / second  %TODO: Confirm units
v.kcat_cel = 0.1 * 2.992e-22 / 4.17775e-20; %mmol converted / mmol enzyme / second  %TODO: Confirm units
v.km_cel = 0.125; %Units = mmol
v.enz_weight = 8.5853183e-20; %grams
v.enzByWeight = false; %if true enzyme units are grams. Otherwise they're molecules

%% Glucose variables
v.km_glc = 0.01; %millimoles.
v.vmax_glc = 10;
v.glcpercellulose = 218.2;%how many units of glucose get released per unit cellulose? 1:1 assumes units are by weight
%Avicel degree of polymerization cite : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3217711/
%v.biomassperglc = 1/(1.447e9); %how many units of biomass are produced for every unit glc consumed?
v.biomassperglc = 6.9109e-4;
%weight of E coli / weight of glc molecule = 433fg / (2.992e-22g) = 4.33e-13 / (2.992e-22) = 1.447e9
v.glcuptakerate = -10; %bound for the GLC import reaction
%This should be higher than required for max enzyme production


%Enzyme variables v 2- set high to check that everything works OK
%v.kat = 100;
%v.km = 0.0001;

%% WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0; %Diffusion constants
v.diff_enz = 0;%1e-6; %Default is 1e-6;
v.diff_glc = 0;%1e-6;
%v.enzdecayrate = 0; %degradation of enzyme over time
%v.deathrate = 0;
%v.enzymeperglc = 0.05; %how many units of enzyme are produced for every unit glc consumed?
v.initcellulose = 1;
v.initglc = 0.001; %needed to bootstrap production
v.initmedia = 10000; %amount of all other media added
v.initialpop = 1e-4;

%% COMETS EXECUTION VARIABLES
v.maxcycles = 2000;
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second. CAUTION: Make sure this agrees with the units in kcat.
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.maxspacebiomass = 1000;
v.spaceWidth = 1; %1cm, so v=1mL and most media are 1 millimolar
v.lograte = 10;

v_default = v;

%% DEFAULT MODEL
%model_default = load('C:\sync\biomes\models\ecoli_core_model_norm.mat'); %model with metabolite names normalized
model_default = load('C:\sync\biomes\models\iJO1366 (1).mat'); %model with metabolite names normalized
model_default = model_default.iJO1366;

%% EXECUTION
%varnames = {'cellulose','alpha','model'};
%celtab = array2table(nan(length(cel_amts) * length(alphas),length(varnames)),'VariableNames',varnames);
nrows = length(mode) * length(cel_amts) * length(alphas) * length(deathrates) * length(enzdecaypersec);
c = zeros(1,nrows);
a = zeros(1,nrows);
m = cell(1,nrows);
idx = 1;
for h = 1:length(mode)
    v = v_default;
    v.mode = mode{h};
    for i = 1:length(cel_amts)
        v.initcellulose = cel_amts(i);
        for j=1:length(alphas)
            alpha = alphas(j);
            v.alpha = alpha;
            for k = 1:length(deathrates)
                for l = 1:length(enzdecaypersec)
                    v.enzdecaypersec = enzdecaypersec(l);
                    v.enzdecayperhour = enzdecayperhour(l);
                    %v.enzhalflife = enzhalflife(l);
                    v.deathrate = deathrates(k) / v.timestep;
                    v.deathperhour = deathrates(k);
                    mname = [v.mode '-a' num2str(alpha) '-c' num2str(v.initcellulose) '-d' num2str(v.deathrate)];
                    model = buildEcoliModel(v,model_default,mname);
                    %model = executeModel(model);
                    c(idx) = cel_amts(i);
                    a(idx) = alphas(j);
                    m{idx} = model;
                    idx = idx + 1;
                end
            end
        end
    end
end

celtab = table(c',a',m','VariableNames',{'cellulose_init','alpha','model'});
end