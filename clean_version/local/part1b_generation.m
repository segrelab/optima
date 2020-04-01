%Create the data for "1_prove_system_works.mlx"
function dat1b = part1b_generation(arg)
%This script creates simulations with decaying enyzmes

p.version = 1;
v = arg.v;
model_default = arg.model;

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0;
v.initcellulose = 0;

v.initenzyme = 1;
v.initglc = 0;
v.initmedia = 1; %amount of all other media added
v.initialpop = 0.1;
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.initH = 1000;
v.deathrate = 0;
v.enzdecayperhour = 0;
v.costfactor = 3;

%%COMETS EXECUTION VARIABLES
v.timestep = 1/3600; %1/(60*60); %1/3600 hours = 1 second.
v.maxcycles = 300; %2000/v.timestep;
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 1000;
v.spaceWidth = 10^(1/3);
v.lograte = 1;

v_default = v;

dat1b = table();


%%Default run
layout_default = buildDat1BLayout(model_default,v);
dat1b = runAndLoad(layout_default,dat1b);

%% vary decay rate
v = v_default;
v.enzdecayperhour = 10;
layout = buildDat1BLayout(model_default,v);
dat1b = runAndLoad(layout,dat1b);

v = v_default;
v.enzdecayperhour = 100;
layout = buildDat1BLayout(model_default,v);
dat1b = runAndLoad(layout,dat1b);

end
%% functions
function layout = buildDat1BLayout(model, v)
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'costfactor',v.costfactor,'maintRxnAsObjective',v.maintRxnAsObjective);
enzmodel.v = v;
layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);
layout = setMedia(layout,'enzyme[e]',v.initenzyme);
layout.params.numExRxnSubsteps = 25;
end

function res = loadResult(layout)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]'};
media = parseMediaLog(layout.params.mediaLogName,medianames);

res.timestep = unique(media.t);
res.t = res.timestep * layout.models{1}.v.timestep; %time in hours

res.glc_amt = media.amt(stridx('glc-D[e]',media.metname));
res.enzyme_amt = media.amt(stridx('enzyme[e]',media.metname));
res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));

res.glc_delta = calcDelta(res.glc_amt);
res.enzyme_delta = calcDelta(res.enzyme_amt);
res.cellulose_delta = calcDelta(res.cellulose_amt);
res.etoh_delta = calcDelta(res.etoh_amt);

biomass = parseBiomassLog(layout.params.biomassLogName);
res.biomass= biomass.biomass;
res.biomass_delta = calcDelta(res.biomass);
end

function t = applyResult(layout, t)
%add the layout and its results to the table
res = loadResult(layout);
tab = table;

tab.layout = {layout};
tab.t = {res.t};
tab.timestep = {res.timestep};
tab.glc_amt = {res.glc_amt};
tab.enzyme_amt = {res.enzyme_amt};
tab.cellulose_amt = {res.cellulose_amt};
tab.etoh_amt = {res.etoh_amt};
tab.biomass = {res.biomass};
tab.glc_delta = {res.glc_delta};
tab.enzyme_delta = {res.enzyme_delta};
tab.cellulose_delta = {res.cellulose_delta};
tab.etoh_delta = {res.etoh_delta};
tab.biomass_delta = {res.biomass_delta};

tab.alpha = layout.models{1}.v.alpha;
tab.decayrate = layout.models{1}.v.enzdecayperhour;
tab.costfactor = layout.models{1}.v.costfactor;

if size(t,1) > 0
    t = [t;tab];
else
    t = tab;
end
end

function t = runAndLoad(layout,t)
curdir = pwd;
cd('C:\sync\biomes\cellulose\optima\temp');
runComets(layout);
t = applyResult(layout,t);
cd(curdir);
end

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
end
