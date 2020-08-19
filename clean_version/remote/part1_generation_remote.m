function runlayouts = part1_generation_remote(arg,remcfg)
%Create the data for "1_prove_system_works.mlx"

%Run a simple 1x1 system for 100 min, without cellulose, showing that you can vary these
%settings:
%  alpha
%  costfactor
%  decayrate

%Version history: 3: change model to yeastGEM.8.3.4.devel(Nov 11, 2019)

p.version = 3;
v = arg.v;

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0;
v.initcellulose = 0;
%v.initglc = 1e-1; %needed to bootstrap production
v.initglc = 1;
v.initmedia = 1; %amount of all other media added
v.initialpop = .0001;
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.deathrate = 0;
v.enzdecayperhour = 0;
v.costfactor = 1;

%%COMETS EXECUTION VARIABLES
%v.timestep = 1/3600; %1/(60*60); %1/3600 hours = 1 second.
v.timestep = 3;
v.maxcycles = 500; %2000/v.timestep;
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 1000;
v.spaceWidth = 10^(1/3);
v.lograte = 1;

v_default = v;

dat1 = table();
runlayouts = cell(0);

%%DEFAULT MODEL
model_default = arg.model;

%%Default run
layout_default = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout_default,dat1);
runlayouts = [runlayouts {layout_default}];

%% vary alpha
v = v_default;
v.alpha = v.alpha * 10;
layout = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout,dat1);
runlayouts = [runlayouts {layout}];

v = v_default;
v.alpha = v.alpha / 10;
layout = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout,dat1);
runlayouts = [runlayouts {layout}];

%% vary decay rate
v = v_default;
v.enzdecayperhour = .10;
layout = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout,dat1);
runlayouts = [runlayouts {layout}];

v = v_default;
v.enzdecayperhour = .0100;
layout = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout,dat1);
runlayouts = [runlayouts {layout}];

%% vary cost factor
v = v_default;
v.costfactor = 3;
layout = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout,dat1);
runlayouts = [runlayouts {layout}];

v = v_default;
v.costfactor = 9;
layout = buildDat1Layout(model_default,v);
%dat1 = runAndLoad(layout,dat1);
runlayouts = [runlayouts {layout}];

%% vary kcat


%% vary KM

%% Get the remote execution files


end

%% functions
function layout = buildDat1Layout(model, v)
%enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
%enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'costfactor',v.costfactor,'maintRxnAsObjective',v.maintRxnAsObjective,'atp',v.atp_per_peptide,'gtp',v.gtp_per_peptide);
enzmodel = prepareYeastGEMModel('model',model,'v',v);
if ~isfield(v,'variableInitEnzyme')
    v.variableInitEnzyme = false;
end
if ~isfield(v,'addEnzymeToMaint')
    v.addEnzymeToMaint = false;
end
enzmodel.v = v;
layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);


% 
% layout = setMedia(layout,'cellulose',v.initcellulose);
% layout = setMedia(layout,'enzyme[e]',v.initenzyme);
% layout = setMedia(layout,'zymst[e]',v.initzymst);
% layout = setMedia(layout,'ergst[e]',v.initergst);
% layout = setMedia(layout,'o2[e]',v.initO2);
% layout = setMedia(layout,'h[e]',v.initH);
% layout = setMedia(layout,'h2o[e]',v.initH2O);
% layout = setMedia(layout,'co2[e]',v.initCO2);
% layout = setMedia(layout,'pi[e]',v.initPi);
% layout = setMedia(layout,'so4[e]',v.initSO4);
% layout = setMedia(layout,'nh4[e]',v.initNH4);
% layout = setMedia(layout,'k[e]',v.initK);
% layout = setMedia(layout,'fe2[e]',v.initfe2);
% layout = setMedia(layout,'oleate[e]',v.initoleate);
% layout = setMedia(layout,'palmitoleate[e]',v.initpalmitoleate);
% 
% layout = setMedia(layout,'glc-D[e]',v.initglc);
% 
% layout = setMedia(layout,'gthox[e]',v.initgthox);
end

function res = loadResult(layout)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose'};% 'etoh[e]'};
media = parseMediaLog(layout.params.mediaLogName,medianames);

res.timestep = unique(media.t);
res.t = res.timestep * layout.models{1}.v.timestep; %time in hours

res.glc_amt = media.amt(stridx('glc-D[e]',media.metname));
res.enzyme_amt = media.amt(stridx('enzyme[e]',media.metname));
res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
%res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));

res.glc_delta = calcDelta(res.glc_amt);
res.enzyme_delta = calcDelta(res.enzyme_amt);
res.cellulose_delta = calcDelta(res.cellulose_amt);
%res.etoh_delta = calcDelta(res.etoh_amt);

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
%tab.etoh_amt = {res.etoh_amt};
tab.biomass = {res.biomass};
tab.glc_delta = {res.glc_delta};
tab.enzyme_delta = {res.enzyme_delta};
tab.cellulose_delta = {res.cellulose_delta};
%tab.etoh_delta = {res.etoh_delta};
tab.biomass_delta = {res.biomass_delta};

tab.alpha = layout.models{1}.v.alpha;
tab.decayrate = layout.models{1}.v.enzdecayperhour;
tab.costfactor = layout.models{1}.v.costfactor;

tab.variableInitEnzyme = layout.models{1}.v.variableInitEnzyme;
tab.addEnzymeToMaint = layout.models{1}.v.addEnzymeToMaint;
tab.maintRxnAsObjective = layout.models{1}.v.maintRxnAsObjective;

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
