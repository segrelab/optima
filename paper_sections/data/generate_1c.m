%function dat3 = generate_1()
%Create the data for "1_prove_system_works.mlx"

%Run a simple 1x1 system for 100 min, cellulose, showing that you can vary these
%settings:
%  kcat, KM

p.version = 1;
v = getOptVars();

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0;
v.initcellulose = 10;
%v.initglc = 1e-1; %needed to bootstrap production
v.initglc = 0;
v.initmedia = 0; %amount of all other media added
v.initialpop = 0;
%v.initO2 = 1000; %oxygen level in the media
%v.initNH4 = 1000; %ammonium level in the media
%v.initPi = 1000;
%v.initSO4 = 1000;
%v.initK = 1000;
%v.initCO2 = 1000;
%v.initH2O = 1000;
v.deathrate = 0;
v.enzdecayperhour = 0;
v.costfactor = 3;
v.initenzyme = .1;

%%COMETS EXECUTION VARIABLES
%v.timestep = 1/3600; %1/(60*60); %1/3600 hours = 1 second.
v.timestep = 1/3600;
v.maxcycles = 3600; %2000/v.timestep;
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 10000;
v.spaceWidth = 10^(1/3);
v.lograte = 1;

v.kcat_cel = .1;
v.km_cel = 0.5;

v_default = v;

dat3 = table();

%%DEFAULT MODEL
model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904.mat'); %Saccaromyces cerevisiae iMM904
model_default = model_default.model;

%%Default run
layout_default = builddat3Layout(model_default,v);
dat3 = runAndLoad(layout_default,dat3,v);

%% vary kcat
v = v_default;
v.kcat_cel = .5;
layout = builddat3Layout(model_default,v);
dat3 = runAndLoad(layout,dat3,v);

v = v_default;
v.kcat_cel = .25;
layout = builddat3Layout(model_default,v);
dat3 = runAndLoad(layout,dat3,v);


%% vary KM
v = v_default;
v.km_cel = 0.75;
layout = builddat3Layout(model_default,v);
dat3 = runAndLoad(layout,dat3,v);

v = v_default;
v.km_cel = 0.25;
layout = builddat3Layout(model_default,v);
dat3 = runAndLoad(layout,dat3,v);


%% functions
function layout = builddat3Layout(model, v)
% enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
% enzmodel.v = v;
% layout = buildCelOptLayout(enzmodel);
% layout = setInitialPop(layout,'1x1',0);
% layout.models = {};
% layout = setMedia(layout,layout.mets,0);
layout = createLayout(model);
layout = setInitialPop(layout,'1x1');
layout = addExternalReaction(layout,'degrade_cellulose',{'cellulose','glc__D_e'},[-1 1],'enzyme', 'enzyme_e', 'k', v.kcat_cel, 'km', v.km_cel);
layout.params.writeMediaLog = true;
layout.params.timeStep = v.timestep;
layout.params.maxCycles = v.maxcycles;
layout = setMedia(layout,'enzyme_e',v.initenzyme);
layout = setMedia(layout,'cellulose',v.initcellulose);
end

function res = loadResult(layout,v)
medianames = {'glc__D_e' 'enzyme_e' 'cellulose'};
media = parseMediaLog(layout.params.mediaLogName,medianames);

res.timestep = unique(media.t);
res.t = res.timestep * v.timestep; %time in hours

res.glc_amt = media.amt(stridx('glc__D_e',media.metname));
res.enzyme_amt = media.amt(stridx('enzyme_e',media.metname));
res.cellulose_amt = media.amt(stridx('cellulose',media.metname));

res.glc_delta = calcDelta(res.glc_amt);
res.enzyme_delta = calcDelta(res.enzyme_amt);
res.cellulose_delta = calcDelta(res.cellulose_amt);


res.cellulolysis_rate = (3600 * v.timestep) * res.cellulose_delta;

end

function t = applyResult(layout, t, v)
%add the layout and its results to the table
res = loadResult(layout,v);
tab = table;

tab.layout = {layout};
tab.t = {res.t};
tab.timestep = {res.timestep};
tab.glc_amt = {res.glc_amt};
tab.enzyme_amt = {res.enzyme_amt};
tab.cellulose_amt = {res.cellulose_amt};
tab.glc_delta = {res.glc_delta};
tab.enzyme_delta = {res.enzyme_delta};
tab.cellulose_delta = {res.cellulose_delta};
tab.cellulolysis_rate = {res.cellulolysis_rate};

tab.alpha = v.alpha;
tab.decayrate = v.enzdecayperhour;
tab.costfactor = v.costfactor;
tab.km = v.km_cel;
tab.kcat = v.kcat_cel;

if size(t,1) > 0
    t = [t;tab];
else
    t = tab;
end
end

function t = runAndLoad(layout,t,v)
curdir = pwd;
cd('C:\sync\biomes\cellulose\optima\temp');
runComets(layout);
t = applyResult(layout,t,v);
cd(curdir);
end

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
end
