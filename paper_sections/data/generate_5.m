%Create 2D layouts for observation of depletion zones around a colony
p.version = 2;
%v2: change model to YeastGEM 8.3
v = getOptVars();


v.km_cel = .03;
v.kcat_cel = 27;


%%WORLD/LAYOUT VARIABLES

initcellulose_gperliter = 20;

v.xdim = 21;
v.ydim = v.xdim;
v.diff_cel = 0;
v.initglc = 0;
v.initmedia = 10; %amount of all other media added in rich media
v.initialpop = 1e-2;
v.initO2 = 10; %oxygen level in the media
v.initNH4 = 10; %ammonium level in the media
v.initPi = 10;
v.initSO4 = 10;
v.initK = 10;
v.initCO2 = 10;
v.initH2O = 10;
v.initzymst = 10;
v.initergst = 10;
%v.initgthox = 0;
v.initfe2 = 100;
v.deathrate = 0;
v.enzdecayperhour = 0;
v.costfactor = 1;

v.glc_per_cellulose = 218.2; %if 1, cellobiose instead of cellulose in media

alphas = [exp(-2) exp(-3) exp(-4) exp(-5)];%exp(-4);

v.initenzyme = v.initialpop * v.alpha;

%%COMETS EXECUTION VARIABLES
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
v.maxcycles = 750/v.timestep; 
%v.maxcycles = 60; %2000/v.timestep;
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 50;
%v.maxspacebiomass = 1000; 
v.spaceWidth = .1^(1/3); %units are cm


initcellulose_mmolperliter = 1000 * initcellulose_gperliter / (180.156 * v.glcpercellulose);

spaceWidth = v.spaceWidth; %
spaceDepth = v.spaceWidth; %just in case I want to change this to a different constant
spaceVolume = v.spaceWidth * v.spaceWidth * spaceDepth; %in cm^3
spaceVolume_L = spaceVolume / 1000;
v.initcellulose = initcellulose_mmolperliter * spaceVolume_L;


v_default = v;
dat5 = table();


%Vary the diffusion rates
dde = v_default.diff_enz;
ddg = v_default.diff_glc;

%enzdiffs = [dde (dde * 10) (dde * 100) (dde / 10) (dde / 100)];
%glcdiffs = [ddg (ddg * 10) (ddg * 100) (ddg / 10) (ddg / 100)];
enzdiffs = dde;
glcdiffs = ddg;

%%DEFAULT MODEL
%model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904.mat'); %Saccaromyces cerevisiae iMM904
model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = model_default.model;

%%set up and launch runs
for k = 1:length(alphas)
    v.alpha = alphas(k);
    for i = 1:length(enzdiffs)
        de = enzdiffs(i);
        for j = 1:length(glcdiffs)
            dg = glcdiffs(j);
            v.diff_enzyme = de;
            v.diff_glc = dg;
            layout = buildDat5Layout(model_default,v);
            dat5 = runAndLoad(layout,dat5);
        end
        v = v_default;
    end
end

%% functions
function layout = buildDat5Layout(model, v)
%enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha);
enzmodel.v = v;
mid = ceil(v.xdim/2);
layout = buildCelOptLayout(enzmodel);
layout = setDims(layout,v.xdim,v.ydim);
layout = setInitialPop(layout,'colonies',[v.initialpop],false);
layout = setInitialMediaInCell(layout,mid,mid,'enzyme[e]',v.initenzyme);
layout = setStaticMedia(layout,1:v.xdim,1:v.ydim,'o2[e]',v.initO2);
%layout = setMedia(layout,'gthox_e',v.initgthox);
%layout = setMedia(layout,'zymst_e',v.initzymst);
%layout = setMedia(layout,'ergst_e',v.initergst);
%layout = setMedia(layout,'zymstest_SC_e',v.initzymst);
%layout = setMedia(layout,'ergstest_SC_e',v.initergst);
end

function res = loadResult(layout)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]'};
media = parseMediaLog(layout.params.mediaLogName,medianames);

media.time = media.t * layout.models{1}.v.timestep; %time in hours

res.media = media; 

%res.glc_amt = media.amt(stridx('glc__D_e',media.metname));
%res.enzyme_amt = media.amt(stridx('enzyme_e',media.metname));
%res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
%res.etoh_amt = media.amt(stridx('etoh_e',media.metname));

%res.glc_delta = calcDelta(res.glc_amt);
%res.enzyme_delta = calcDelta(res.enzyme_amt);
%res.cellulose_delta = calcDelta(res.cellulose_amt);
%res.etoh_delta = calcDelta(res.etoh_amt);

biomass = parseBiomassLog(layout.params.biomassLogName);
res.biomass= biomass;
%res.biomass_delta = calcDelta(res.biomass);

res.diff_enz = layout.models{1}.v.diff_enz;
res.diff_glc = layout.models{1}.v.diff_glc;

end

function t = applyResult(layout, t)
%add the layout and its results to the table
res = loadResult(layout);
tab = table;

tab.layout = {layout};
tab.media = {res.media};
tab.biomass = {res.biomass};

% 
% tab.t = {res.t};
% tab.timestep = {res.timestep};
% tab.glc_amt = {res.glc_amt};
% tab.enzyme_amt = {res.enzyme_amt};
% tab.cellulose_amt = {res.cellulose_amt};
% tab.etoh_amt = {res.etoh_amt};
% tab.biomass = {res.biomass};
% tab.glc_delta = {res.glc_delta};
% tab.enzyme_delta = {res.enzyme_delta};
% tab.cellulose_delta = {res.cellulose_delta};
% tab.etoh_delta = {res.etoh_delta};
% tab.biomass_delta = {res.biomass_delta};

tab.alpha = layout.models{1}.v.alpha;
tab.log_alpha = log(layout.models{1}.v.alpha);
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
