function layouts = part3_generation_remote(arg)

%Find the alpha that produces the maximum biomass yield across a range of
%initial biomass concentrations

%script derived from generate_2.m

p.version = 2;
%v1: 10/22/18: initial version derived from generate_2.m v4
%v2: 11/28/18: change model to YeastGEM v8.3

alphas = -(5:.25:13);
alphas = -(2:.25:4.75);
%alphas = -(9.25:.25:10);
%alphas = -(6:.1:9);
%alphas = -(1:12)
alphas = -(.25:.25:10);
alphas = -(.02:.02:5);
alphas = [-0.02 -(.1:.1:5)];
% alphas = -2.75;
% alphas = -[2.7:.01:2.8];
%celranges = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1 2 3 4 5 6 7 8 9 10] / 20;
%celranges = [.1 .2 .4 .6 .8 1 2 4 6 8 10] / 20;
celranges = [1/8 1/4 1/2 1 2 4 8];
%celranges = [1/8 8];
%celranges = [1/4 4];
%celranges = [1/2 1 2];
% celranges = [1 .9 .99 1.1 1.25];
% celranges = [.95:.025:1.05];
%celranges = 8;

v = arg.v;

variableInitEnzyme = arg.variableInitEnzyme;
v.enzInMaint = arg.enzInMaint;

% v = getOptVars();
% 
% %%ENZYME VARIABLES: taken from result of part 2
% v.km_cel = .03;
% v.kcat_cel = 27;
% 
% %%WORLD/LAYOUT VARIABLES
% 
% %INITIAL CELLULOSE
% %2g Avicel, 450mL stock, 100mL used
% mmol_glc = 1000 * 2 / (180.156 - 17); 
% mmol_avicel = mmol_glc / v.glcpercellulose;
% v.initcellulose = mmol_avicel * 100/450;
% 
% 
% % 
% % v.km_cel = 0.007;
% % v.kcat_cel = 3.2;
% 
% v.xdim = 1;
% v.ydim = 1;
% v.diff_cel = 0;
% v.initglc = 0;
% v.initmedia = 100; %amount of all other media added in rich media
% v.initO2 = 100; %oxygen level in the media
% v.initNH4 = 100; %ammonium level in the media
% v.initPi = 100;
% v.initSO4 = 100;
% v.initK = 100;
% v.initCO2 = 100;
% v.initH2O = 100;
% v.initergst = 0;%0.3967; %396.65 g/mol, .01g/L, 100mL
% v.initzymst = v.initergst;
% v.initgthox = 0;
% v.initfe2 = 100;
% v.deathrate = 0;
% v.enzdecayperhour = 0.01;
% v.costfactor = 1;
% 
% %GLC UPTAKE PARAMETERS
% v.km_glc = v.km_glc / 10;
% v.vmax_glc = v.vmax_glc * 1;
% v.glcuptakerate = v.vmax_glc;
% 
% %INITIAL POPULATION
% %"Cells were cultured to approximately 2x10^5 cells/mL
% % %Volume is 100mL, so 2x10^7 cells at 5e-11 g/cell = 1mg
% % v.initialpop = 1e-3;
% %max observed cell density was .27 g/L GDW. Use .027/100mL
% %max observed cell count is reported as 4e9
% %.027g/100mL / 4e9cells/ML
% weightpercell = .027/(4e9); %6.75e-12
% v.initialpop = 2e7 * weightpercell;
% 
% 
% 
% 
% v.alpha = 1e-1;
%v.initenzyme = v.initialpop * v.alpha;

%%COMETS EXECUTION VARIABLES
%v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
%v.maxcycles = 1000;  % 500/v.timestep; 
%v.maxcycles = 60; %2000/v.timestep;
%v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
%v.maxspacebiomass = 9999; %
%v.spaceWidth = 100^(1/3);
%v.lograte = 1;

v_default = v;

dat3 = table();
layouts = cell(0);

%%DEFAULT MODEL
%model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904
%model_default = load('C:\sync\biomes\models\yeast-GEM\ModelFiles\mat\yeastGEM_renamed.mat');
% model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = load('C:\sync\biomes\models\yeast-GEM-8.3.5-dev\ModelFiles\mat\yeastGEM_8-3-5-dev');
if isfield(model_default,'model')
    model_default = model_default.model;
end



v = v_default;



%% vary alpha and initcel
for i = 1:length(alphas)
    for j = celranges
        v = v_default;
        v.alpha = 10^(alphas(i));
        %v.initenzyme = v.initialpop * v.alpha;
        v.initcellulose = j * v.initcellulose;
        v.cellulosescale = j;
        if variableInitEnzyme
            v.initenzyme = v.initialpop * v.alpha;
        end
        layout = builddat3Layout(model_default,v);
        %dat3 = runAndLoad(layout,dat3);
        layouts = [layouts {layout}];
    end
end
end
%% functions
function layout = builddat3Layout(model, v)
% enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
% enzmodel.v = v;

enzmodel = prepareYeastGEMModel('model',model,'v',v);
enzmodel.v = v;


layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;

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
res.biomass_max = max(biomass.biomass);

res.initcel = layout.models{1}.v.initcellulose;
res.celscale = layout.models{1}.v.cellulosescale;
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
tab.kcat_cel = layout.models{1}.v.kcat_cel;
tab.km_cel = layout.models{1}.v.km_cel;
tab.vmax_glc = layout.models{1}.v.vmax_glc;
tab.km_glc = layout.models{1}.v.km_glc;

tab.initcel = res.initcel;
tab.biomass_max = res.biomass_max;
tab.celscale = res.celscale;

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
