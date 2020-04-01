function dat = generate_cheater_1x1_fixed(dat)
%GENERATE_CHEATER_1X1_FIXED Create a data table which includes an enzyme
%producer and a non-producer

%Version History:
%   1: Initial Form
%   2: randomOrder-> True
%   3: restore ATP maintenance reaction demand to Producer model
%   4: create addMaintenanceSubreaction function to allow 'nonviable'
%   models to take up glc and partially meet energy demand (6/20/18)
%   5: 9/6/2018 Remove addMaintenanceSubreaction. Change enzyme stoich to be based on average of enzymes from
%   BRENDA. Alter media names to use [e] format, switched to '_renamed'
%   model. Removed initial glucose, replaced with initial enzyme = alpha *
%   initialpop
%   6: Switch to YeastGEM v8.3 model, use kinetics for optimal fitting to
%   denHann experiment

if exist('v','var')
    clear(v);
end
if exist('p','var')
    clear(p);
end

p.version = 6;

p.repititions = 4; %How many times should an experiment be repeated with a different random seed?
%p.timestep = [.05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .7 .8 .9 1 1.1 1.25];
p.timestep = .5;

v = getOptVars();

%p.alpha = [.0001 .00005 .0005 .001 .005 .01 .015 .02 .03 .04 .06 .1 .13 .175 .2];
%p.alpha = [p.alpha exp(-6.8) exp(-6.6) exp(-6.4) exp(-6.2) exp(-6) exp(-5.8) exp(-5.6) exp(-5.4) exp(-5.2) exp(-5) exp(-8)];
%p.alpha = [exp(-7) exp(-6.6) exp(-6.3) exp(-6) exp(-5.6) exp(-5.3) exp(-5) exp(-4.5) exp(-4) exp(-3.5) exp(-3) exp(-2.5) exp(-2) exp(-1.5) exp(-1)];
alphas = [-1 -2 -3 -3.25 -4 -5 -6 -7 -8];
p.alpha = exp(alphas);
%p.alpha = exp(-5);
%p.alpha = [exp(-7) exp(-6) exp(-5)];
%p.alpha = [exp(-6) exp(-5)];

%p.enzdecayperhour = [0 1/2400 1/1200 1/600 1/4800 1/12000 1/240 1/120];
%p.enzdecayperhour = [0 1/1200 1/600 1/2400 1/240];
p.enzdecayperhour = .01;

%p.deathrate = [0 .1 .2 .05 .01];
%p.deathrate = [0 0.01];
p.deathrate = 0;
%p.deathrate = [.1 .2];
p.dilutionrate = 0; %fraction of media replaced per hour
%p.costfactor = [1 .5 2];
%p.costfactor = [1 3];
p.costfactor = 1;
v.richmedia = false; %provide a bunch of everything

p.maintenance = true; %keep the ATP maintenance flux?
p.useMaintenanceHack = [false]; %add the reaction which allows partial maintenance
v.useMaintenanceHack = p.useMaintenanceHack(1);
p.deathrate = 0;%1/720;

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0;
v.initcellulose = 0.5;
v.initglc = 0; 
v.initmedia = 1000; %amount of all other media added
v.initialpop = 1e-3;
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.initFe2 = 1000;
v.deathrate = p.deathrate;

%%COMETS EXECUTION VARIABLES
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
v.maxcycles = floor(800 / v.timestep); %1500h ~ 2 months
p.maxcycles = v.maxcycles;
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 10000;
v.spaceWidth = 10^(1/3);
v.lograte = 1;

%print values
p = p
v.p = p

%%DEFAULT MODEL
% model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904
model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = model_default.model;
%S cerevisiae model ref: http://bigg.ucsd.edu/models/iMM904 https://www.ncbi.nlm.nih.gov/pubmed/19321003

if ~p.maintenance %remove the ATP maintenance flux from the model
    model_default.lb(model_default.lb > 0) = 0;
end

if nargin < 1
    dat = table();
end
[nrows, ncols] = size(dat);
n = nrows + 1;

curdir = pwd;
tdir = 'C:\sync\biomes\cellulose\optima\temp';
cd(tdir);
for maint = 1:length(p.useMaintenanceHack)
for rep = 1:p.repititions
for f = 1:length(p.timestep)
for i = 1:length(p.alpha)
    for j = 1:length(p.enzdecayperhour)
        for k = 1:length(p.costfactor)
            if ncols == 0 || ~any(dat.alpha == p.alpha(i) & dat.decayrate == p.enzdecayperhour(j) & dat.costfactor == p.costfactor(k))
                v.useMaintenanceHack = p.useMaintenanceHack(maint);
                v.timestep = p.timestep(f);
                v.maxcycles = floor(p.maxcycles / v.timestep);
                
                v.alpha = p.alpha(i);
                v.enzdecayperhour = p.enzdecayperhour(j);
                v.costfactor = p.costfactor(k);
                layout = buildCheaterLayout(model_default,v);
                
                seed = floor(random('Uniform',0,1000));
                layout.params.randomSeed = seed;
                
                output = runCheaterLayout(layout);
                res = loadCheaterResult(layout);
                
                dat.alpha(n) = v.alpha;
                dat.decayrate(n) = v.enzdecayperhour;
                dat.costfactor(n) = v.costfactor;
                dat.timestep{n} = res.timestep;
                dat.t{n} = res.t;
                dat.glc_amt{n} = res.glc_amt;
                dat.glc_delta{n} = res.glc_delta;
                dat.enzyme_amt{n} = res.enzyme_amt;
                dat.enzyme_delta{n} = res.enzyme_delta;
                dat.cellulose_amt{n} = res.cellulose_amt;
                dat.cellulose_delta{n} = res.cellulose_delta;
                dat.etoh_amt{n} = res.etoh_amt;
                dat.etoh_delta{n} = res.etoh_delta;
                dat.biomass_producer{n} = res.biomass_producer;
                dat.biomass_cheater{n} = res.biomass_cheater;
                dat.layout{n} = cleanupLayout(layout);
                dat.v{n} = v;
                dat.sim_complete(n) = (res.cellulose_amt(end) + res.glc_amt(end)) <= 1e-9;
                dat.maint_replaced{n} = v.useMaintenanceHack;
                dat.ts(n) = v.timestep;
                dat.random_seed(n) = seed;
                n = n + 1;
            end
        end
    end
end
end
end
end
dat.maint_replaced = cellfun(@(x) logical(x),dat.maint_replaced);
cd(curdir);
end

function layout = buildCheaterLayout(model, v)
if v.useMaintenanceHack %apply the Mainteance Subreaction method to allow partial satisfaction of the ATP demand
    model = addMaintenanceSubreaction(model);
end
%enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha);
enzmodel.v = v;
layout = buildCelOptLayout(enzmodel);
cheatermodel = prepareYeastGEMCheater(model);
layout = addModel(layout,cheatermodel);
layout = setInitialPop(layout,'1x1',[v.initialpop,v.initialpop]);
layout = setMedia(layout,'enzyme[e]',v.alpha * v.initialpop);
layout.params.randomOrder = true;
end

function res = runCheaterLayout(layout)
res = runComets(layout);
end

function res = loadCheaterResult(layout)
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
res.biomass_producer = biomass.biomass(biomass.model == 0);
res.biomass_cheater = biomass.biomass(biomass.model == 1);

end

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
end

function layout = cleanupLayout(layout)
layout.models = {};
end

function model = prepareYeastGEMCheater(model)
model.modifications = cell(0);
model.S = sparse(model.S);
model.modifications = [model.modifications {'Made S matrix sparse'}];

lb = -10;
model.modifications = [model.modifications {['Change ''glutathione disulfide exchange'' lb from ' num2str(model.lb(stridx('glutathione disulfide exchange',model.rxnNames))) ' to ' lb]}];
model.lb(stridx('glutathione disulfide exchange',model.rxnNames)) = lb;

model.mets{268} = 'ac[e]';
model.mets{288} = 'adp[c]';
model.mets{311} = 'nh4[e]';
model.mets{324} = 'atp[c]';
model.mets{345} = 'co2[e]';
model.mets{431} = 'glc-D[e]';
model.mets{512} = 'ergst[e]';
model.mets{520} = 'etoh[e]';
model.mets{556} = 'gdp[c]';
model.mets{572} = 'gthox[e]';
model.mets{592} = 'gtp[c]';
model.mets{612} = 'h2o[e]';
model.mets{720} = 'fe2[e]';
model.mets{748} = 'ala-L[c]';
model.mets{758} = 'arg-L[c]';
model.mets{762} = 'asn-L[c]';
model.mets{766} = 'asp-L[c]';
model.mets{774} = 'cys-L[c]';
model.mets{784} = 'glu-L[c]';
model.mets{791} = 'gln-L[c]';
model.mets{795} = 'gly[c]';
model.mets{798} = 'his-L[c]';
model.mets{807} = 'ile-L[c]';
model.mets{812} = 'leu-L[c]';
model.mets{816} = 'lys-L[c]';
model.mets{820} = 'met-L[c]';
model.mets{823} = 'phe-L[c]';
model.mets{826} = 'pro-L[c]';
model.mets{830} = 'ser-L[c]';
model.mets{836} = 'thr-L[c]';
model.mets{839} = 'trp-L[c]';
model.mets{842} = 'tyr-L[c]';
model.mets{847} = 'val-L[c]';
model.mets{1005} = 'o2[e]';
model.mets{1037} = 'pi[e]';
model.mets{1055} = 'k[e]';
model.mets{1135} = 'so4[e]';
model.mets{1215} = 'zymst[e]';

model.modifications = [model.modifications {'Replaced the following ids in mets field: ac[e] adp[c] nh4[e] atp[c] co2[e] glc-D[e] ergst[e] etoh[e] gdp[c] gthox[e] gtp[c] h2o[e] fe2[e] ala-L[c] arg-L[c] asn-L[c] asp-L[c] cys-L[c] glu-L[c] gln-L[c] gly[c] his-L[c] ile-L[c] leu-L[c] lys-L[c] met-L[c] phe-L[c] pro-L[c] ser-L[c] thr-L[c] trp-L[c] tyr-L[c] val-L[c] o2[e] pi[e] k[e] so4[e] zymst[e]'}];

end