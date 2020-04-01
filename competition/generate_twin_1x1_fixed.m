function dat = generate_twin_1x1_fixed(dat)
%GENERATE_CHEATER_1X1_FIXED Create a data table which includes two identical copies of a model

%Version History:
%   1: Initial Form. Copied from generate_cheater_1x1_fixed V5

if exist('v','var')
    clear(v);
end
if exist('p','var')
    clear(p)
end

p.version = 1;

p.repititions = 1; %How many times should an experiment be repeated with a different random seed?
p.timestep = [.05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .7 .8 .9 1 1.1 1.25];
p.timestep = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.1 1.25];

v = getOptVars();

%p.deathrate = [0 .1 .2 .05 .01];
%p.deathrate = [0 0.01];
p.deathrate = 0;
%p.deathrate = [.1 .2];
p.dilutionrate = 0; %fraction of media replaced per hour
v.richmedia = false; %provide a bunch of everything

p.maintenance = true; %keep the ATP maintenance flux?
p.useMaintenanceHack = [true]; %add the reaction which allows partial maintenance
v.useMaintenanceHack = p.useMaintenanceHack(1);
p.deathrate = 0;%1/720;

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0;
v.initcellulose = 1;
v.initglc = 10;
v.initmedia = 1; %amount of all other media added
v.initialpop = 1e-6;
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.deathrate = p.deathrate;

%%COMETS EXECUTION VARIABLES
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
v.maxcycles = floor(100 / v.timestep); %1500h ~ 2 months
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
model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904
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
            if ncols == 0 || ~any(dat.alpha == p.alpha(i) & dat.decayrate == p.enzdecayperhour(j) & dat.costfactor == p.costfactor(k))
                v.useMaintenanceHack = p.useMaintenanceHack(maint);
                v.timestep = p.timestep(f);
                v.maxcycles = floor(p.maxcycles / v.timestep);
                
                layout = buildTwinLayout(model_default,v);
                
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
                dat.biomass1{n} = res.biomass1;
                dat.biomass2{n} = res.biomass2;
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
dat.maint_replaced = cellfun(@(x) logical(x),dat.maint_replaced);
cd(curdir);
end
function layout = buildTwinLayout(model, v)
if v.useMaintenanceHack %apply the Mainteance Subreaction method to allow partial satisfaction of the ATP demand
    model = addMaintenanceSubreaction(model);
end
model2 = model;
model2.description = 'Copy 2';
model.v = v;
layout = buildCelOptLayout(model);
layout = addModel(layout,model2);
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
res.biomass1 = biomass.biomass(biomass.model == 0);
res.biomass2 = biomass.biomass(biomass.model == 1);

end

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
end

function layout = cleanupLayout(layout)
layout.models = {};
end