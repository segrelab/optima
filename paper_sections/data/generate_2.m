%function dat1 = generate_1()
%Create the data for "2_fit_in_vivo_data.mlx"

p.version = 5;
%v2: 7/12/2018: add ergosterol and zymosterol exchange rxns to buildSCerevisiaeModel so that anaerobic sims can be performed
%v3: 8/1/2018: Change amino acide recipe for enzyme to reflect average
%costs across 3 species
%v4: 8/2: remove static media components from buildceloptlayout
%v5: 11/18: change model to yeastGEM8.3.2

alphas = -(1:8);
% alphas = -(1:.5:8);
% alphas = -5.7;
% alphas = [-3 -6];
%alphas = -2.6;
%alphas = -3.3;
alphas = -3.2;

varykineticparams = true;
varyalphas = false;
varyallparams = false;

v = getOptVars();

%%ENZYME VARIABLES

% turn this up to see what a powerful response looks like
% v.kcat_cel = 10000 * v.kcat_cel;
% v.km_cel = v.km_cel / 10000;


%%WORLD/LAYOUT VARIABLES

v.initcellulose = .5; %carbon source for this experiment is 5mM p-nitrophenyl-beta-D-glucopyranoside, in 100mL culture
%v.glcpercellulose = 1;

v.km_cel = .03;
v.kcat_cel = 27;
% 
% v.km_cel = 0.007;
% v.kcat_cel = 3.2;

v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0;
v.initmedia = 0; %amount of all other media added in rich media
v.initO2 = 0; %oxygen level in the media
v.initNH4 = 100; %ammonium level in the media
v.initPi = 100;
v.initSO4 = 100;
v.initK = 100;
v.initCO2 = 100;
v.initH2O = 100;
v.initergst = 0.3967; %396.65 g/mol, .01g/L, 100mL
v.initzymst = 0; %v.initergst;
v.initglc = 0;
v.initgthox = 10;
v.initfe2 = 100;
v.deathrate = 0;
v.enzdecayperhour = 0.01;
v.costfactor = 1;

%GLC UPTAKE PARAMETERS
v.km_glc = v.km_glc / 10;
v.vmax_glc = v.vmax_glc * 1;
v.glcuptakerate = v.vmax_glc;

%INITIAL POPULATION
%"Cells were cultured to approximately 2x10^5 cells/mL
% %Volume is 100mL, so 2x10^7 cells at 5e-11 g/cell = 1mg
% v.initialpop = 1e-3;
%max observed cell density was .27 g/L GDW. Use .027/100mL
%max observed cell count is reported as 4e9
%.027g/100mL / 4e9cells/ML
weightpercell = .027/(4e9); %6.75e-12
v.initialpop = 2e7 * weightpercell;


%INITIAL CELLULOSE
%2g Avicel, 450mL stock, 100mL used
mmol_glc = 1000 * 2 / (180.156 - 17); 
mmol_avicel = mmol_glc / v.glcpercellulose;
v.initcellulose = mmol_avicel * 100/450;


v.alpha = 1e-1;
%v.initenzyme = v.initialpop * v.alpha;
v.initenzyme = v.initialpop * exp(-8);


%%COMETS EXECUTION VARIABLES
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
v.maxcycles = 250;  % 500/v.timestep; 
%v.maxcycles = 60; %2000/v.timestep;
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 9999; %
v.spaceWidth = 100^(1/3);
v.lograte = 1;

v_default = v;

dat2 = table();

%%DEFAULT MODEL
%model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904
%model_default = load('C:\sync\biomes\models\yeast-GEM\ModelFiles\mat\yeastGEM_renamed.mat');
%model_default = load('C:\sync\biomes\models\Scerevisiae_iIN800_renamed.mat');
model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = model_default.model;


%% Initial Enzyme Concentration
%In order to satisfy the ATP maintenance flux and still have some glucose
%left over to produce biomass and enzyme, the system must be in a state
%where delta(glc) is positive. That is, v_catalysis > v_consumption

%v_consumption = lb_glc_exchange =biomass * vmax_glc_uptake * glc / (km_glc_uptake + glc)
%v_catalysis = [enzyme] * kcat * [cellulose] / (km + [cellulose])
%therefore [enzyme] > ( lb_glc_exchange ) / ( kcat * glcpercellulose *
%[cellulose] / (km + [cellulose]) )

% % bound_glc = getMinGlucoseFlux(model_default);
% % bound_glc = v.timestep * bound_glc;
% % overshot_factor = 1.01; %how much larger is the initial glc delta than necessary?
% % kcat_per_timestep = 3600 * v.kcat_cel * v.timestep; %kcat is in seconds, 
% % v.initenzyme = overshot_factor * bound_glc * v.initialpop/ ...
% %     (kcat_per_timestep * v.initcellulose / ...
% %     (v.km_cel + v.initcellulose));

%let's give it some glucose to start off as well, to skip the initial lag
%v.initglc = bound_glc * v.initialpop;
%disable initial glc
v.initglc = 0;

v.initcellulose = v.initcellulose - (v.initglc/v.glcpercellulose);


%% Default run
%layout_default = buildDat2Layout(model_default,v);
%dat2 = runAndLoad(layout_default,dat2);

%% Trial 1, Den Haan, Rose, Lynd and van Zyl 2006
%experiment conditions:
%   -anaerobic
%   -substrate: amorphous phosphoric acid swollen cellulose (PASC)
%   -substrate recipe: 100ml YP-PASC medium supplemented with 0.01gl egosterol and 0.42g/l Tween 80
%   -cultures were inoculated to approcimately 2x10^5 cells/ml 

%organism details:
%   -origin strain: S.cerevisiae Y294 (ATCC 201160)
%   -endoglucanase cel7B from Trichoderma reesei
%   -beta-glucosidase BGL1 from Saccharomycopsis fubuligera
%   -NO exoglucanase
%   -plasmid expression under control of the yENO1 promoter

%key observations:
%   -one unit of enzyme activity is defined as the amount of enzyme required to produce 1 umol of reducing sugar per minute under the assay conditions
%   -maximum endoglucanase activity in Y294[CEL5] was 0.30 U/mg*DCW
%   -maximum beta-glucosidase activity in Y294[CEL5] was 0.48 U/mg*DCW
%   -maximum cell density observed was 3.88e7 cells/mL
%   -ethanol production observed up to 1.0 g/L

%assumptions
%   -due to lack of exoglucanase, assume the bottleneck in cellulose breakdown is endoglucanase activity
%   -endoglucanase and beta-glucosidase are produced at the same rate and have the same cost
%   -assuming a yeast cell weighs 5e-11 g, or 5e-8 mg 

% current issue/blocker: conversion from cell dry weight <-> number of cells
% for the initial approach, don't worry about mapping the data exactly, normalize the curves to have the same maximum


%v = v_default;
%cellweight = 5e-8; %mg
%v.initialpop = cellweight * 2e5;

%parameters for endoglucanase from T reesei
%v.kcat = .01; %source: BRENDA brenda-enzymes.org/enzyme.php?ecno=3.2.1.4#TURNOVER%20NUMBER%20[1/s]
%v.km = .125; %source: BRENDA brenda-enzymes.org/enzyme.php?ecno=3.2.1.4#KM%20VALUE%20[mM]

% layout = buildDat2Layout(model_default,v);
% %layout = setMedia(layout,'o2[e]',0); %anaerobic
% dat2 = runAndLoad(layout,dat2);


%% vary alpha
%alphas = -3:.1:-1;
v_default = v;
if (varyalphas)
    for i = 1:length(alphas)
        v = v_default;
        v.alpha = exp(alphas(i));
        %v.initenzyme = v.initialpop * v.alpha;
        %v.initenzyme = 1000000;
        %     v.initglc = 1000;
        layout = buildDat2Layout(model_default,v);
        dat2 = runAndLoad(layout,dat2);
    end
end
%% vary kinetic params
%ROUND 1: alpha = exp(-2.6)
% kcats = power(10,-(.1:.2:1.5));
% kms = power(10,-(0:.5:4));

%ROUND 2: alpha = exp(-2.6)
%kcats = [.1 .5 1 3 5 8 12 15 20 25 30 35 40 45 50];
kcats = 20;
kms = .005:.03:.515;

if (varykineticparams)
    dat2_kin = dat2;
    for a = alphas
        for i = kcats
            for j = kms
                v = v_default;
                v.alpha = exp(a);
                v.initenzyme = v.initialpop * v.alpha;
                v.kcat_cel = i;
                %v.kcat_cel = v.kcat_cel * i;
                v.km_cel = v.km_cel * j;
                %v.km_Cel = j;
                layout = buildDat2Layout(model_default,v);
                dat2_kin = runAndLoad(layout,dat2_kin);
            end
        end
    end
end

% %ROUND 3
% %search-type kinetic variability
% lastkm = v.km_cel;
% lastkcat = v.kcat_cel;
% coolrate = .25;
% stepsize = 1;
% runningsearch = true;
% iterations = 3;
% laststep_km = 1; %was the last step an increase?
% laststep_kcat = 1;
% lastimproved_km = 1;
% lastimproved_kcat = 1;
% lastrmse = 99;
% if (varykineticparams)
%     dat2_kin = dat2;
%     replicate_in_vivo_data;
%     for a = alphas
%         while runningsearch
%             v = v_default;
%             v.alpha = exp(a);
%             v.initenzyme = v.initialpop * v.alpha;
% 
%             %run the varied KM attempt
%             v.kcat_cel = lastkcat;
%             kmstep = lastkm * stepsize * laststep_km * lastimproved_km;
%             v.km_cel = lastkm + kmstep;
%             layout = buildDat2Layout(model_default,v);
%             dat2_kin = runAndLoad(layout,dat2_kin);
%             %process the varied KM attempt
%             if v.km_cel > lastkm
%                 laststep_km = 1;
%             else
%                 laststep_km = -1;
%             end
%             [rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2_kin.t{end},dat2_kin.biomass{end},3,'log10');
%             if rmse <= lastrmse
%                 lastrmse = rmse;
%                 lastimproved_km = 1;
%                 lastkm = v.km_cel;
%             else
%                 lastimproved_km = -1;
%             end
%                 
%             %run the varied KCAT attempt
%             v.km_cel = lastkm;
%             kcatstep = lastkcat * stepsize * laststep_kcat * lastimproved_kcat;
%             v.kcat_cel = lastkcat + kcatstep;
%             layout = buildDat2Layout(model_default,v);
%             dat2_kin = runAndLoad(layout,dat2_kin);
%             %process the varied KCAT attempt
%             if v.kcat_cel > lastkcat
%                 laststep_kcat = 1;
%             else
%                 laststep_kcat = -1;
%             end
%             [rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2_kin.t{end},dat2_kin.biomass{end},3,'log10');
%             if rmse <= lastrmse
%                 lastrmse = rmse;
%                 lastimproved_kcat = 1;
%                 lastkcat = v.kcat_cel;
%             else
%                 lastimproved_kcat = -1;
%             end
%             
%             iterations = iterations - 1;
%             if iterations <= 0
%                 runningsearch = false;
%             end
%             
%             stepsize = stepsize * (1-coolrate);
%         end
%     end
% end


%search-type kinetic variability, also walking along alpha
lastkm = v.km_cel;
lastkcat = v.kcat_cel;
lastalpha = exp(alphas(1));
coolrate = .25;
stepsize = 1;
runningsearch = true;
iterations = 10;
laststep_km = 1; %was the last step an increase?
laststep_kcat = 1;
laststep_alpha = 1;
lastimproved_km = 1;
lastimproved_kcat = 1;
lastimproved_alpha = 1;
lastrmse = 99;
if (varyallparams)
    dat2_kin = dat2;
    replicate_in_vivo_data;
    while runningsearch
        v = v_default;
        
        %run the varied KM attempt
        v.alpha = lastalpha;
        %v.initenzyme = v.initialpop * v.alpha;
        v.kcat_cel = lastkcat;
        kmstep = lastkm * stepsize * laststep_km * lastimproved_km;
        v.km_cel = lastkm + kmstep;
        layout = buildDat2Layout(model_default,v);
        dat2_kin = runAndLoad(layout,dat2_kin);
        %process the varied KM attempt
        if v.km_cel > lastkm
            laststep_km = 1;
        else
            laststep_km = -1;
        end
        [rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2_kin.t{end},dat2_kin.biomass{end},3,'log10');
        if rmse <= lastrmse
            lastrmse = rmse;
            lastimproved_km = 1;
            lastkm = v.km_cel;
        else
            lastimproved_km = -1;
        end
        
        %run the varied KCAT attempt
        v.alpha = lastalpha;
        %v.initenzyme = v.initialpop * v.alpha;
        v.km_cel = lastkm;
        kcatstep = lastkcat * stepsize * laststep_kcat * lastimproved_kcat;
        v.kcat_cel = lastkcat + kcatstep;
        layout = buildDat2Layout(model_default,v);
        dat2_kin = runAndLoad(layout,dat2_kin);
        %process the varied KCAT attempt
        if v.kcat_cel > lastkcat
            laststep_kcat = 1;
        else
            laststep_kcat = -1;
        end
        [rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2_kin.t{end},dat2_kin.biomass{end},3,'log10');
        if rmse <= lastrmse
            lastrmse = rmse;
            lastimproved_kcat = 1;
            lastkcat = v.kcat_cel;
        else
            lastimproved_kcat = -1;
        end
        
        %run the varied ALPHA attempt
        v.kcat_cel = lastkcat;
        v.km_cel = lastkm;
        alphastep = lastalpha * stepsize * laststep_alpha * lastimproved_alpha;
        v.alpha = lastalpha + alphastep;
        %v.initenzyme = v.initialpop * v.alpha;
        
        layout = buildDat2Layout(model_default,v);
        dat2_kin = runAndLoad(layout,dat2_kin);
        %process the varied KCAT attempt
        if v.alpha > lastalpha
            laststep_alpha = 1;
        else
            laststep_alpha = -1;
        end
        [rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2_kin.t{end},dat2_kin.biomass{end},3,'log10');
        if rmse <= lastrmse
            lastrmse = rmse;
            lastimproved_alpha = 1;
            lastalpha = v.alpha;
        else
            lastimproved_alpha = -1;
        end
        
        iterations = iterations - 1;
        if iterations <= 0
            runningsearch = false;
        end
        
        stepsize = stepsize * (1-coolrate);
    end
end



%% functions
function layout = buildDat2Layout(model, v)
%enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha);
enzmodel.v = v;

% enzmodel.lb(stridx('EX_zymst',enzmodel.rxns)) = -1e6;
% enzmodel.lb(stridx('EX_ergst',enzmodel.rxns)) = -1e6;
% enzmodel.lb(stridx('EX_co2',enzmodel.rxns)) = -1e6;

%enzmodel = addReaction(enzmodel,'EX_co2[e]',{'co2[e]'},-1,'lb',-10,'ub',10000);

layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);

layout = setMedia(layout,layout.mets,0);

layout = setMedia(layout,'cellulose',v.initcellulose);
layout = setMedia(layout,'enzyme[e]',v.initenzyme);
layout = setMedia(layout,'zymst[e]',v.initzymst);
layout = setMedia(layout,'ergst[e]',v.initergst);
% layout = setMedia(layout,'zymstest_SC[e]',v.initzymst);
% layout = setMedia(layout,'ergstest_SC[e]',v.initergst);
layout = setMedia(layout,'o2[e]',v.initO2);
layout = setMedia(layout,'h2o[e]',v.initH2O);
layout = setMedia(layout,'co2[e]',v.initCO2);
layout = setMedia(layout,'pi[e]',v.initPi);
layout = setMedia(layout,'so4[e]',v.initSO4);
layout = setMedia(layout,'nh4[e]',v.initNH4);
layout = setMedia(layout,'k[e]',v.initK);
layout = setMedia(layout,'fe2[e]',v.initfe2);

layout = setMedia(layout,'glc-D[e]',v.initglc);

layout = setMedia(layout,'gthox[e]',v.initgthox);

% layout.params.defaultVmax = 1000;
% layout.params.defaultKm = 1e-8;

% %%add some O2 dependent metabolites
% %dependentmets = {'pa_SC[c]' 'pc_SC[c]' 'pe_SC[c]' 'ps_SC[c]' 'triglyc_SC[c]'};
% dependentmets = {'triglyc_SC[c]'};
% %layout.models{1} = addExchangeRxn(layout.models{1},dependentmets);
% layout.models{1} = addReaction(layout.models{1},'EX_triglyc_SC[c]',{'triglyc_SC[c]'},-1,'lowerBound',-99999);
% dep_amt = 10;
% layout = setMedia(layout,dependentmets,dep_amt);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
%temporary
% layout = setMedia(layout,'glc-D[e]',5);
% layout = setMedia(layout,'o2[e]',5);
end

function res = loadResult(layout)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]' 'gthox[e]'};
media = parseMediaLog(layout.params.mediaLogName,medianames);

res.timestep = unique(media.t);
res.t = res.timestep * layout.models{1}.v.timestep; %time in hours

res.glc_amt = media.amt(stridx('glc-D[e]',media.metname));
res.enzyme_amt = media.amt(stridx('enzyme[e]',media.metname));
res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));
res.gthox_amt = media.amt(stridx('gthox[e]',media.metname));

res.glc_delta = calcDelta(res.glc_amt);
res.enzyme_delta = calcDelta(res.enzyme_amt);
res.cellulose_delta = calcDelta(res.cellulose_amt);
res.etoh_delta = calcDelta(res.etoh_amt);
res.gthox_delta = calcDelta(res.gthox_amt);

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
tab.gthox_amt = {res.gthox_amt};
tab.biomass = {res.biomass};
tab.glc_delta = {res.glc_delta};
tab.enzyme_delta = {res.enzyme_delta};
tab.cellulose_delta = {res.cellulose_delta};
tab.etoh_delta = {res.etoh_delta};
tab.gthox_delta = {res.gthox_delta};
tab.biomass_delta = {res.biomass_delta};

tab.alpha = layout.models{1}.v.alpha;
tab.decayrate = layout.models{1}.v.enzdecayperhour;
tab.costfactor = layout.models{1}.v.costfactor;
tab.kcat_cel = layout.models{1}.v.kcat_cel;
tab.km_cel = layout.models{1}.v.km_cel;
tab.vmax_glc = layout.models{1}.v.vmax_glc;
tab.km_glc = layout.models{1}.v.km_glc;

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
