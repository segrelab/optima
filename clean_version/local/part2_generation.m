function dat2 = part2_generation(args)
%PAET2_GENERATION Summary of this function goes here
%   Detailed explanation goes here



v = args.v;
v_default = v;
iterations = args.iterations;
dat2 = table();

model_default = args.model;

%search-type kinetic variability, also walking along alpha
lastkm = v.km_cel;
lastkcat = v.kcat_cel;
lastalpha = v.alpha;
coolrate = .25;
stepsize = 1;
if isfield(args,'stepsize')
    stepsize = args.stepsize;
end
runningsearch = iterations > 0;
laststep_km = 1; %was the last step an increase?
laststep_kcat = 1;
laststep_alpha = 1;
lastimproved_km = 1;
lastimproved_kcat = 1;
lastimproved_alpha = -1;
roundalpha = false;
if isfield(args,'roundalpha') %round logalpha to the nearest 0.1?
    roundalpha = args.roundalpha;
end
roundkm = true;
if isfield(args,'roundkm')
    roundkm = args.roundkm;
end

permuteKM = true;
permuteKcat = true;
permuteAlpha = true;
variableInitEnzyme = false;
if isfield(args,'permuteKM')
    permuteKM = args.permuteKM;
end
if isfield(args,'permuteKcat')
    permuteKcat = args.permuteKcat;
end
if isfield(args,'permuteAlpha')
    permuteAlpha = args.permuteAlpha;
end
if isfield(args,'variableInitEnzyme')
    variableInitEnzyme = args.variableInitEnzyme;
end

sigfigs_alpha = 1;
sigfigs_km = 2;

replicate_in_vivo_data;

layout = buildDat2Layout(model_default, v);
dat2 = runAndLoad(layout,dat2,ivd);
lastrmse = dat2.rmse(1);

while runningsearch
    v = v_default;
    
    if permuteKM
        %run the varied KM attempt
        v.alpha = lastalpha;
        v.kcat_cel = lastkcat;
        kmstep = 2^(stepsize * laststep_km * lastimproved_km);
        nextkm = lastkm * kmstep;

        if roundkm
            minstep_km = 10^(-sigfigs_km);
            if nextkm > lastkm
                nextkm = max(round(nextkm,sigfigs_km),lastkm + minstep_km);
            else
                nextkm = min(round(nextkm,sigfigs_km),lastkm - minstep_km);
            end
            if nextkm < minstep_km
                nextkm = minstep_km;
            end
        end
        v.km_cel = nextkm;
        if isNovel(v,dat2)
            layout = buildDat2Layout(model_default,v);
            dat2 = runAndLoad(layout,dat2,ivd);
            %process the varied KM attempt
            if v.km_cel > lastkm
                laststep_km = 1;
            else
                laststep_km = -1;
            end
            rmse = dat2.rmse(end);
            if rmse < lastrmse
                lastrmse = rmse;
                lastimproved_km = 1;
                lastkm = v.km_cel;
            else
                lastimproved_km = -1;
            end
        end
    end
    
    if permuteKcat
        %run the varied KCAT attempt
        v.alpha = lastalpha;
        v.km_cel = lastkm;
        kcatstep = 2^(stepsize * laststep_kcat * lastimproved_kcat);
        v.kcat_cel = lastkcat + kcatstep;
        if isNovel(v,dat2)
            layout = buildDat2Layout(model_default,v);
            dat2 = runAndLoad(layout,dat2,ivd);
            %process the varied KCAT attempt
            if v.kcat_cel > lastkcat
                laststep_kcat = 1;
            else
                laststep_kcat = -1;
            end
            rmse = dat2.rmse(end);
            if rmse < lastrmse
                lastrmse = rmse;
                lastimproved_kcat = 1;
                lastkcat = v.kcat_cel;
            else
                lastimproved_kcat = -1;
            end
        end
    end
    
    if permuteAlpha
        %run the varied ALPHA attempt
        v.kcat_cel = lastkcat;
        v.km_cel = lastkm;
        alphastep = 2^(stepsize * laststep_alpha * lastimproved_alpha);
        v.alpha = lastalpha * alphastep;
        if roundalpha
            logalpha = log2(v.alpha);
            v.alpha = pow2(round(logalpha,sigfigs_alpha));
            if v.alpha == lastalpha
                loglastalpha = log2(lastalpha);
                minstep_alpha = 10^(-sigfigs_alpha);
                if alphastep > 0
                    v.alpha = pow2(loglastalpha + minstep_alpha);
                else
                    v.alpha = pow2(loglastalpha - minstep_alpha);
                end
            end
        end
        
        %Bug/TODO: Step sizes > 1 can sometimes result in negative alphas
        
        %v.initenzyme = v.initialpop * v.alpha;
        if (v.alpha > 0) && isNovel(v,dat2)
            if variableInitEnzyme
                v.initenzyme = v.initialpop * v.alpha;
            end
            layout = buildDat2Layout(model_default,v);
            dat2 = runAndLoad(layout,dat2,ivd);
            %process the varied ALPHA attempt
            if v.alpha > lastalpha
                laststep_alpha = 1;
            else
                laststep_alpha = -1;
            end
            rmse = dat2.rmse(end);
            if rmse < lastrmse
                lastrmse = rmse;
                lastimproved_alpha = 1;
                lastalpha = v.alpha;
            else
                lastimproved_alpha = -1;
            end
        end
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
enzmodel = prepareYeastGEMModel('model',model,...
    'v',v);
%     'alpha',v.alpha,...
%     'vmax_glc',v.vmax_glc,...
%     'maintRxnAsObjective',v.maintRxnAsObjective,...
%     'atp',v.atp_per_peptide,...
%     'gtp',v.gtp_per_peptide,...
%     'addEnzymeToMaint',v.addEnzymeToMaint,...
%     'anaerobic',v.anaerobic,...
%     'vmax_gly',v.vmax_gly,...
%     'vmax_man',v.vmax_man);
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
layout = setMedia(layout,'h[e]',v.initH);
layout = setMedia(layout,'h2o[e]',v.initH2O);
layout = setMedia(layout,'co2[e]',v.initCO2);
layout = setMedia(layout,'pi[e]',v.initPi);
layout = setMedia(layout,'so4[e]',v.initSO4);
layout = setMedia(layout,'nh4[e]',v.initNH4);
layout = setMedia(layout,'k[e]',v.initK);
layout = setMedia(layout,'fe2[e]',v.initfe2);
layout = setMedia(layout,'oleate[e]',v.initoleate);
layout = setMedia(layout,'palmitoleate[e]',v.initpalmitoleate);
layout = setMedia(layout,'glc-D[e]',v.initglc);
layout = setMedia(layout,'gthox[e]',v.initgthox);
layout = setMedia(layout,'man-D[e]',v.initman);


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

function res = loadResult(layout,ivd)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]'};
media = parseMediaLog(layout.params.mediaLogName,medianames);

res.timestep = unique(media.t);
res.t = res.timestep * layout.models{1}.v.timestep; %time in hours

res.glc_amt = media.amt(stridx('glc-D[e]',media.metname));
res.enzyme_amt = media.amt(stridx('enzyme[e]',media.metname));
res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));
%res.gthox_amt = media.amt(stridx('gthox[e]',media.metname));
%res.glycogen_amt = media.amt(stridx('glycogen[c]',media.metname));

res.glc_delta = calcDelta(res.glc_amt);
res.enzyme_delta = calcDelta(res.enzyme_amt);
res.cellulose_delta = calcDelta(res.cellulose_amt);
res.etoh_delta = calcDelta(res.etoh_amt);
%res.gthox_delta = calcDelta(res.gthox_amt);
%res.glycogen_delta = calcDelta(res.glycogen_amt);

biomass = parseBiomassLog(layout.params.biomassLogName);
res.biomass= biomass.biomass;
res.biomass_delta = calcDelta(res.biomass);

[rmse,shift] = findFitQualityWithLag(ivd.denhaan.ctdat(1:8),ivd.denhaan.gdw(1:8),res.t,res.biomass,3,'log10')
res.rmse = rmse;
res.shift = shift;
res.shift_hours = shift * layout.models{1}.v.timestep;
end

function t = applyResult(layout, t, ivd)
%add the layout and its results to the table
res = loadResult(layout, ivd);
tab = table;

tab.layout = {layout};
tab.t = {res.t};
tab.timestep = {res.timestep};
tab.glc_amt = {res.glc_amt};
tab.enzyme_amt = {res.enzyme_amt};
tab.cellulose_amt = {res.cellulose_amt};
tab.etoh_amt = {res.etoh_amt};
%tab.gthox_amt = {res.gthox_amt};
%tab.glycogen_amt = {res.glycogen_amt};
tab.biomass = {res.biomass};
tab.glc_delta = {res.glc_delta};
tab.enzyme_delta = {res.enzyme_delta};
tab.cellulose_delta = {res.cellulose_delta};
tab.etoh_delta = {res.etoh_delta};
%tab.gthox_delta = {res.gthox_delta};
tab.biomass_delta = {res.biomass_delta};
%tab.glycogen_delta = {res.glycogen_delta};

tab.alpha = layout.models{1}.v.alpha;
tab.decayrate = layout.models{1}.v.enzdecayperhour;
tab.costfactor = layout.models{1}.v.costfactor;
tab.kcat_cel = layout.models{1}.v.kcat_cel;
tab.km_cel = layout.models{1}.v.km_cel;
tab.vmax_glc = layout.models{1}.v.vmax_glc;
tab.km_glc = layout.models{1}.v.km_glc;
tab.anaerobic = layout.models{1}.v.anaerobic;
%tab.initgly_pct = layout.models{1}.v.initgly_pct;

tab.rmse = res.rmse;
tab.shift = res.shift;

if size(t,1) > 0
    t = [t;tab];
else
    t = tab;
end
end

function t = runAndLoad(layout,t,ivd)
curdir = pwd;
cd('C:\sync\biomes\cellulose\optima\temp');
runComets(layout);
t = applyResult(layout,t,ivd);
cd(curdir);
end

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
end

function res = isNovel(v,t)
%is the parameter set not yet represented in the table t?
found = false;
if ~isempty(t)
    for i = 1:length(t.layout)
        tv = t.layout{i}.models{1}.v;
        match = isequaln(tv,v);
        found = found | match;
    end
end
res = ~found;
end