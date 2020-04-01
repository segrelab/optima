function layout = part2_litParams_generation_remote(args)
%PAET2_GENERATION Summary of this function goes here
%   Detailed explanation goes here

v = args.v;

model_default = args.model;

if args.variableInitEnzyme
    v.initenzyme = v.initialpop * v.alpha;
end

layout = buildDat2Layout(model_default,v);

% layout.params.biomassLogRate = 10;
% layout.params.mediaLogRate = 10;

% dat2 = runAndLoad(layout,dat2);
% 
% datadir = 'C:\sync\biomes\cellulose\optima\clean_version\data';
% if ~exist('ivd','var')
%     load([datadir '\ivd.mat']);
% end
% 
% [rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2.t{end},dat2.biomass{end},3,'log10');
% dat2.rmse = rmse;
end



%% functions
function layout = buildDat2Layout(model, v)
%enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
enzmodel = prepareYeastGEMModel('model',model,'v',v);
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
layout = setMedia(layout,'man-D[e]',v.initman);

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


