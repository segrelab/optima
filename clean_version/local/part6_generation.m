function [res,layouts] = part6_generation(arg)
%PART6_GENERATION 
%Creates a 2d layout with a producer and a cheater colony
%The layout is 11 x 60
%The producer colony is at position (6,6)
%The cheater is at position (6,6+x)


%TODO: For each biomass/media record, create two colummns. One has the
%total amount in the plate, the other is x:y:t matrix showing concentration
%at each point in each timestep

xdim = 11;
ydim = 30;
distances = arg.distances; %distance between the cheater and producer, in cm
%distances = [0,2];
%distances = 5;
enzposn = [6,6]; %location of the enzyme producer

v = arg.v;
model = arg.model;
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'vmax_glc',v.vmax_glc,'maintRxnAsObjective',v.maintRxnAsObjective);
enzmodel.v = v;

distances_cells = floor(distances / v.spaceWidth); %distance in cells instead of cm

cheatmodel = prepareCheaterModel(model,v.vmax_glc,v.maintRxnAsObjective);

layout = createLayout();
layout = setDims(layout,xdim,ydim);

drainAtEdges = v.drainatedges; %should enzyme and glc be removed at edges (as opposed to being allowed to bounce back)
if drainAtEdges
    mask = zeros(xdim,ydim);
    mask(1:xdim,1) = 1;
    mask(1:xdim,ydim) = 1;
    mask(1,1:ydim) = 1;
    mask(xdim,1:ydim) = 1;
    layout = applyStaticMediaMask(layout,mask,{'enzyme[e]','glc-D[e]'},[0 0]);
end
layout = addModel(layout,enzmodel);
layout = addModel(layout,cheatmodel);

layout.initial_pop(1,enzposn(1),enzposn(2)) = v.initialpop;

layout = setMedia(layout,layout.mets,0);

layout = setMedia(layout,'cellulose',v.initcellulose);
layout = setInitialMediaInCell(layout,enzposn(1),enzposn(2),'enzyme[e]',v.initenzyme);
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
layout = setMedia(layout,'glc-D[e]',v.initglc);
layout = setMedia(layout,'gthox[e]',v.initgthox);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.maxSpaceBiomass = v.maxspacebiomass;
layout.params.writeBiomassLog = v.writebiomasslog;
layout.params.writeMediaLog = v.writemedialog;
layout.params.writeFluxLog = v.writefluxlog;
layout.params.biomassLogRate = v.lograte;
layout.params.mediaLogRate = v.lograte;
layout.params.fluxLogRate = v.lograte;
layout.params.biomassLogName = [v.filepath '\log_biomass.m'];
layout.params.mediaLogName = [v.filepath '\log_media.m'];
layout.params.fluxLogName = [v.filepath '\log_flux.m'];
layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.pixelScale = 10;
layout.params.spaceWidth = v.spaceWidth;
layout.params.deathRate = v.deathrate;


layout = addExternalReaction(layout,'degrade_cellulose',{'cellulose','glc-D[e]'},[-1 v.glcpercellulose],'enzyme', 'enzyme[e]', 'k', v.kcat_cel/v.glcpercellulose, 'km', v.km_cel);

if (isfield(v,'enzdecayperhour') && ~isfield(v,'enzdecaypersec'))
    v.enzdecaypersec = v.enzdecayperhour/3600;
end
%create the enzyme decay reaction
if isfield(v,'enzdecaypersec')
    if v.enzdecaypersec ~= 0
        layout = addExternalReaction(layout,'enzyme_decay','enzyme[e]',-1,'vmax',v.enzdecaypersec);
    end
end
if isfield(v,'numExRxnSubsteps')
    layout.params.numExRxnSubsteps = v.numExRxnSubsteps;
end

%modify diffusion rates
layout = setDiffusion(layout,'cellulose',v.diff_cel);
layout = setDiffusion(layout,'enzyme[e]',v.diff_enz);
layout = setDiffusion(layout,'glc-D[e]',v.diff_glc);

res = table();
layouts = cell(0);
curdir = pwd;
cd 'C:\sync\biomes\cellulose\optima\temp';
for i = distances_cells
    l = layout;
    l.initial_pop(2,enzposn(1),enzposn(2)+i) = v.initialpop;
    runComets(l);
    r = loadResult(l);
    r.distance(1) = i;
    res = vertcat(res,r);
    layouts = [layouts {l}];
end
cd(curdir);



end

%% functions

function res = loadResult(layout)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]' 'gthox[e]'};
media = parseMediaLog(layout.params.mediaLogName,medianames);
biomass = parseBiomassLog(layout.params.biomassLogName);
res = table();
res.timestep = {unique(media.t)};
res.t = {res.timestep{1} * layout.models{1}.v.timestep}; %time in hours
[cellulose_total,cellulose_coord] = processMedia(media(strcmp('cellulose',media.metname),:));
res.cellulose_total(1) = {cellulose_total};
res.cellulose_coord(1) = {cellulose_coord};
[glc_total,glc_coord] = processMedia(media(strcmp('glc-D[e]',media.metname),:));
res.glc_total(1) = {glc_total};
res.glc_coord(1) = {glc_coord};
[enz_total,enz_coord] = processMedia(media(strcmp('enzyme[e]',media.metname),:));
res.enz_total(1) = {enz_total};
res.enz_coord(1) = {enz_coord};
[etoh_total,etoh_coord] = processMedia(media(strcmp('etoh[e]',media.metname),:));
res.etoh_total(1) = {etoh_total};
res.etho_coord(1) = {etoh_coord};

[producer_total,producer_coord] = processBiomass(biomass(biomass.model==0,:));
res.producer_total(1) = {producer_total};
res.producer_coord(1) = {producer_coord};
[cheater_total,cheater_coord] = processBiomass(biomass(biomass.model==1,:));
res.cheater_total(1) = {cheater_total};
res.cheater_coord(1) = {cheater_coord};

res.cellulose_delta(1) = {calcDelta(res.cellulose_total{1})};
res.glc_delta(1) = {calcDelta(res.glc_total{1})};
res.enz_delta(1) = {calcDelta(res.enz_total{1})};
res.etoh_delta(1) = {calcDelta(res.etoh_total{1})};

res.alpha(1) = layout.models{1}.v.alpha;
res.decayrate(1) = layout.models{1}.v.enzdecayperhour;
res.costfactor(1) = layout.models{1}.v.costfactor;
res.kcat_cel(1) = layout.models{1}.v.kcat_cel;
res.km_cel(1) = layout.models{1}.v.km_cel;
res.vmax_glc(1) = layout.models{1}.v.vmax_glc;
res.km_glc(1) = layout.models{1}.v.km_glc;
res.diff_enz(1) = layout.models{1}.v.diff_enz;
res.diff_glc(1) = layout.models{1}.v.diff_glc;

end


% function t = applyResult(layout, t)
% %add the layout and its results to the table
% res = loadResult(layout);
% tab = table;
% 
% tab.layout = {layout};
% tab.t = {res.t};
% tab.timestep = {res.timestep};
% tab.glc_amt = {res.glc_amt};
% tab.enzyme_amt = {res.enzyme_amt};
% tab.cellulose_amt = {res.cellulose_amt};
% tab.etoh_amt = {res.etoh_amt};
% tab.gthox_amt = {res.gthox_amt};
% tab.biomass = {res.biomass};
% tab.glc_delta = {res.glc_delta};
% tab.enzyme_delta = {res.enzyme_delta};
% tab.cellulose_delta = {res.cellulose_delta};
% tab.etoh_delta = {res.etoh_delta};
% tab.gthox_delta = {res.gthox_delta};
% tab.biomass_delta = {res.biomass_delta};
% 
% tab.alpha = layout.models{1}.v.alpha;
% tab.decayrate = layout.models{1}.v.enzdecayperhour;
% tab.costfactor = layout.models{1}.v.costfactor;
% tab.kcat_cel = layout.models{1}.v.kcat_cel;
% tab.km_cel = layout.models{1}.v.km_cel;
% tab.vmax_glc = layout.models{1}.v.vmax_glc;
% tab.km_glc = layout.models{1}.v.km_glc;
% 
% if size(t,1) > 0
%     t = [t;tab];
% else
%     t = tab;
% end
% end

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

function model = prepareCheaterModel(model,vmax_glc,maintRxnAsObjective)
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
model.mets{337} = 'biomass[c]';
model.lb(stridx('D-glucose exchange' ,model.rxnNames)) = -vmax_glc;

if maintRxnAsObjective
    model = maintenanceAsObjective(model);
end
end

function [total,coord,ut] = processMedia(tab)
%given just the part of the media table pertaining to a single metabolite
ut = unique(tab.t);
n = length(ut);
ux = unique(tab.x);
uy = unique(tab.y); 
total = zeros(n,1);
coord = zeros(length(ux),length(uy),n);
for j = 1:length(ut)
    i = ut(j);
    subtab = tab(tab.t == i,:);
    idx = find(ut==i);
    for k = 1:length(ux)
        x = ux(k);
        for l = 1:length(uy)
            y = uy(l);
            row = subtab(subtab.x == x,:);
            row = row(row.y == y,:);
            row = row(row.t == i,:);
            coord(x,y,idx) = row.amt;
        end
    end
    total(idx) = sum(subtab.amt(subtab.t == i));
end
end

function [total,coord,ut] = processBiomass(tab)
%given just the part of the media table pertaining to a single model
ut = unique(tab.t);
n = length(ut);
ux = unique(tab.x);
uy = unique(tab.y); 
total = zeros(n,1);
coord = zeros(length(ux),length(uy),n);
for j = 1:length(ut)
    i = ut(j);
    subtab = tab(tab.t == i,:);
    idx = find(ut==i);
    for k = 1:length(ux)
        x = ux(k);
        for l = 1:length(uy)
            y = uy(l);
            row = subtab(subtab.x == x,:);
            row = row(row.y == y,:);
            row = row(row.t == i,:);
            coord(x,y,idx) = row.biomass;
        end
    end
    total(idx) = sum(subtab.biomass(subtab.t == i));
end
end
