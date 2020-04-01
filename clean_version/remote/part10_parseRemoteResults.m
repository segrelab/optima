function dat = part10_parseRemoteResults(layouts,batchid)
%DAT6_PARSEREMOTERESULTS load the results obtained by executing COMETS on
%the remote server
%   Detailed explanation goes here

nrows = length(layouts);
curdir = pwd;
workingdir = ['C:\sync\biomes\cellulose\optima\temp\' num2str(batchid)];
cd(workingdir);
dat = table();
for i = 1:nrows
    layout = layouts{i};
    [~,biofilename,biofileext] = fileparts(layout.params.biomassLogName);
    disp(['Loading Biomass ./log_' num2str(i) '/' biofilename biofileext]);
    biomass = parseBiomassLog(['./log_' num2str(i) '/' biofilename biofileext]);
    res = table();
    res.timestep = {unique(biomass.t)};
    res.t = {res.timestep{1} * layout.models{1}.v.timestep}; %time in hours
    if layout.params.writeMediaLog
        medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]' 'gthox[e]'};
        [~,mediafilename,mediafileext] = fileparts(layout.params.mediaLogName);
        disp(['Loading Media ./log_' num2str(i) '/' mediafilename mediafileext]);
        media = parseMediaLog(['./log_' num2str(i) '/' mediafilename mediafileext],medianames);
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
        res.etoh_coord(1) = {etoh_coord};
    end
    [producer_total,producer_coord] = processBiomass(biomass(biomass.model==0,:));
    res.biomass_total(1) = {producer_total};
    res.biomass_coord(1) = {producer_coord};
    
    mu = 0;
    if length(producer_total) > 1
        delta = calcDelta(producer_total);
        delta = delta ./ producer_total(1:end-1);
        t = res.t{1};
        ts = t(2) - t(1);
        mu = delta / ts;
    end
    res.biomass_growthrate(1) = {mu};
        
    if layout.params.writeMediaLog
        res.cellulose_delta(1) = {calcDelta(res.cellulose_total{1})};
        res.glc_delta(1) = {calcDelta(res.glc_total{1})};
        res.enz_delta(1) = {calcDelta(res.enz_total{1})};
        res.etoh_delta(1) = {calcDelta(res.etoh_total{1})};
    end
    
    res.alpha(1) = layout.models{1}.v.alpha;
    res.decayrate(1) = layout.models{1}.v.enzdecayperhour;
    res.costfactor(1) = layout.models{1}.v.costfactor;
    res.kcat_cel(1) = layout.models{1}.v.kcat_cel;
    res.km_cel(1) = layout.models{1}.v.km_cel;
    res.vmax_glc(1) = layout.models{1}.v.vmax_glc;
    res.km_glc(1) = layout.models{1}.v.km_glc;
    res.diff_enz(1) = layout.models{1}.v.diff_enz;
    res.diff_glc(1) = layout.models{1}.v.diff_glc;
    res.enz_init(1) = layout.models{1}.v.initenzyme;
    
    if isfield(layout.models{1}.v,'drainatedges')
        res.drainAtEdges(1) = layout.models{1}.v.drainatedges;
    else
        res.drainAtEdges(1) = false;
    end
    
    dat = vertcat(dat,res);
end

cd(curdir);
end
%% functions

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
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
