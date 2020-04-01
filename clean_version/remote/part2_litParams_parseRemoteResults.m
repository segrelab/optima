function dat2 = part2_litParams_parseRemoteResults(args,layouts,batchid)
%PAET2_GENERATION Summary of this function goes here
%   Detailed explanation goes here

v = args.v;
v_default = v;
dat2 = table();

model_default = args.model;

v = v_default;

datadir = 'C:\sync\biomes\cellulose\optima\clean_version\data';
if ~exist('ivd','var')
    load([datadir '\ivd.mat']);
end

curdir = pwd;
cd(['C:\sync\biomes\cellulose\optima\temp\' num2str(batchid)]);
for i = 1:length(layouts)
    dat2 = applyResult(layouts{i},dat2,i);
    
    %[rss,r2,V,K,rmse] = findFitQuality(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2.t{end},dat2.biomass{end},3,'log10');
    [rmse,shift] = findFitQualityWithLag(ivd.denhaan.ctdat,ivd.denhaan.gdw,dat2.t{end},dat2.biomass{end},3,'log10');
    dat2.rmse(i) = rmse;
    dat2.shift(i) = shift;
end

cd(curdir);
end


%% functions
function res = loadResult(layout,i)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]' 'gthox[e]'};
[~,mediafilename,mediafileext] = fileparts(layout.params.mediaLogName);
media = parseMediaLog(['./log_' num2str(i) '/' mediafilename mediafileext],medianames);

res.timestep = unique(media.t);
res.t = res.timestep * layout.models{1}.v.timestep; %time in hours

res.glc_amt = media.amt(stridx('glc-D[e]',media.metname));
res.enzyme_amt = media.amt(stridx('enzyme[e]',media.metname));
res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));
res.gthox_amt = media.amt(stridx('gthox[e]',media.metname));
%res.glycogen_amt = media.amt(stridx('glycogen[c]',media.metname));

res.glc_delta = calcDelta(res.glc_amt);
res.enzyme_delta = calcDelta(res.enzyme_amt);
res.cellulose_delta = calcDelta(res.cellulose_amt);
res.etoh_delta = calcDelta(res.etoh_amt);
res.gthox_delta = calcDelta(res.gthox_amt);
%res.glycogen_delta = calcDelta(res.glycogen_amt);

[~,biofilename,biofileext] = fileparts(layout.params.biomassLogName);
biomass = parseBiomassLog(['./log_' num2str(i) '/' biofilename biofileext]);
res.biomass= biomass.biomass;
res.biomass_delta = calcDelta(res.biomass);
end

function t = applyResult(layout, t, i)
%add the layout and its results to the table
res = loadResult(layout,i);
tab = table;

tab.layout = {layout};
tab.t = {res.t};
tab.timestep = {res.timestep};
tab.glc_amt = {res.glc_amt};
tab.enzyme_amt = {res.enzyme_amt};
tab.cellulose_amt = {res.cellulose_amt};
tab.etoh_amt = {res.etoh_amt};
tab.gthox_amt = {res.gthox_amt};
%tab.glycogen_amt = {res.glycogen_amt};
tab.biomass = {res.biomass};
tab.glc_delta = {res.glc_delta};
tab.enzyme_delta = {res.enzyme_delta};
tab.cellulose_delta = {res.cellulose_delta};
tab.etoh_delta = {res.etoh_delta};
tab.gthox_delta = {res.gthox_delta};
tab.biomass_delta = {res.biomass_delta};
%tab.glycogen_delta = {res.glycogen_delta};

tab.alpha = layout.models{1}.v.alpha;
tab.decayrate = layout.models{1}.v.enzdecayperhour;
tab.costfactor = layout.models{1}.v.costfactor;
tab.kcat_cel = layout.models{1}.v.kcat_cel;
tab.km_cel = layout.models{1}.v.km_cel;
tab.vmax_glc = layout.models{1}.v.vmax_glc;
tab.km_glc = layout.models{1}.v.km_glc;
%tab.glycogen_pct_weight = layout.models{1}.v.initgly_pct;

tab.rmse = 9999;
tab.shift = 9999;

if size(t,1) > 0
    t = [t;tab];
else
    t = tab;
end
end

function d = calcDelta(x)
%find the difference between every step
d = x(2:end) - x(1:end-1);
end
