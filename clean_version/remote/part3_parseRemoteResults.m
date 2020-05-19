function dat = part3_parseRemoteResults(layouts,batchid)

curdir = pwd;
cd(['C:\sync\biomes\cellulose\optima\temp\' num2str(batchid)]);
dat = table();
for i = 1:length(layouts)
    layout = layouts{i};
    dat = applyResult(layout,dat,i);
    
end



cd(curdir);
fclose all;
end

function res = loadResult(layout,i)
medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]'};
[~,mediafilename,mediafileext] = fileparts(layout.params.mediaLogName);
disp(['Loading Media File ./log_' num2str(i) '/' mediafilename mediafileext]);
media = parseMediaLog(['./log_' num2str(i) '/' mediafilename mediafileext],medianames);
fclose all;

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

[~,biofilename,biofileext] = fileparts(layout.params.biomassLogName);
biomass = parseBiomassLog(['./log_' num2str(i) '/' biofilename biofileext]);
res.biomass= biomass.biomass;
res.biomass_delta = calcDelta(res.biomass);
res.biomass_max = max(biomass.biomass);

res.initcel = layout.models{1}.v.initcellulose;
res.celscale = layout.models{1}.v.cellulosescale;
end

function t = applyResult(layout, t, i)
%add the layout and its results to the table
res = loadResult(layout,i);
tab = table;

if size(res,1) > 0
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
end

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
