function dat = part9_parseRemoteResults(layouts,batchid)

curdir = pwd;
cd(['C:\sync\biomes\cellulose\optima\temp\' num2str(batchid)]);
dat = table();
for i = 1:length(layouts)
    layout = layouts{i};
    
    [~,biofilename,biofileext] = fileparts(layout.params.biomassLogName);
    biopath = ['./log_' num2str(i) '/' biofilename biofileext];
    if exist(biopath,'file')
        dat = applyResult(layout,dat,i);
    else
        biopath2 = ['./log_' num2str(i) '/biomass.m'];
        if exist(biopath2,'file')
            layout.params.biomassLogName = 'biomass.m';
            layout.params.mediaLogName = 'media.m';
            dat = applyResult(layout,dat,i);
        else
            disp (['Biomass file not found: ' biopath]);
        end
    end
    
end

disp(['Loaded ' num2str(size(dat,1)) ' of ' num2str(length(layouts)) ' records.']);

cd(curdir);
end

function res = loadResult(layout,i)
[~,biofilename,biofileext] = fileparts(layout.params.biomassLogName);
disp(['Loading Biomass File ./log_' num2str(i) '/' biofilename biofileext]);
biomass = parseBiomassLog(['./log_' num2str(i) '/' biofilename biofileext]);
res.timestep = unique(biomass.t);
res.t = res.timestep * layout.models{1}.v.timestep; %time in hours
res.biomass = biomass.biomass;
res.biomass_delta = calcDelta(res.biomass);
res.biomass_max = max(biomass.biomass);
if layout.params.writeMediaLog
    medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose' 'etoh[e]'};
    [~,mediafilename,mediafileext] = fileparts(layout.params.mediaLogName);
    disp(['Loading Media File ./log_' num2str(i) '/' mediafilename mediafileext]);
    media = parseMediaLog(['./log_' num2str(i) '/' mediafilename mediafileext],medianames);
    
    res.glc_amt = media.amt(stridx('glc-D[e]',media.metname));
    res.enzyme_amt = media.amt(stridx('enzyme[e]',media.metname));
    res.cellulose_amt = media.amt(stridx('cellulose',media.metname));
    res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));
    
    res.glc_delta = calcDelta(res.glc_amt);
    res.enzyme_delta = calcDelta(res.enzyme_amt);
    res.cellulose_delta = calcDelta(res.cellulose_amt);
    res.etoh_delta = calcDelta(res.etoh_amt);
end
res.knockoutIdx = layout.models{1}.v.knockoutIdx;
end

function t = applyResult(layout, t, i)
%add the layout and its results to the table
res = loadResult(layout,i);
tab = table;

if size(res,1) > 0
    tab.layout = {layout};
    tab.t = {res.t};
    tab.timestep = {res.timestep};
    tab.biomass = {res.biomass};
    tab.biomass_delta = {res.biomass_delta};
    tab.biomass_max = res.biomass_max;
    if layout.params.writeMediaLog
        tab.glc_amt = {res.glc_amt};
        tab.enzyme_amt = {res.enzyme_amt};
        tab.cellulose_amt = {res.cellulose_amt};
        tab.etoh_amt = {res.etoh_amt};
        tab.glc_delta = {res.glc_delta};
        tab.enzyme_delta = {res.enzyme_delta};
        tab.cellulose_delta = {res.cellulose_delta};
        tab.etoh_delta = {res.etoh_delta};
    end
    
    tab.alpha = layout.models{1}.v.alpha;
    tab.decayrate = layout.models{1}.v.enzdecayperhour;
    tab.costfactor = layout.models{1}.v.costfactor;
    tab.kcat_cel = layout.models{1}.v.kcat_cel;
    tab.km_cel = layout.models{1}.v.km_cel;
    tab.vmax_glc = layout.models{1}.v.vmax_glc;
    tab.km_glc = layout.models{1}.v.km_glc;
    
    tab.knockoutIdx = res.knockoutIdx;
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
