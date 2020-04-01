function dat = part1_parseRemoteResults(layouts,batchid)

nrows = length(layouts);
curdir = pwd;
workingdir = ['C:\sync\biomes\cellulose\optima\temp\' num2str(batchid)];
cd(workingdir);
dat = table();
for i = 1:nrows
    layout = layouts{i};
    res = table();
    medianames = {'glc-D[e]' 'enzyme[e]' 'cellulose'};% 'etoh[e]'};
    [~,mediafilename,mediafileext] = fileparts(layout.params.mediaLogName);
    media = parseMediaLog(['./log_' num2str(i) '/' mediafilename mediafileext],medianames);
    
    res.timestep{1} = unique(media.t);
    res.t{1} = res.timestep{1} * layout.models{1}.v.timestep; %time in hours
    
    res.glc_amt{1} = media.amt(stridx('glc-D[e]',media.metname));
    res.enzyme_amt{1} = media.amt(stridx('enzyme[e]',media.metname));
    res.cellulose_amt{1} = media.amt(stridx('cellulose',media.metname));
    %res.etoh_amt = media.amt(stridx('etoh[e]',media.metname));
    
    res.glc_delta{1} = calcDelta(res.glc_amt{1});
    res.enzyme_delta{1} = calcDelta(res.enzyme_amt{1});
    res.cellulose_delta{1} = calcDelta(res.cellulose_amt{1});
    %res.etoh_delta = calcDelta(res.etoh_amt);
    
    [~,biofilename,biofileext] = fileparts(layout.params.biomassLogName);
    biomass = parseBiomassLog(['./log_' num2str(i) '/' biofilename biofileext]);
    res.biomass{1} = biomass.biomass;
    res.biomass_delta{1} = calcDelta(res.biomass{1});
    
    res.decayrate{1} = layout.models{1}.v.enzdecayperhour;
    
    res.enzyme_in_maint = layout.models{1}.v.addEnzymeToMaint;
    res.variable_init_enzyme = layout.models{1}.v.variableInitEnzyme;
    
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
