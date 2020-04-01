function celtab = refreshCelOptDataTable(celtab)
%REFRESHCELOPTDATATABLE Execute and fill out rows of an existing "celtab"

for i = 1:size(celtab,1)
    row = celtab(i,:);
    %execute if it hasn't already been done
    if ~isfield(row.model{1},'biomass')
        m = executeModel(row.model{1});
        celtab.model{i} = m;
    end

end

final_biomass = cellfun(@(x) x.biomass(end),celtab.model);
celtab.biomass_final = final_biomass;

final_enzyme = cellfun(@(x) x.enzyme_amt(end),celtab.model);
celtab.enzyme_final = final_enzyme;

final_cel = cellfun(@(x) x.cellulose_amt(end),celtab.model);
celtab.cellulose_final = final_cel;

%get the max growth rate, and other model-wise variables
g = zeros(size(celtab,1),1);
deathrates = zeros(size(celtab,1),1);
deathsperhour = zeros(size(celtab,1),1);
for i = 1:length(celtab.model)
    %get the rate a few steps before the curve levels off
    model = celtab.model{i};
    deltas = log(model.biomass(2:end)) - log(model.biomass(1:end-1));
    g(i) = max(deltas);
    if isfield(model.v,'deathrate')
        deathrates(i) = model.v.deathrate;
        deathsperhour = model.v.deathrate / model.v.timestep;
    end
end
celtab.vgrowth_max = g;
celtab.deathrate = deathrates;
%celtab.deathperhour = deathsperhour;

modes = cellfun(@(x) x.enzExpressionMode,celtab.model, 'UniformOutput', false);
celtab.mode = modes;

km_cel = cellfun(@(x) x.v.km_cel,celtab.model);
kcat_cel = cellfun(@(x) x.v.kcat_cel,celtab.model);
km_glc = cellfun(@(x) x.v.km_glc,celtab.model);
vmax_glc = cellfun(@(x) x.v.vmax_glc,celtab.model);
celtab.km_cel = km_cel;
celtab.kcat_cel = kcat_cel;
celtab.km_glc = km_glc;
celtab.vmax_glc = vmax_glc;

xdim = cellfun(@(x) x.v.xdim,celtab.model);
ydim = cellfun(@(x) x.v.ydim,celtab.model);
celtab.xdim = xdim;
celtab.ydim = ydim;

end