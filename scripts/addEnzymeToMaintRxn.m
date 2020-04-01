function model = addEnzymeToMaintRxn(model,alpha)
%add enzyme production steps to the maintenance reaction of the model.
%Alpha should be given in non-log form.
%Assumes the yeast-GEM model
%ALSO CONVERTS ALPHA IN THE GROWTH REACTION TO SCALE BASED ON GLUCOSE INFLUX
%requires the Cobra Toolbox

%Find fluxes in certain conditions when enzyme is absent
m = model;
growthid = m.rxns{stridx('growth',m.rxnNames,false)};
m = changeRxnMets(m,{'enzyme[c]','enzyme[e]'},{'enzyme[c]','enzyme[e]'},growthid,[0;0]);
idx.growth = stridx('growth',m.rxnNames,false);
idx.maint = stridx('non-growth associated maintenance reaction',m.rxnNames,false);
idx.glc = stridx('D-glucose exchange',m.rxnNames,false);
maintid = m.rxns{idx.maint};
%first, block all unnecessary uptakes
exchanges = findExchRxns(m);
uptakes = (exchanges & (m.lb < 0));
m.lb(uptakes) = 0;
needed = {'ammonium exchange' 'D-glucose exchange' 'ergosterol exchange'...
    'H+ exchange' 'palmitoleate exchange' 'phosphate exchange'...
    'sulphate exchange' 'oleate exchange'};
for i = 1:length(needed)
    m.lb(stridx(needed{i},m.rxnNames,false)) = model.lb(stridx(needed{i},model.rxnNames,false));
end    
%Scale to glc uptake == 1
m.lb(stridx('D-glucose exchange',m.rxnNames,false)) = -1;

%A: glc flux for maintenance only
m.c(idx.growth) = -1;
m.c(idx.maint) = 0;
m.lb(idx.maint) = m.ub(idx.maint);
opt = optimizeCbModel(m);
glc.maint = opt.x(idx.glc);
growth.maint = opt.x(idx.growth);
maint.maint = opt.x(idx.maint);
%B: glc flux for growth with maintenance
m.c(idx.growth) = 1;
opt = optimizeCbModel(m);
glc.full = opt.x(idx.glc);
growth.full = opt.x(idx.growth);
maint.full = opt.x(idx.maint);
%C: glc flux for growth without maintenance
m.lb(idx.maint) = 0;
opt = optimizeCbModel(m);
glc.grow = opt.x(idx.glc);
growth.grow = opt.x(idx.growth);
maint.grow = opt.x(idx.maint);

%-----Trial 1
%if e.g. growth.grow = 0.1704, then when alpha = 0.02 the enz export
%should be 0.1704 * 2% = 0.0034. So if growth.full = 0.1644, which is
%96.48% of growth.grow, the enz export from growth is 0.0034*.9648 and
%the enz export from maintenance is 0.0034*.0352
pctgrowthflux = growth.full / growth.grow;
a.grow = alpha * growth.grow;
a.full = a.grow * pctgrowthflux;
a.maint = a.grow * (1 - pctgrowthflux);

%add the calculated alphas in to the model
model = changeRxnMets(model,{'enzyme[c]','enzyme[e]'},{'enzyme[c]','enzyme[e]'},growthid,[-a.full;a.full]);
model = changeRxnMets(model,{'enzyme[c]','enzyme[e]'},{'enzyme[c]','enzyme[e]'},maintid,[-a.maint;a.maint]);
model.modifications = [model.modifications {'Altered alpha to scale to max glucose uptake in addEnzymeToMaintRxn()'}];
end


% % %-------Storing the old version, which just added an amount of enzyme
% % %equal to alpha ----------
% % maintrxns = find(model.lb > 0);
% % if ~isempty(maintrxns)
% %     for idx = maintrxns
% %        model = changeRxnMets(model,{'enzyme[c]','enzyme[e]'},{'enzyme[c]','enzyme[e]'},model.rxns{idx},[-alpha;alpha]);
% %     end
% %     if ~isfield(model,'modifications')
% %         model.modifications = {};
% %     end
% %     model.modifications = [model.modifications {['Added enzyme to ' num2str(sum(model.lb>0)) ' maintenance reaction(s)']}];
% % end
% % end