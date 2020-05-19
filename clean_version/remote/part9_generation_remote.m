function layouts = part9_generation_remote(arg)
%Generate layouts that knock out each reaction in the model
p.version = 2;
%v2: add disableUptakeOnlyForExch: if true, still allow secretion of the
%exchange metabolites. If false, reaction is completely blocked. V1 had the
%behavior of having this set to false

v = arg.v;
v.variableInitEnzyme = arg.variableInitEnzyme;
v.enzInMaint = arg.enzInMaint;
v.disableUptakeOnlyForExch = arg.disableUptakeOnlyForExch;
v.writemedialog = arg.recordMedia;
v.writebiomasslog = true;
v.writefluxlog = false;


v_default = v;

model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = model_default.model;

lowerIdx = arg.lowerIdx;
upperIdx = min(arg.upperIdx, length(model_default.rxns));

if v_default.variableInitEnzyme
    v_default.initenzyme = v_default.initialpop * v_default.alpha;
end
layout = builddat9layout(model_default,v_default);

indexes = lowerIdx:upperIdx;
if strcmpi(arg.mode,'set')
    indexes = arg.idxSet;
end

layouts = cell(length(indexes),1);

exchanges = findExchRxns(layout.models{1});


%do the work
for i = 1:length(indexes)
    idx = indexes(i);
    templayout = layout;
    isExchange = exchanges(idx);
    if (~isExchange || ~v.disableUptakeOnlyForExch)
        templayout.models{1}.lb(idx) = 0;
        templayout.models{1}.ub(idx) = 0;
    else %is exchange. So only disable uptake
        directionIsPos = sum(templayout.models{1}.S(:,idx)) <= 0;
        if directionIsPos
            templayout.models{1}.ub(idx) = 0;
        else
            templayout.models{1}.lb(idx) = 0;
        end
    end
    templayout.models{1}.v.knockoutIdx = idx;
    templayout.models{1}.modifications = [templayout.models{1}.modifications {['Knock out rxn ' num2str(idx)]}];
    layouts{i} = templayout;
end
layout.models{1}.v.knockoutIdx = 0;
layouts = [{layout}; layouts];
end
%% functions
function layout = builddat9layout(model,v)
%enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'vmax_glc',v.vmax_glc,'maintRxnAsObjective',v.maintRxnAsObjective,'atp',v.atp_per_peptide,'gtp',v.gtp_per_peptide,'addEnzymeToMaint',v.enzInMaint);
enzmodel = prepareYeastGEMModel('model',model,'v',v);
enzmodel.v = v;

layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
end