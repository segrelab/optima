function layouts = part8_generation_remote(arg)
%Generate layouts for simulation that include a range of alphas and a range
%of enzyme decay rates
p.version = 1;

alphas = arg.alphas;
decayscales = arg.decayscales;
v = arg.v;
v.variableInitEnzyme = arg.variableInitEnzyme;
v.enzInMaint = arg.enzInMaint;

v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;

v_default = v;

layouts = cell(0);

model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = model_default.model;

%do the work
for i = 1:length(alphas)
    for j = 1:length(decayscales)
        v = v_default;
        v.alpha = alphas(i);
        v.enzdecayperhour = v.enzdecayperhour * decayscales(j);
        v.decayscale = decayscales(j);
        if v.variableInitEnzyme
            v.initenzyme = v.initialpop * v.alpha;
        end
        layout = builddat8layout(model_default,v);
        layouts = [layouts {layout}];
    end
end
end

%% functions
function layout = builddat8layout(model,v)
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'vmax_glc',v.vmax_glc,'maintRxnAsObjective',v.maintRxnAsObjective,'atp',v.atp_per_peptide,'gtp',v.gtp_per_peptide,'addEnzymeToMaint',v.enzInMaint);
enzmodel.v = v;

layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);
layout = setMedia(layout,layout.mets,0);
layout = setMedia(layout,'cellulose',v.initcellulose);
layout = setMedia(layout,'enzyme[e]',v.initenzyme);
layout = setMedia(layout,'zymst[e]',v.initzymst);
layout = setMedia(layout,'ergst[e]',v.initergst);
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
layout = setMedia(layout,'oleate[e]',v.initoleate);
layout = setMedia(layout,'palmitoleate[e]',v.initpalmitoleate);
layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
end