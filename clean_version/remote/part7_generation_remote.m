function layouts = part7_generation_remote(arg)
%PART7_GENERATION_REMOTE
%Generates the layouts to perform sensitivity analysis on
%cellulose-degrading COMETS simulations

v = arg.v;
% if isfield(arg,'variableInitEnzyme')
%     v.variableInitEnzyme = arg.variableInitEnzyme;
% else
%     v.variableInitEnzyme = 0;
% end
% if isfield(arg,'addEnzymeToMaint')
%     v.addEnzymeToMaint = arg.addEnzymeToMaint;
% else
%     v.addEnzymeToMaint = 0;
% end
v_default = v;
factors = arg.factors;
factors = factors(factors ~= 1);
model = arg.model;

layouts = cell(0);

varInitEnz = [0 1];
enzInMaint = [0 1];

for vie = varInitEnz
    for eim = enzInMaint
        v_default.variableInitEnzyme = vie;
        v_default.addEnzymeToMaint = eim;
        v = v_default;
        
        %build the base case (no permutations)
        v.permuted = 'none';
        v.permutation_scale = 1;
        layouts = addcase(layouts,model,v);
        
        %build the permuted cases
        if arg.dims.alpha == 1
            for f = factors
                v = v_default;
                v.alpha = v.alpha * f;
                v.permuted = 'alpha';
                v.permutation_scale = f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.initenz == 1
            for f = factors
                v = v_default;
                v.permuted = 'initial enzyme';
                v.permutation_scale = f;
                v.initenzyme = v.initenzyme * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.alpha_initenz == 1
            for f = factors
                v = v_default;
                v.permuted = 'alpha plus initial enzyme';
                v.permutation_scale = f;
                v.alpha = v.alpha * f;
                v.initenzyme = v.initenzyme * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.kcat_cel == 1
            for f = factors
                v = v_default;
                v.permuted = 'kcat_cel';
                v.kcat_cel = v.kcat_cel * f;
                v.permutation_scale = f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.km_cel == 1
            for f = factors
                v = v_default;
                v.permuted = 'km_cel';
                v.km_cel = v.km_cel * f;
                v.permutation_scale = f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.vmax_glc == 1
            for f = factors
                v = v_default;
                v.permuted = 'vmax_glc';
                v.permutation_scale = f;
                v.vmax_glc = v.vmax_glc * f;
                v.glcuptakerate = v.vmax_glc;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.km_glc == 1
            for f = factors
                v = v_default;
                v.permuted = 'km_glc';
                v.permutation_scale = f;
                v.km_glc = v.km_glc * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.decayrate == 1
            for f = factors
                v = v_default;
                v.permuted = 'decayrate';
                v.permutation_scale = f;
                v.enzdecayperhour = v.enzdecayperhour * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.costfactor == 1
            for f = factors
                v = v_default;
                v.permuted = 'costfactor';
                v.permutation_scale = f;
                v.costfactor = v.costfactor * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.atpperpeptide == 1
            for f = factors
                v = v_default;
                v.permuted = 'atp per peptide';
                v.permutation_scale = f;
                v.atp_per_peptide = v.atp_per_peptide * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.gtpperpeptide == 1
            for f = factors
                v = v_default;
                v.permuted = 'gtp per peptide';
                v.permutation_scale = f;
                v.gtp_per_peptide = v.gtp_per_peptide * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.energyperpeptide == 1
            for f = factors
                v = v_default;
                v.permuted = 'energy per peptide';
                v.permutation_scale = f;
                v.atp_per_peptide = v.atp_per_peptide * f;
                v.gtp_per_peptide = v.gtp_per_peptide * f;
                layouts = addcase(layouts,model,v);
            end
        end
    end
end


end

function newlayouts = addcase(oldlayouts,model,v)
    newlayouts = [oldlayouts {builddat7layout(model,v)}];
end

function layout = builddat7layout(model, v)
% enzmodel = buildScerevisiaeModel(v,model,['SC_A' num2str(v.alpha) 'D' num2str(v.enzdecayperhour)],'fixed');
% enzmodel.v = v;

enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'vmax_glc',v.vmax_glc,...
    'maintRxnAsObjective',v.maintRxnAsObjective,'atp',v.atp_per_peptide,...
    'gtp',v.gtp_per_peptide,'addEnzymeToMaint',v.addEnzymeToMaint,'costfactor',v.costfactor);
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
layout = setMedia(layout,'oleate[e]',v.initoleate);
layout = setMedia(layout,'palmitoleate[e]',v.initpalmitoleate);

layout = setMedia(layout,'glc-D[e]',v.initglc);

layout = setMedia(layout,'gthox[e]',v.initgthox);
% layout.params.defaultVmax = 1000;
% layout.params.defaultKm = 1e-8;

% %%add some O2 dependent metabolites
% %dependentmets = {'pa_SC[c]' 'pc_SC[c]' 'pe_SC[c]' 'ps_SC[c]' 'triglyc_SC[c]'};
% dependentmets = {'triglyc_SC[c]'};
% %layout.models{1} = addExchangeRxn(layout.models{1},dependentmets);
% layout.models{1} = addReaction(layout.models{1},'EX_triglyc_SC[c]',{'triglyc_SC[c]'},-1,'lowerBound',-99999);
% dep_amt = 10;
% layout = setMedia(layout,dependentmets,dep_amt);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.biomassLogRate = 5;
layout.params.mediaLogRate = 5;
%temporary
% layout = setMedia(layout,'glc-D[e]',5);
% layout = setMedia(layout,'o2[e]',5);
end