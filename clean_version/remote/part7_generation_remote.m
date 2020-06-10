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
varInitEnz = 0;
enzInMaint = 0;


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
        if arg.dims.initman == 1
            for f = factors
                v = v_default;
                v.permuted = 'initial mannose';
                v.permutation_scale = f;
                v.mannose_pct = v.mannose_pct * f;
                v.initman = v.initman * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.ngam == 1
            for f = factors
                v = v_default;
                v.permuted = 'NGAM';
                v.permutation_scale = f;
                v.NGAM = v.NGAM * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.gam == 1
            for f = factors
                v = v_default;
                v.permuted = 'GAM';
                v.permutation_scale = f;
                v.GAM = v.GAM * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.ngam_and_gam == 1
            for f = factors
                v = v_default;
                v.permuted = 'NGAM & GAM';
                v.permutation_scale = f;
                v.NGAM = v.NGAM * f;
                v.GAM = v.GAM * f;
                layouts = addcase(layouts,model,v);
            end
        end
        if arg.dims.initialpop == 1
            for f = factors
                v = v_default;
                v.permuted = 'Initial Pop';
                v.permutation_scale = f;
                v.initialpop = v.initialpop * f;
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

%enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'vmax_glc',v.vmax_glc,...
%    'maintRxnAsObjective',v.maintRxnAsObjective,'atp',v.atp_per_peptide,...
%    'gtp',v.gtp_per_peptide,'addEnzymeToMaint',v.addEnzymeToMaint,'costfactor',v.costfactor);
enzmodel = prepareYeastGEMModel('model',model,'v',v);
enzmodel.v = v;


layout = buildCelOptLayout(enzmodel);
layout = setInitialPop(layout,'1x1',[v.initialpop]);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.biomassLogRate = 5;
layout.params.mediaLogRate = 5;

end