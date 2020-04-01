function model = buildSpeciesModel(v,defaultmodel,description,mode,objname)
model = defaultmodel;

%add biomass pseudometabolite to the biomass reaction
%model = addMetabolitesToRxn(model,'biomass_pseudo',find(model.c),1);
%mode is one of: {'constitutive' 'proportional' 'demand' 'fixed'}, or the mode's first letter
%create the enzyme reaction
%include ala__d_C?
aminos = {'ala__L_c' 'arg__L_c' 'asn__L_c' 'asp__L_c' 'cys__L_c' 'gln__L_c' 'glu__L_c'...
    'gly_c' 'his__L_c' 'ile__L_c' 'leu__L_c' 'lys__L_c' 'met__L_c' 'phe__L_c' 'pro__L_c' ...
    'ser__L_c' 'thr__L_c' 'trp__L_c' 'tyr__L_c' 'val__L_c'}; %array of amino acid names
% aminos = {'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'gln-L[c]' 'glu-L[c]'...
%    'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' 'pro-L[c]' ...
%    'ser-L[c]' 'thr-L[c]' 'trp-L[c]' 'tyr-L[c]' 'val-L[c]'}; %array of amino acid names

%process the 'mode' argument
if nargin < 4
    if isfield(v,'mode')
        mode = v.mode;
    else
        warning('buildSpeciesModel missing ''mode'' argument. Assuming ''fixed'' mode.');
        mode = 'fixed';
    end
end
mode = lower(mode);
if startsWith(mode,'c')
    mode = 'constitutive';
elseif startsWith(mode,'p')
    mode = 'proportional';
elseif startsWith(mode,'d')
    mode = 'demand';
elseif startsWith(mode,'f')
    mode = 'fixed';
else
    error('The ''mode'' argument to buildSpeciesModel must be one of {''constitutive'', ''proportional'', ''demand'', ''fixed''}');
end


model.v = v;
model.description = description;
model.enzExpressionMode = mode;

%create the enzyme decay reaction
if ~isfield(v,'enzdecaypersec')
    v.enzdecaypersec = 0;
    v.enzdecayperhour = 0;
end

if ~isfield(v,'costfactor') %scale of the enzyme cost to biomass cost
    v.costfactor = 1;
end
costfactor = v.costfactor;

%Create amino acid exchange reactions
createAAExchanges = false;
if createAAExchanges
    examinos = {'ala__L_e' 'arg__L_e' 'asn__L_e' 'asp__L_e' 'cys__L_e' 'gln__L_e' 'glu__L_e'...
        'gly_e' 'his__L_e' 'ile__L_e' 'leu__L_e' 'lys__L_e' 'met__L_e' 'phe__L_e' 'pro__L_e' ...
        'ser__L_e' 'thr__L_e' 'trp__L_e' 'tyr__L_e' 'val__L_e'};
    for i = 1:length(aminos)
        model = addReaction(model,['transport_' aminos{i}(1:end-3)],{aminos{i},examinos{i}},[1,-1],1);
    end
    model = addExchangeRxn(model,examinos);
end


%Make sure rates are capped in a sensible way, in terms of vmax and alpha
glcex = stridx('EX_glc__D_e',model.rxns);
model.lb(glcex) = -abs(v.glcuptakerate);
model = setKm(model,glcex,v.km_glc);
model = setVmax(model,glcex,v.vmax_glc);

%Removed 6/4/18
%%fix positive lb issue
%model.lb(model.lb>0) = 0;

%%%%%   CONSTITUTIVE MODE
%Enzyme production is maximized up to a given rate (denoted by alpha)
%before any resources are allocated to growth.
%Enforced by setting Ex_enzyme upper bound = alpha * maxGlcRate
if strcmp(mode,'constitutive')
    
%     %Set AA makeup for enzyme equal to biomass
%     aminoidxs = zeros(1,length(aminos));
%     for i = 1:length(aminos)
%         aminoidxs(i) = stridx(aminos{i},model.mets,false);
%     end
%     aminocoeffs = model.S(aminoidxs,find(model.c));
%     
%     %rates = aminocoeffs * v.maxEnzymeRate * v.alpha;
%     rates = aminocoeffs;
%     rates = rates * costfactor; %apply cost scaling factor
%     rates(length(rates)+1) = 1;
%     reactants = [aminos 'enzyme_c'];
    [rates, reactants] = getEnzymeStoich(v);
    model = addReaction(model,'produce_enzyme',reactants, rates, 0);
    
    model = addReaction(model,'secrete_enzyme',{'enzyme_c' 'enzyme_e'},[-1 1],0);
    
    model = addExchangeRxn(model,{'enzyme_e'});
    
    %set dual objective reactions:
    objidx = find(model.c);
    model.c(objidx) = 2;
    model.c(stridx('secrete_enzyme',model.rxns)) = 1;
    
    model = setBiomassRxn(model,objidx);

    
    %set the bottleneck
    glcidx = stridx('EX_glc__D_e',model.rxns);
    gb = model.lb(glcidx);
    
    enzidx = stridx('EX_enzyme_e',model.rxns);
    model.ub(enzidx) = v.alpha * -gb; %units are millimoles per GDW per hour
    model.lb(enzidx) = 0; %don't let it run in reverse
    
end

%%%%    PROPORTIONAL MODE
%Enzyme production is rate is a fixed proportion of instantateous growth
%rate. Alpha denotes how many grams of enzyme are produced
if strcmp(mode,'proportional')
    
    enzstoich = v.alpha; %units are in millimoles
    bstoich = (1 - v.alpha); %units are in s^-1
    
    %add a Biomass Pseudometabolite to the biomass reaction
    obj = stridx(objname,model.rxns);
    model = addMetabolitesToRxn(model,{'biomass_pseudo'},model.rxns{obj},bstoich);
    model = setBiomassRxn(model,obj);
%     %create enzyme with Amino Acid cost per gram equal to that of biomass
%     %Set AA makeup for enzyme equal to biomass
%     aminoidxs = zeros(1,length(aminos));
%     for i = 1:length(aminos)
%         aminoidxs(i) = stridx(aminos{i},model.mets,false);
%     end
%     aminocoeffs = model.S(aminoidxs,obj);
%     rates = aminocoeffs;
%     rates = rates * costfactor; %apply cost scaling factor
%     rates(length(rates)+1) = 1;
%     reactants = [aminos 'enzyme_c'];
    [rates, reactants] = getEnzymeStoich(v);
    model = addReaction(model,'produce_enzyme',reactants, rates, 0);
    
    %enzyme export reaction forces the proportionality bt biomass & enz
    
%     if isfield(v,'enzByWeight')
%         if v.enzByWeight
%             enzstoich = enzstoich * v.enz_weight;
%             bstoich = bstoich * v.enz_weight;
%         end
%     end
    model = addReaction(model,'secrete_enzyme',{'enzyme_c','biomass_pseudo','enzyme_e'},[-enzstoich,-bstoich,enzstoich],0);
    model = addExchangeRxn(model,{'enzyme_e'});
    
end

%%%%     DEMAND MODE
% Percent of Carbon consumed that goes towards enzyme production decreases
% linearly as a function of the carbon uptake rate, approaching 0 as glc
% influx approaches its max value * alpha

if strcmp(mode,'demand')
    
%     %create enzyme with Amino Acid cost per gram equal to that of biomass
%     %Set AA makeup for enzyme equal to biomass
%     obj = stridx(objname,model.rxns);
%     aminoidxs = zeros(1,length(aminos));
%     for i = 1:length(aminos)
%         aminoidxs(i) = stridx(aminos{i},model.mets,false);
%     end
%     aminocoeffs = model.S(aminoidxs,obj);
%     rates = aminocoeffs;
%     rates = rates * costfactor; %apply cost scaling factor
%     rates(length(rates)+1) = 1;
%     reactants = [aminos 'enzyme_c'];
    [rates, reactants] = getEnzymeStoich(v);
    model = addReaction(model,'produce_enzyme',reactants, rates, 0);
    model = addReaction(model,'secrete_enzyme',{'enzyme_c','enzyme_e'},[-1,1],0);
    model = addExchangeRxn(model,{'enzyme_e'});
    
    %add a pseudometabolite to the glc transport reactions
    exrxn1 = stridx('GLCtex_copy1',model.rxns);
    model = addMetabolitesToRxn(model,{'demand_pseudo'},{'GLCtex_copy1'},1/v.alpha);
    %exrxn2 = stridx('GLCtex_copy2',model.rxns);
    %model = addMetabolitesToRxn(model,{'demand_pseudo'},exrxn2,1/v.alpha);
    
    %add a pseudometabolite to the biomass rxn
    obj = stridx(objname,model.rxns);
    model = addMetabolitesToRxn(model,{'biomass_pseudo'},model.rxns(obj),1);
    
    %add a reaction to have the pseudos cancel out
    model = addReaction(model,'limit_growth',{'demand_pseudo' 'biomass_pseudo'},[-1 -1],0);
    
    %add the sink reaction for excess demand
    %model = addReaction(model,'sink_demand_pseudo',{'demand_pseudo'},-1,false);
    model = addExchangeRxn(model,{'demand_pseudo'},-1000,0);
    
    %set the new objectives
    %TODO: Move this somewhere other than model.c !
    model.c = zeros(1,length(model.rxns));
    model.c(obj) = 1;
    model.c(stridx('EX_enzyme_e',model.rxns)) = 2;
    %model.c(stridx('EX_demand_pseudo',model.rxns)) = 3;
    
    
end

%%%% FIXED MODE
% alpha denotes the mmol enzyme produced per g growth.
% enzyme secretion produces a pseudometabolite which cancels out with the
% biomass pseudometabolite
if strcmp(mode,'fixed')
%     
%     bstoich = 1; %units are in s^-1
%     
%     %add a Biomass Pseudometabolite to the biomass reaction
%     obj = stridx(objname,model.rxns);
%     model = addMetabolitesToRxn(model,{'biomass_pseudo'},model.rxns{obj},bstoich);
%     model = setBiomassRxn(model,obj);
%     
%     [rates, reactants] = getEnzymeStoich(v);
%     model = addReaction(model,'produce_enzyme',reactants, rates, 0);
%     
%     model = addReaction(model,'secrete_enzyme',{'enzyme_c','enzyme_e','enzyme_pseudo'},[-1,1,1],0);
%     model = addReaction(model,'sink_pseudometabolites',{'enzyme_pseudo','biomass_pseudo'},[-v.alpha,-1],0);
%     model = addExchangeRxn(model,{'enzyme_e'});
%     
    
    %the biomass reaction requires the transport of enzyme_c->enzyme_e
    [rates, reactants] = getEnzymeStoich(v);
    model = addReaction(model,'produce_enzyme',reactants, rates, 0);
    
    %bstoich = 1;
    enzstoich = v.alpha;
    
    %obj = stridx(objname,model.rxns);
    model = addMetabolitesToRxn(model,{'enzyme[c]';'enzyme[e]'},objname,[-enzstoich;enzstoich]);
    model = addExchangeRxn(model,{'enzyme[e]'});
    
end

%the yeastGEM model doesn't include a proper name for glycine, so add this
if any(strcmp('cbs_407[c]',model.mets))
    model = addReaction(model,'renamed_glycine',{'gly[c]','cbs_407[c]'},[-1,1],'reversible',1);
end

% %added 6/20/18
% %Allow models to be 'viable' even if maintenance demand is not met, but
% %still prioritize maintenence
% %This must be run after an explicit biomass reaction is set!
% if any(model.lb > 0)
%     model = addMaintenanceSubreaction(model);
% end
% %Put this method on ice for now -8/27/18

model = refreshCsense(model);
end