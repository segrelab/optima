function model = addMaintenanceSubreaction(model)
%ADDMAINTENANCESUBREACTION modify the handling of maintenance reactions so
%that the model will still be viable if maintenance demand cannot be met
%
%Arguments:
%   model: a COBRA model
%Returns:
%   model: a COBRA model
%
%Approach:
%	• If not already set, assign the curren objective as the Biomass Rxn
%	• Create a new Maintenance Reaction "MR" with the same stoichiometry as R, 
%       plus a "maintenance_pseudometabolite" product because COBRA will 
%       not let you create redundant reactions
% 	• Create a sink reaction for the maintenance_pseudometabolite
% 	• Set the maximum flux of MR to the minimum flux of R
% 	• Set the minimum flux of MR to 0
% 	• Set the objective of the model so that the first objective is maximization of MR
% 	• Set the maximum flux of  R to its original maximum minus the maximum of MR
% 	• Set the minimum flux of R to 0
% %
% $Author: mquintin $	$Date: 2018/6/20 $	$Revision: 1 $
% Last edit: mquintin 2018/6/20
% Copyright: Daniel Segrè, Boston University Bioinformatics Program 2018

maintidxs = find(model.lb > 0);%find all reactions that qualify as maintenance
%maintidxs = [];

if ~isfield(model,'objective')
    model.objective = model.c;
end

if ~isfield(model,'biomassRxn')
    %assume it's the maximized reaction with highest priority
    bio = model.objective == min(model.objective(model.objective > 0));
    model.biomassRxn = bio;
end

for i = 1:length(maintidxs)
    idx = maintidxs(i);
    
    %create the new reaction
    stoich = model.S(:,idx);
    stoich = stoich(find(stoich)); %just keep the nonzero elements
    metidxs = find(model.S(:,idx));
    metids = model.mets(metidxs);
    %pseudometname =['maint_pseudo_' num2str(idx)];
    %metids = [metids; {pseudometname}];
    %stoich = [stoich; 1];
    bound = model.lb(idx);
    rxnname = ['maintenance_' num2str(idx)];
    rev = 0;
    if isfield(model,'rev')
        rev = model.rev(idx);
    end
    model = addReaction(model,rxnname,'metaboliteList',metids,'stoichCoeffList',stoich,'reversible',rev,'upperBound',bound);
    
    %create a sink for the psuedometabolite
    %model = addReaction(model,['SINK_' pseudometname],'metaboliteList',{pseudometname},'stoichCoeffList',-1,'reversible',false);
    
    %shift all existing objectives' priority
    positiveObjectives = find(model.objective > 0);
    model.objective(positiveObjectives) = model.objective(positiveObjectives) + 1;
    negativeObjectives = find(model.objective < 0);
    model.objective(negativeObjectives) = model.objective(negativeObjectives) - 1;
    %update the model objective
    rxnidx = stridx(rxnname,model.rxns);
    model.c(rxnidx) = 1;
    model.objective(rxnidx) = 1;
    
    %update the flux bounds of the original reaction
    model.lb(idx) = 0;
    model.ub(idx) = model.ub(idx) - bound;
    
end

end

