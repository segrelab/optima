function model = maintenanceAsObjective(model)
%Remove any lower bounds greater than 0
%then set the objective field to maximize those reactions
%and set the biomass reaction to the original objective
%update 9/18/19: add a maintenance pseudometabolite so you can track
%whether maintenance is being satisfied
maintrxns = find(model.lb > 0);
if ~isempty(maintrxns)
    model.ub(maintrxns) = model.lb(maintrxns); %probably not necessary, but in case these were not ==min
    model.lb(maintrxns) = 0;
    model.description = [model.description '_noMaintRxn'];
    if ~isfield(model,'modifications')
        model.modifications = cell(0);
    end
    model.modifications = [model.modifications {'Removed maintenance reaction lower bounds where lb>0'}];
    
    baseobjs = find(model.c);
    %set the maint rxns as highest priority
    model.objective = zeros(1,length(model.rxns));
    model.objective(maintrxns) = 1:length(maintrxns);
    %set the original objectives as lower priority (high absolute value is low priority)
    model.objective(baseobjs) = model.c(baseobjs) * (length(maintrxns) + 1);
    model.modifications = [model.modifications {'Added maintenance reaction as primary objective function'}];
    
    %set the biomass reaction
    if ~isfield(model,'biomassRxn')
        model.biomassRxn = baseobjs;
    end
end

end

