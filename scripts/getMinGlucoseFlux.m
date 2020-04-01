function res = getMinGlucoseFlux(model)
%GETMINGLUCOSEFLUX Find the minimum glucose uptake rate (in mmol/gCDW*h)
%required to make the given Cobra model viable

%Change the objective
%Remember that exchange reactions are set so that flux < 0 is uptake. So
%minimizing uptake means maximizing the reaction
model.c(1:end) = 0;
glcex = stridx('EX_glc__D_e',model.rxns);
model.c(glcex) = 1;

if ~any(model.c)
    error('Could not find the glucose exchange reaction (''EX_glc__D_e'').');
end

%optimize
changeCobraSolver('gurobi');
opt = optimizeCbModel(model);

%find glc flux
res = -opt.x(glcex);

end

