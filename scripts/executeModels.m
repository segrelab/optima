function models = executeModels(models)
%if models haven't been run yet, run them and attach the biomass log
for i = 1:length(models)
    if ~isfield(models{i},'biomass')
        models{i} = executeModel(models{i});
    end
end
end