function model = changeCompartmentFromBracketsToUnderscore(model)
%CHANGECOMPARTMENTFROMBRACKETSTOUNDERSCORE changes how the compartment is
%designated in the given model's 'mets' field


mets = model.mets;
for i = 1:length(mets)
    met = mets{i};
    if (~isempty(regexp(met,'[\[\(]', 'once')))
        %Greedy search, so there should only be two tokens
        %regardless of the number of brackets in the string
        [tokens,tmp] = regexp(met,'(.+)\[(.+)]','tokens','match');
        %Now replace the name in underscore format
        mets{i} = [tokens{1}{1} '_' tokens{1}{2}];
    end
end

model.mets = mets;



end

