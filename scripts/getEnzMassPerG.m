function w = getEnzMassPerG(model)
%GETENZMASSPERG calculate how many grams of enzyme are produced per gram
%biomass
%   Detailed explanation goes here

%determine alpha
met = stridx('enzyme[c]',model.mets);
alpha = -min(model.S(met,:));

%get enzyme stoich
rxn = stridx('produce_enzyme',model.rxns);
aminos = {'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'gln-L[c]' 'glu-L[c]'...
    'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' 'phe-L[c]' 'pro-L[c]' ...
    'ser-L[c]' 'thr-L[c]' 'trp-L[c]' 'tyr-L[c]' 'val-L[c]'};
aaidx = zeros(20,1);
for i = 1:20
    aaidx(i) = stridx(aminos{i},model.mets,false);
end
stoich = model.S(aaidx,rxn);

%do the math
aaweights = [71.08 156.2 114.1 115.1 103.1 129.1 128.1 57.05 137.1 113.2 113.2 128.2 131.2 147.2 94.12 87.08 101.1 186.2 163.2 99.13];
w = stoich .* (aaweights' / 1000);
w = w * alpha ;
w = -sum(w);
w = full(w);
end

