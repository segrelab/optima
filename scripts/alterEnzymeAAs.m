function [model, netchange] = alterEnzymeAAs(model,aaname,foldchange,varargin)
%ALTERENZYMEAAS change the AA composition of the enzyme in the model
%without changing the total number of AAs needed
% foldchange = number to multiply the given AA's stoich by
% The second returned argument is the sum of the AA stoichiometries (for
% testing purposes)

trnas = cell(0);
if length(varargin) < 1
    aminos = {'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'gln-L[c]' 'glu-L[c]'...
    'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' 'phe-L[c]' 'pro-L[c]' ...
    'ser-L[c]' 'thr-L[c]' 'trp-L[c]' 'tyr-L[c]' 'val-L[c]'};
else
    aminos = varargin{1};
    trnas = varargin{2};
end

if ~ismember(aaname,aminos)
    error(['AA name not recognized. Should be one of ' aminos{:}]);
end

aminoidxs = zeros(1,20);
for i = 1:20
    aminoidxs(i) = stridx(aminos{i},model.mets,false);
end

targetAAidx = stridx(aaname,model.mets,false);

rxn = stridx('produce_enzyme',model.rxns,false);

sum0 = full(sum(model.S(aminoidxs,rxn)));

s0 = model.S(targetAAidx,rxn);
sdelta = s0 * (foldchange - 1);

fs = zeros(1,20);
for i = 1:20
    if ~strcmp(aminos{i},aaname)
        f = abs(model.S(aminoidxs(i),rxn) / (sum0 - s0));
        newval = model.S(aminoidxs(i),rxn) - (f * sdelta);
        model.S(aminoidxs(i),rxn) = newval;
        if ~isempty(trnas)
            trnaidx = stridx(trnas{i},model.mets,false);
            model.S(trnaidx,rxn) = -newval;
        end
        fs(i) = f;
    end
end

model.S(targetAAidx,rxn) = s0 * foldchange;
if ~isempty(trnas)
    trnaidx = stridx(trnas{stridx(aaname,aminos,false)},model.mets,false);
    model.S(trnaidx,rxn) = -model.S(targetAAidx,rxn);
end

sumf = full(sum(model.S(aminoidxs,rxn)));

%testing: netchange should be 0, ft should be 1
netchange = sumf - sum0;
ft = sum(fs);
end

