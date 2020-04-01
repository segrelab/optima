function model = buildScerevisiaeModel(v,defaultmodel,description,mode)
%Builds a yeast model that secretes cellulase
objname = 'BIOMASS_SC5_notrace';

if isfield(defaultmodel,'modelID') && strcmp(defaultmodel.modelID(1:min(length(defaultmodel.modelID),8)),'yeastGEM')
    objname = defaultmodel.rxns{find(defaultmodel.c)};
end
if isfield(defaultmodel,'modelID') && strcmp(defaultmodel.modelID,'iIN800')
    objname = defaultmodel.rxns{find(defaultmodel.c)};
end
if isfield(defaultmodel,'modelID') && strcmp(defaultmodel.modelID,'yeastGEM_v8.3.2')
    objname = defaultmodel.rxns{find(defaultmodel.c)};
end

model = buildSpeciesModel(v,defaultmodel,description,mode,objname);

% add ergosterol and zymosterol exchange to enable anaerobic growth
%lb = -10;
% model = addReaction(model,'EX_zymst_e',{'zymst[e]'},-1,'lowerBound',lb);
% model = addReaction(model,'EX_zymstest_SC_e',{'zymstest_SC[e]'},-1,'lowerBound',lb);
% model = addReaction(model,'EX_ergst_e',{'ergst[e]'},-1,'lowerBound',lb);
% model = addReaction(model,'EX_ergstest_SC_e',{'ergstest_SC[e]'},-1,'lowerBound',lb);
%model = changeRxnBounds(model,{'EX_zymst_e' 'EX_zymstest_SC_e' 'EX_ergst_e', 'EX_ergstest_SC_e'},lb,'l');


%model = addExchangeRxn(model,{'zymst[e]','zymstest_SC[e]','ergst[e]','ergstest_SC[e]'});

end