%alter the Amino Acid composition for the enzyme, and get the change in
%max objective flux

foldchanges = 0:.1:5; 
foldchanges = [0 1/5 1/3 .5 1.5 2 3 5]; 

aminos = {'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'gln-L[c]' 'glu-L[c]'...
    'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' 'phe-L[c]' 'pro-L[c]' ...
    'ser-L[c]' 'thr-L[c]' 'trp-L[c]' 'tyr-L[c]' 'val-L[c]'};

aminos = {'Ala-tRNA(Ala) [cytoplasm]' 'Arg-tRNA(Arg) [cytoplasm]' ...
    'Asn-tRNA(Asn) [cytoplasm]' 'Asp-tRNA(Asp) [cytoplasm]' 'Cys-tRNA(Cys) [cytoplasm]' 'Gln-tRNA(Gln) [cytoplasm]' ...
    'Glu-tRNA(Glu) [cytoplasm]' 'Gly-tRNA(Gly) [cytoplasm]' 'His-tRNA(His) [cytoplasm]' 'Ile-tRNA(Ile) [cytoplasm]' ...
    'Leu-tRNA(Leu) [cytoplasm]' 'Lys-tRNA(Lys) [cytoplasm]' 'Met-tRNA(Met) [cytoplasm]' 'Phe-tRNA(Phe) [cytoplasm]' ...
    'Pro-tRNA(Pro) [cytoplasm]' 'Ser-tRNA(Ser) [cytoplasm]' 'Thr-tRNA(Thr) [cytoplasm]' 'Trp-tRNA(Trp) [cytoplasm]' ...
    'Tyr-tRNA(Tyr) [cytoplasm]' 'Val-tRNA(Val) [cytoplasm]'};

trna_free_names = {'tRNA(Ala) [cytoplasm]' 'tRNA(Arg) [cytoplasm]' 'tRNA(Asn) [cytoplasm]' ...
    'tRNA(Asp) [cytoplasm]' 'tRNA(Cys) [cytoplasm]' 'tRNA(Gln) [cytoplasm]' 'tRNA(Glu) [cytoplasm]' ...
    'tRNA(Gly) [cytoplasm]' 'tRNA(His) [cytoplasm]' 'tRNA(Ile) [cytoplasm]' 'tRNA(Leu) [cytoplasm]' ...
    'tRNA(Lys) [cytoplasm]' 'tRNA(Met) [cytoplasm]' 'tRNA(Phe) [cytoplasm]' 'tRNA(Pro) [cytoplasm]' ...
    'tRNA(Ser) [cytoplasm]' 'tRNA(Thr) [cytoplasm]' 'tRNA(Trp) [cytoplasm]' 'tRNA(Tyr) [cytoplasm]' ...
    'tRNA(Val) [cytoplasm]'};

%model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904 
model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model_default = model_default.model;

p.alpha = exp(-3.25); 
p.mode = {'Fixed'};%{'Constitutive','Proportional','Fixed};
p.costfactor = 1;

%%STRATEGY VARIABLES
v.defaultbound = 10; %default exchange upper bound
%v.maxEnzymeRate = 1; %maximum enzyme:biomass ratio by weight for fixed production. 1:1 w/ alpha=1 should mean half of all amino acids go into enzyme
v.alpha = p.alpha; %fixed percentage of the maximum enzyme rate being expressed

%%ENZYME VARIABLES. See https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.4
%values used: wild-type T.ressei on 4-methylumbelliferyl cellobioside. https://www.brenda-enzymes.org/literature.php?e=3.2.1.4&r=737187
v.kcat_cel = 27; %mmol converted / mmol enzyme / second  %TODO: Confirm units
v.km_cel = 0.03; %Units = mmoles/cm^3
v.enz_weight =  8.5853183e-20; %grams
v.enzByWeight = false;
v.atp_per_peptide = 2; %used to calculate energy costs in getEnzymeStoich
v.gtp_per_peptide = 2; %see https://www.ncbi.nlm.nih.gov/books/NBK224633/
v.diff_enz = 0;

%%Glucose variables
v.km_glc = 0.01; %millimoles/cm^3.
v.vmax_glc = 10;
v.glcpercellulose = 218.2;%how many units of glucose get released per unit cellulose? 1:1 assumes units are by weight
%Avicel degree of polymerization cite : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3217711/
%v.biomassperglc = 1/(1.447e9); %how many units of biomass are produced for every unit glc consumed?
%weight of E coli / weight of glc molecule = 433fg / (2.992e-22g) = 4.33e-13 / (2.992e-22) = 1.447e9
v.glcuptakerate = -10; %bound for the GLC import reaction
%This should be higher than required for max enzyme production

%Enzyme variables v 2- set high to check that everything works OK
%v.kat = 100;
%v.km = 0.0001;

v.costfactor = p.costfactor;

%print values
p = p;
v_default = v;

basestoich = getEnzymeStoich(v);
basestoich = basestoich(1:20);
seqlength = -sum(basestoich);

%model = buildScerevisiaeModel(v_default,model_default,'AA_model','fixed');
model = prepareYeastGEMModel('model',model_default,'alpha',v.alpha);

enzIdx = stridx('produce_enzyme',model.rxns);

aaIdx = zeros(1,20);
for i = 1:20
    aaIdx(i) = stridx(aminos{i},model.mets,false);
end

models = cell((20 * length(foldchanges)) + 1,1);
objtab = table();
nrows = 20 * length(foldchanges);
objtab.model = cell(nrows,1);
objtab.foldchange = zeros(nrows,1);
objtab.AA = cell(nrows,1);
objtab.costfactor = zeros(nrows,1);
objtab.max_obj_flux = zeros(nrows,1);
objtab.alpha = zeros(nrows,1);
objtab.log_alpha = zeros(nrows,1);

changeCobraSolver('gurobi');

%increase each AA
currow = 0;
for i = 1:20
    for foldchange = foldchanges
        currow = currow + 1;
        m = alterEnzymeAAs(model,aminos{i},foldchange,aminos,trna_free_names);
        objtab.model{currow} = m;
        objtab.foldchange(currow) = foldchange;
        objtab.AA(currow) = aminos(i);
        objtab.costfactor(currow) = p.costfactor;
        opt = optimizeCbModel(m);
        objtab.max_obj_flux(currow) = opt.f;
        objtab.alpha(currow) = p.alpha;
        objtab.log_alpha(currow) = log(p.alpha);
    end
end

%last row is the default
currow = currow + 1;
objtab.model{currow} = model;
objtab.foldchange(currow) = 1;
objtab.AA{currow} = 'none';
objtab.costfactor(currow) = p.costfactor;
opt = optimizeCbModel(model);
objtab.max_obj_flux(currow) = opt.f;
objtab.alpha(currow) = p.alpha;
objtab.log_alpha(currow) = log(p.alpha);

default_flux = objtab.max_obj_flux(end);
objtab.flux_ratio = objtab.max_obj_flux(1:end) / default_flux;
