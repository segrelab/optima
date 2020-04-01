%alter the Amino Acid composition for the enzyme, and get the change in
%yield

%Last updat: 10/25/18: Change params to mactch results from s2

foldchange = [0 1/5 1/3 .5 1.5 2 3 5]; 

% aminos_default = {'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'gln-L[c]' 'glu-L[c]'...
%     'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' 'phe-L[c]' 'pro-L[c]' ...
%     'ser-L[c]' 'thr-L[c]' 'trp-L[c]' 'tyr-L[c]' 'val-L[c]'};
aminos_default = {'Ala-tRNA(Ala) [cytoplasm]' 'Arg-tRNA(Arg) [cytoplasm]' ...
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

p.enzdecayperhour = 0.01;%1/600;
p.alpha = exp(-3.25); %find where alpha collapses
p.deathrate = 0;
p.dilutionrate = 0; %fraction of media replaced per hour
p.mode = {'Fixed'};%{'Constitutive','Proportional','Fixed};
p.layouttype = {'1x1'}; %{'1x1','chemostat'};
p.diff_enz = 1e-6; %Default is 1e-6;
p.costfactor = 1;
v.richmedia = false; %provide a bunch of everything

%%STRATEGY VARIABLES
v.defaultbound = 10; %default exchange upper bound
%v.maxEnzymeRate = 1; %maximum enzyme:biomass ratio by weight for fixed production. 1:1 w/ alpha=1 should mean half of all amino acids go into enzyme
v.alpha = p.alpha; %fixed percentage of the maximum enzyme rate being expressed

%%ENZYME VARIABLES. See https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.4
%values used: wild-type T.ressei on 4-methylumbelliferyl cellobioside. https://www.brenda-enzymes.org/literature.php?e=3.2.1.4&r=737187
v.kcat_cel = 27; %mmol converted / mmol enzyme / second  
v.km_cel = 0.03; %Units = mmoles/cm^3
v.enz_weight =  8.5853183e-20; %grams
v.enzByWeight = false;
v.atp_per_peptide = 2; %used to calculate energy costs in getEnzymeStoich
v.gtp_per_peptide = 2; %see https://www.ncbi.nlm.nih.gov/books/NBK224633/
v.diff_enz = 0;

%%Glucose variables
v.km_glc = v.km_glc / 10;
v.vmax_glc = v.vmax_glc * 1;
v.glcpercellulose = 218.2;%how many units of glucose get released per unit cellulose? 1:1 assumes units are by weight
%Avicel degree of polymerization cite : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3217711/
%v.biomassperglc = 1/(1.447e9); %how many units of biomass are produced for every unit glc consumed?
%weight of E coli / weight of glc molecule = 433fg / (2.992e-22g) = 4.33e-13 / (2.992e-22) = 1.447e9
v.glcuptakerate = v.vmax_glc; %bound for the GLC import reaction
%This should be higher than required for max enzyme production

%Enzyme variables v 2- set high to check that everything works OK
%v.kat = 100;
%v.km = 0.0001;

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0; %Diffusion constants
%v.diff_enz = 1e-6; %Default is 1e-6;
v.diff_glc = 1e-6;
%v.enzdecayrate = 0; %degradation of enzyme over time
%v.deathrate = 0;
%v.enzymeperglc = 0.05; %how many units of enzyme are produced for every unit glc consumed?
v.initglc = 0; %needed to bootstrap production
v.initmedia = 1; %amount of all other media added
v.initO2 = 100; %oxygen level in the media
v.initNH4 = 100; %ammonium level in the media
v.initPi = 100;
v.initSO4 = 100;
v.initK = 100;
v.initCO2 = 100;
v.initH2O = 100;
v.initergst = 0.3967; %396.65 g/mol, .01g/L, 100mL
v.initzymst = v.initergst;
v.initgthox = 0;
v.initfe2 = 100;
v.initenzyme = v.initialpop * v.alpha;
v.format = p.layouttype{1};
v.continuous = startsWith(v.format,'c');

mmol_glc = 1000 * 2 / (180.156 - 17); 
mmol_avicel = mmol_glc / v.glcpercellulose;
v.initcellulose = mmol_avicel * 100/450;

weightpercell = .027/(4e9); %6.75e-12
v.initialpop = 2e7 * weightpercell;

%%COMETS EXECUTION VARIABLES
v.maxcycles = 250;
v.timestep = 1; %1/(60*60); %1/3600 hours = 1 second.
v.filepath = 'C:\sync\biomes\cellulose\optima\temp';
v.writemedialog = true;
v.writebiomasslog = true;
v.writefluxlog = false;
v.lograte = 1;
v.maxspacebiomass = 1000;
v.spaceWidth = 10^(1/3);

v.costfactor = p.costfactor;

%print values
p = p;
v_default = v;

basestoich = getEnzymeStoich(v);
basestoich = basestoich(1:20);
seqlength = -sum(basestoich);


%model = buildScerevisiaeModel(v_default,model_default,'AA_model','fixed');
model = prepareYeastGEMModel('model',model_default,'alpha',v.alpha);
model.v = v;

baseLayout = buildCelOptLayout(model);
baseLayout.params.writeBiomassLog = true;
baseLayout.params.writeMediaLog = true;

baseLayout = setInitialMedia(baseLayout,'enzyme[e]',v.initenzyme);
baseLayout = setMedia(baseLayout,'zymst[e]',v.initzymst);
baseLayout = setMedia(baseLayout,'ergst[e]',v.initergst);
baseLayout = setMedia(baseLayout,'cellulose',v.initcellulose);
baseLayout = setMedia(baseLayout,'pi[e]',v.initPi);
baseLayout = setMedia(baseLayout,'so4[e]',v.initSO4);
baseLayout = setMedia(baseLayout,'nh4[e]',v.initNH4);
baseLayout = setMedia(baseLayout,'glc-D[e]',v.initglc);
baseLayout = setMedia(baseLayout,'fe2[e]',v.initfe2);
baseLayout = setMedia(baseLayout,'o2[e]',v.initO2);

enzIdx = stridx('produce_enzyme',baseLayout.models{1}.rxns);

aaIdx = zeros(1,20);
for i = 1:20
    aaIdx(i) = stridx(aminos_default{i},baseLayout.models{1}.mets,false);
end

alltabs = table();

currdir = pwd;
cd 'C:\sync\biomes\cellulose\optima\temp';
changeCobraSolver('gurobi');

for j = 1:length(foldchange)
    fc = foldchange(j);
    
    layouts = cell(21,1);
    tab = table();
    
    %increase each AA
    for i = 1:20
        layout = baseLayout;
        layout.models{1} = alterEnzymeAAs(layout.models{1},aminos_default{i},fc,aminos_default,trna_free_names);
        layouts{i} = layout;
    end
    
    %row 21 is the default
    layouts{21} = baseLayout;
    
    tab.layout = layouts;
    aminos = [aminos_default {'none'}];
    tab.AA = aminos';
    tab.foldchange(1:21) = fc;
    tab.costfactor(1:21) = p.costfactor;
    tab.max_obj_flux = zeros(21,1);
        
    for i = 1:21
        disp(['Executing layout ' num2str(i) ' of 21.'])
        l = tab.layout{i};
        runComets(l);
        biomass = parseBiomassLog(l.params.biomassLogName);
        tab.biomass{i} = biomass;
        %tab.biomass_max(i) = max(biomass);
        
        media = parseMediaLog(l.params.mediaLogName,{'cellulose','glc-D[e]','enzyme[e]'});
        cel = media(stridx('cellulose',media.metname),:);
        glc = media(stridx('glc-D[e]',media.metname),:);
        enz = media(stridx('enzyme[e]',media.metname),:);
        tab.cellulose{i} = cel.amt;
        tab.glucose{i} = glc.amt;
        tab.enzyme{i} = enz.amt;
        tab.timestep{i} = cel.t;
        tab.t{i} = cel.t * l.models{1}.v.timestep;
        opt = optimizeCbModel(l.models{1});
        tab.max_obj_flux(i) = opt.f;
    end
    
    [nr,nc] = size(alltabs)
    if nr == 0
        alltabs = tab;
    else
        alltabs = vertcat(alltabs,tab);
    end
end
cd(currdir);

tab = alltabs;