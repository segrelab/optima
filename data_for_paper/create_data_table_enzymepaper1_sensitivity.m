function dat = create_data_table_enzymepaper1_sensitivity(varargin)
%%
% Creates or updates a data table ("dat") which contains E coli models that
% have completed COMETS runs and have their data stored
%
% Data columns:
%
% * version- increment this every time you modify the underlying generation
% scripts so you know the result is valid
% * alpha
% * mode- eg 'proportional', 'constitutive', 'demand', etc
% * layouttype- '1x1', 'chemostat',
% * deathrate- for biomass
% * decayrate- for enzyme
% * layout- the layout, including the model, whose execution was called to create the result
% * costfactor- scale of enzyme production cost relative to biomass
% stoichiometry
%%
%%SCRIPT CONTROL VARIABLES
refreshmodel = true; %load the default model from file instead of using the 'model_default' in the workspace
updateOldRows = false; %If the version on an existing record is below the current, should it be regenerated?
purgeModelFields = true; %Clear out most of the info in the model structs to save space
extend = false; %if the cellulose is not completely cleared out, should the time be extended and the sim rerun?
p.version = 18;
%Version History
%1: Initial Setup
%2: Increase initial O2 and NH4 in media from 100mmol to 1000mmol
%3: Switch organism from E coli to S cerevisiae
%4: change AA composition of enzyme to match Punctularia strigosozonata
%       CBHI(1) http://www.uniprot.org/uniprot/R7S4L5
%5: change proportional model so that enzyme amount is # of molecules,
%not grams of enzyme
%6: add ethanol to the set of tracked metabolites. set O2 to be static
%7: set H2O and CO2 to be static. Increased amount of other media mets
%to 1000mmol
%8: add ATP and GTP costs to enzyme production
%9: Switch to "fixed" mode so alpha denotes mmol enzyme / g growth
%10: add pi_c as a product to enzyme production stoichiometry
%11: restore ATP maintenance reaction lower bound = 1
%12: 6/8/2018: add continuous/chemostat culture format to buildCelOptLayout
%13: 7/12/2018: add ergosterol and zymosterol exchange rxns to buildSCerevisiaeModel so that anaerobic sims can be performed
%       Note:these are already present in iMM904. When feeding these to iMM904, use the metabolites zymstest_SC_e and ergstest_SC_e
%14: change Amino Acid recipe for enzyme to be based on average enzymes
%across 3 species. See "getEnzymeScoith.m" and "s4.2_AA_comp_species.mlx"
%15: 8/30/18. Alter media names to use [e] format, switched to '_renamed' model. Removed initial glucose, replaced with initial enzyme = alpha * initialpop
%16: No major changes, just changing some default values
%17: Change default kinetics to match best fit for DenHaan data
%18: Change model to the YeastGEM v8.3 model

%%NEW MODEL PARAMETERS
% p.alpha = .1;
% p.enzdecayperhour = .2;
% p.deathrate = .1;
% p.dilutionrate = 1/12; %fraction of media replaced per hour
% p.mode = {'Proportional'};
% p.layouttype = {'1x1'};
% p.diff_enz = 1e-6; %Default is 1e-6;

%p.alpha = [.1 .175 .25 .375 .5 0.625 0.75 1 1.25 1.5 1.75 2];
%p.alpha = [.05 .15 .2 .3 .35 .4 .45 .65 .85 1.25 1.75 .6:.1:1];
%p.alpha = [.001 .0025 .005 .01 .025 .05 .1 .15 .2 .25 .3 .35 .4 .5 .6 .7 .8 .9 1];
%p.alpha = [.0001 .00025 .0005 .00075 .1:.01:.3 .99 .05:.05:.95];
%p.alpha = [.005 .01 .015 .02 .025 .03 .05 .1 .15 .2:.1:.9];
%p.alpha = [.0025 .025 .25];
%p.alpha = 0.15;

%1 millimole = v.enzweight*6.02e23 / 1000 grams
%set an alpha range from 1-30% protein weight
v.enz_weight = 8.5853e-20;
mmweight = v.enz_weight * 6.02214e23 / 1000; %g/mmol
pcts = 1:5:30;
p.alpha = zeros(1,length(pcts));
for i = 1:length(pcts)
    %pct/100 = weight of x mmol
    pct = pcts(i);
    p.alpha(i) = (pct/100) / mmweight;
end

v = getOptVars();

%p.alpha = [.1 .01 .001 .0005 .005 .05 .0075 .075 .75];
%p.alpha = [p.alpha 0.2 .01:.01:1 .15 .0075 .0025];
%p.alpha = [p.alpha .02 .03 .06 .0025 .75 .05 .01 .015 .0075 .005 .001 .2 .15 .04 .4 .75 1];
%p.alpha = [.0001 .00005 .0005 .001 .005 .01 .015 .02 .03 .04 .06 .1 .13 .175 .2];
p.alpha = [exp(-6.8) exp(-6.6) exp(-6.4) exp(-6.2) exp(-6) exp(-5.8) exp(-5.6) exp(-5.4) exp(-5.2) exp(-5) exp(-8)];
p.alpha = [exp(-8) exp(-7.5) exp(-7) exp(-6.5) exp(-6) exp(-5.5) exp(-5) exp(-4.5) exp(-4) exp(-3.5) exp(-3) exp(-2.5) exp(-2) exp(-1.5) exp(-1) exp(-.5) exp(0)];
p.alpha = [p.alpha exp(-8.5) exp(-9) exp(-9.5) exp(-10)];
p.alpha = exp(-3.25);
%p.alpha = [exp(-3.5*.9) exp(-3.5*1.1) exp(-5.25) exp(-1.75)];
%p.enzdecayperhour = [1/24 1/12 1/6 .1];
%p.enzdecayperhour = [0 .1 .2 .05 .01];
p.enzdecayperhour = [0 1/2400 1/1200 1/600 1/4800 1/12000 1/240 1/180 1/120 1/360 1/480];
%p.enzdecayperhour = [1/6 1/12 1/18 1/24 1/36 1/48];
%p.enzdecayperhour = [1/18 1/24 1/36 1/48];
%p.enzdecayperhour = [0 1/1000 1/100]
%p.enzdecayperhour = [1/500 1/200 1/100 1/1000 0];
p.enzdecayperhour = .01;
v.enzdecayperhour = p.enzdecayperhour(1);
%p.enzdecayperhour = [1/240 .9/240 1.1/240 1.5/240 0.5/240];

%p.deathrate = [0 .1 .2 .05 .01];
%p.deathrate = [0 0.01];
%p.deathrate = [0.00125 (0.00125 * .9) (0.00125 * 1.1) (0.00125 - 0.000625) (0.00125 + 0.000625)];
p.deathrate = 0.001;
p.deathrate = 0;
v.deathrate = p.deathrate(1);
%p.deathrate = [.1 .2];
p.dilutionrate = 0; %fraction of media replaced per hour
p.mode = {'Fixed'};%{'Constitutive','Proportional','Fixed};
p.layouttype = {'1x1'}; %{'1x1','chemostat'};
p.diff_enz = 1e-6; %Default is 1e-6;
%p.costfactor = [1 .5 2];
p.costfactor = 1;
v.costfactor = p.costfactor(1);
v.richmedia = false; %provide a bunch of everything

%%STRATEGY VARIABLES
v.defaultbound = 10; %default exchange upper bound
%v.maxEnzymeRate = 1; %maximum enzyme:biomass ratio by weight for fixed production. 1:1 w/ alpha=1 should mean half of all amino acids go into enzyme
v.alpha = p.alpha(1); %fixed percentage of the maximum enzyme rate being expressed

%%ENZYME VARIABLES. See https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.4
%values used: wild-type T.ressei on 4-methylumbelliferyl cellobioside. https://www.brenda-enzymes.org/literature.php?e=3.2.1.4&r=737187
v.kcat_cel = 27; %mmol converted / mmol enzyme / second  
v.km_cel = 0.03; %Units = mmoles/cm^3v.enz_weight =  8.5853183e-20; %grams
v.enzByWeight = false;
v.atp_per_peptide = 2; %used to calculate energy costs in getEnzymeStoich
v.gtp_per_peptide = 2; %see https://www.ncbi.nlm.nih.gov/books/NBK224633/

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

v.initglc = 0;
v.initmedia = 10; %amount of all other media added
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.initFe2 = 1000;
v.initergst = 0.3967; %396.65 g/mol, .01g/L, 100mL
v.initzymst = v.initergst;
v.format = p.layouttype{1};
v.continuous = startsWith(v.format,'c');

mmol_glc = 1000 * 2 / (180.156 - 17); 
mmol_avicel = mmol_glc / v.glcpercellulose;
v.initcellulose = mmol_avicel * 100/450;

weightpercell = .027/(4e9); %6.75e-12
v.initialpop = 2e7 * weightpercell;

v.initenzyme = v.initialpop * v.alpha;

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
v.lograte = 1;

%print values
p = p
v_default = v

%%DEFAULT MODEL
if refreshmodel
    %model_default = load('C:\sync\biomes\models\ecoli_core_model_norm.mat');
    %model_default = load('C:\sync\biomes\models\iJO1366 (1).mat');
    %model_default = model_default.iJO1366;
%     model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904
%     model_default = model_default.model;
    %S cerevisiae model ref: http://bigg.ucsd.edu/models/iMM904 https://www.ncbi.nlm.nih.gov/pubmed/19321003
    model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
    model_default = model_default.model;
end

%build the list of new models
models = {};
models{end+1} = build(model_default,p,v_default,p.deathrate,p.enzdecayperhour,p.alpha,v.initcellulose,v.kcat_cel,v.km_cel,v.timestep,v.maxcycles,p.costfactor);
multipliers = [2 4 8 16 32 64 1/2 1/4 1/6 1/8 1/16 1/32 1/64];
%multipliers = 1;
%multipliers = 2;
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    deathrate = deathrate * multipliers(i);
    v.deathrate = deathrate;
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    decayrate = decayrate * multipliers(i);
    v.enzdecayperhour = decayrate;
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    alpha = alpha * multipliers(i);
    v.alpha = alpha;
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    initcell = initcell * multipliers(i);
    v.initcellulose = initcell;
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    kcat = kcat * multipliers(i);
    v.kcat_cel = kcat;
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    km = km * multipliers(i);
    v.km_cel = km;
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end
% for i = 1:length(multipliers)
%     v = v_default;
%     deathrate = v_default.deathrate;
%     decayrate = v_default.enzdecayperhour;
%     alpha = v_default.alpha;
%     initcell = v_default.initcellulose;
%     kcat = v_default.kcat_cel;
%     km = v_default.km_cel;
%     ts = v_default.timestep;
%     maxcycles = v_default.maxcycles;
%     cost = v_default.costfactor;
%     
%     ts = ts * multipliers(i);
%     v.timestep = ts;
%     
%     models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
% end
% for i = 1:length(multipliers)
%     v = v_default;
%     deathrate = v_default.deathrate;
%     decayrate = v_default.enzdecayperhour;
%     alpha = v_default.alpha;
%     initcell = v_default.initcellulose;
%     kcat = v_default.kcat_cel;
%     km = v_default.km_cel;
%     ts = v_default.timestep;
%     maxcycles = v_default.maxcycles;
%     cost = v_default.costfactor;
%     
%     maxcycles = maxcycles * multipliers(i);
%     v.maxcycles = maxcycles;
%     
%     models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
% end
for i = 1:length(multipliers)
    v = v_default;
    deathrate = v_default.deathrate;
    decayrate = v_default.enzdecayperhour;
    alpha = v_default.alpha;
    initcell = v_default.initcellulose;
    kcat = v_default.kcat_cel;
    km = v_default.km_cel;
    ts = v_default.timestep;
    maxcycles = v_default.maxcycles;
    cost = v_default.costfactor;
    
    cost = cost * multipliers(i);
    v.costfactor = cost;
    
    models{end+1} = build(model_default,p,v,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost);
end


% Cell Death Rate, 0 : 0.0025 (0.00125)
% Enzyme Decay Rate, 0 : 1/120 (1/240)
% Alpha: exp(-7) : exp(0) ( exp(-3) )
% Initial Cellulose Amount: .001mmol : 2mmol (1 mmol)
% Enzyme Kcat: .001 : 10 (5)
% Enzyme KM .01 : 1 (0.5)
% Timestep Length: .05 : 2 (1)
% Total Run Time: 500 : 2000 (1250)
% Cost Factor: 0 : 2 (1)





%if given a table, we'll be updating it. Otherwise, create one
if ~isempty(varargin)
    dat = varargin{1};
    %check if each model matches something in dat, add it if not
    for i = 1:length(models)
        model = models{i};
        rowid = find(model.v.alpha == dat.alpha & ...
            strcmp(model.v.mode,dat.mode) & ...
            strcmp(model.v.layouttype,dat.layouttype) & ...
            model.v.deathrate == dat.deathrate & ...
            model.v.enzdecayperhour == dat.decayrate & ...
            model.v.dilutionrate == dat.dilutionrate & ...
            model.v.diff_enz == dat.diff_enz & ...
            model.v.costfactor == dat.costfactor);
        if any(rowid)
            for r = 1:length(rowid) %should only be one, but loop to be safe
                currowid = rowid(r);
                if dat.version(currowid) < p.version %the layout/model is out of date
                    layout = buildlayout(model,v);
                    if strcmp(model.v.layouttype,'chemostat')
                        layout = makeChemostat(layout,deathrate,1);
                    end
                    layout = setMedia(layout,'enzyme[e]',model.v.alpha * model.v.initialpop);
                    dat.layout{currowid} = layout;
                end
            end
        else %if ~any(rowid)
            layout = buildlayout(model);
            deathrate = model.v.deathrate;
            if strcmp(model.v.layouttype,'chemostat')
                layout = makeChemostat(layout,deathrate,1);
            end
            version = 0;
            layout = setMedia(layout,'enzyme[e]',model.v.alpha * model.v.initialpop);
            alpha = model.v.alpha;
            mode = {model.v.mode};
            layouttype = {model.v.layouttype};
            decayrate = model.v.enzdecayperhour;
            diff_enz = model.v.diff_enz;
            biomass_final = 0;
            dilutionrate = model.v.dilutionrate;
            costfactor = model.v.costfactor;
            layout = {layout};
            tab = table(version,alpha,mode,layouttype,deathrate,decayrate,dilutionrate,layout,diff_enz,costfactor,biomass_final);
            dat = [dat;tab];
            
        end
    end
else
    for i = 1:length(models)
        model = models{i};
        layout = buildlayout(model);
        deathrate = model.v.deathrate;
        if strcmp(model.v.layouttype,'chemostat')
            layout = makeChemostat(layout,deathrate,1);
        end
        layout = setMedia(layout,'enzyme[e]',model.v.alpha * model.v.initialpop);
        version = 0;
        alpha = model.v.alpha;
        mode = {model.v.mode};
        layouttype = {model.v.layouttype};
        decayrate = model.v.enzdecayperhour;
        diff_enz = model.v.diff_enz;
        biomass_final = 0;
        dilutionrate = model.v.dilutionrate;
        costfactor = model.v.costfactor;
        
        if ~exist('dat','var')
            %create dat
            
            diff_enz = model.v.diff_enz;
            biomass_final = 0;
            
            model = {model};
            layout = {layout};
            dat = table(version,alpha,mode,layouttype,deathrate,decayrate,dilutionrate,layout,diff_enz,costfactor,biomass_final);
            %if purgeModelFields
            %    dat = purgeModelFieldsInDat(dat);
            %end
        else
            model = {model};
            layout = {layout};
            tab = table(version,alpha,mode,layouttype,deathrate,decayrate,dilutionrate,layout,diff_enz,costfactor,biomass_final);
            %if purgeModelFields
            %    tab = purgeModelFieldsInDat(tab);
            %end
            dat = [dat;tab];
        end
    end
end

%check for old versions of records, and rebuild the layout
[nrows, ncols] = size(dat);
if updateOldRows
    for i = 1:nrows
        row = dat(i,:);
        if row.version < p.version
            v = v_default;
            
            v.alpha = row.alpha;
            v.mode = row.mode;
            v.layouttype = row.layouttype;
            v.diff_enz = row.diff_enz;
            v.costfactor = row.costfactor;
            v.deathrate = row.deathrate;
            v.enzdecayperhour = row.decayrate;
            desc = [v.mode '-A' num2str(v.alpha) '-DR' num2str(v.deathrate) '-ED' num2str(v.enzdecayperhour) '-C' num2str(v.costfactor)];
            if strcmpi('chemostat',row.layouttype)
                v.enzdecayperhour = 0;
                v.deathrate = 0;
                v.dilutionrate = row.dilutionrate;
                desc = [v.mode '-A' num2str(v.alpha) '-D' num2str(v.dilutionrate) 'C' num2str(v.costfactor)];
            end
            newmodel = buildScerevisiaeModel(v,model_default,desc,v.mode);
        else
            v.enzdecayperhour = row.decayrate;
            v.deathrate = row.deathrate;
            v.dilutionrate = 0;
            desc = [v.mode '-A' num2str(v.alpha) '-DR' num2str(v.deathrate) '-ED' num2str(v.enzdecayperhour) '-C' num2str(v.costfactor)];
            newmodel = buildScerevisiaeModel(v,model_default,desc,v.mode);
        end
        v.alpha = row.alpha;
        newlayout = buildlayout(newmodel);
        dat.layout{i} = newlayout;
    end
end

%execute and replace the layout for every model whose version is lower than
%the current version
curdir = pwd;
tdir = 'C:\sync\biomes\cellulose\optima\temp';
cd(tdir);
[nrows, ~] = size(dat);
for i = 1:nrows
    disp(['Evaluating record ' num2str(i) ' of ' num2str(nrows)]);
    if dat.version(i) < p.version
        runlayout = executeDatLayout(dat.layout{i},extend);
        dat.layout(i) = {runlayout};
        dat.version(i) = p.version;
        dat.biomass_final(i) = runlayout.models{1}.biomass(end);
    end
end
cd(curdir);

if purgeModelFields
    dat = purgeModelFieldsInDat(dat);
end
end
%%
function model = build(model_default,p,v_default,deathrate,decayrate,alpha,initcell,kcat,km,ts,maxcycles,cost)


v = v_default;
v.alpha = alpha;
v.mode = p.mode{1};
v.layouttype = p.layouttype{1};
v.diff_enz = p.diff_enz(1);
v.costfactor = cost;
v.enzdecayperhour = decayrate;
v.deathrate = deathrate;
v.dilutionrate = 0;
v.kcat_cel = kcat;
v.km_cel = km;
v.maxcycles = maxcycles;
v.timestep = ts;
v.initinitcellulose = initcell;
desc = ['YeastGEM-' v.mode '-A' num2str(v.alpha) '-DR' num2str(v.deathrate) '-ED' num2str(v.enzdecayperhour) '-C' num2str(v.costfactor)];
%model = buildScerevisiaeModel(v,model_default,desc,v.mode);
model = prepareYeastGEMModel('model',model_default,'alpha',v.alpha);
model.v = v;
model.description = desc;
model.enzExpressionMode = 'proportional';
end

function layout = executeDatLayout(layout,extend)
runComets(layout);

tab = parseBiomassLog('log_biomass.m');
layout.models{1}.biomass = tab.biomass; %no need to sum if there's no biomass diffusion
layout.models{1}.timestep = tab.t;

if layout.models{1}.v.writemedialog
    media = parseMediaLog([layout.models{1}.v.filepath '\log_media.m'],{'enzyme[e]', 'cellulose', 'glc-D[e]', 'ac[e]', 'etoh[e]'});
    layout.models{1}.enzyme_amt = zeros(length(layout.models{1}.timestep),1);
    layout.models{1}.cellulose_amt = zeros(length(layout.models{1}.timestep),1);
    layout.models{1}.glc_amt = zeros(length(layout.models{1}.timestep),1);
    layout.models{1}.acetate_amt = zeros(length(layout.models{1}.timestep),1);
    for i = 1:length(layout.models{1}.timestep)
        time = layout.models{1}.timestep(i);
        enzrows = (media.t == time) & strcmp(media.metname,'enzyme[e]');
        celrows = (media.t == time) & strcmp(media.metname,'cellulose');
        glcrows = (media.t == time) & strcmp(media.metname,'glc-D[e]');
        acrows = (media.t == time) & strcmp(media.metname,'ac[e]');
        erows = (media.t == time) & strcmp(media.metname,'etoh[e]');
        subtab = media(enzrows,:);
        layout.models{1}.enzyme_amt(i) = sum(subtab.amt);
        subtab = media(celrows,:);
        layout.models{1}.cellulose_amt(i) = sum(subtab.amt);
        subtab = media(glcrows,:);
        layout.models{1}.glc_amt(i) = sum(subtab.amt);
        subtab = media(acrows,:);
        layout.models{1}.acetate_amt(i) = sum(subtab.amt);
        subtab = media(erows,:);
        layout.models{1}.etoh_amt(i) = sum(subtab.amt);
    end
    
    %if the model didn't digest all the carbon, extend the time and run
    %it again
    if (layout.models{1}.cellulose_amt(end) + layout.models{1}.glc_amt(end)) > 1e-8
        if (extend && (layout.models{1}.v.timestep * layout.models{1}.v.maxcycles) < 13140) %1/2 year = 4380
            disp('Model failed to complete cellulose digestion. Extending time and rerunning...');
            layout.params.timeStep = layout.params.timeStep * 2;
            layout.params.maxCycles = layout.params.maxCycles * 1.5;
            layout.models{1}.v.timestep = layout.models{1}.v.timestep * 2;
            layout.models{1}.v.maxcycles = layout.models{1}.v.maxcycles * 1.5;
            layout = executeDatLayout(layout,extend);
        end
    end
end
end

function newlayout = buildlayout(model)
newlayout = buildCelOptLayout(model);
newlayout = setMedia(newlayout,'cellulose',model.v.initcellulose);
newlayout = setMedia(newlayout,'zymst[e]',model.v.initzymst);
newlayout = setMedia(newlayout,'ergst[e]',model.v.initergst);
newlayout = setMedia(newlayout,'zymstest_SC[e]',model.v.initzymst);
newlayout = setMedia(newlayout,'ergstest_SC[e]',model.v.initergst);
newlayout = setMedia(newlayout,'enzyme[e]',model.v.alpha * model.v.initialpop);
newlayout = setMedia(newlayout,'o2[e]',model.v.initO2);
newlayout = setMedia(newlayout,'pi[e]',model.v.initPi);
newlayout = setMedia(newlayout,'so4[e]',model.v.initSO4);
newlayout = setMedia(newlayout,'nh4[e]',model.v.initNH4);
newlayout = setMedia(newlayout,'fe2[e]',model.v.initFe2);
% newlayout = setMedia(newlayout,'ocdca[e]',model.v.initergst * 42);
% newlayout = setMedia(newlayout,'ocdcya[e]',model.v.initergst * 42);
% newlayout = setMedia(newlayout,'ocdcea[e]',model.v.initergst * 42);
% newlayout = setMedia(newlayout,'hdcea[e]',model.v.initergst * 42);
% newlayout = setMedia(newlayout,'hdca[e]',model.v.initergst * 42);
% newlayout = setMedia(newlayout,'ttdca[e]',model.v.initergst * 42);


%newlayout = setMedia(newlayout,'glc-D[e]',model.v.initergst * 42);

%newlayout.models{1} = addReaction(newlayout.models{1},'EX_triglyc_SC[c]',{'triglyc_SC[c]'},-1,'lowerBound',-99999);
dep_amt = 10;
%newlayout = setMedia(newlayout,'triglyc_SC[c]',dep_amt);
end