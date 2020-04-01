function dat = create_data_table_enzymepaper1(varargin)
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
p.version = 15;
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

%p.alpha = [.1 .01 .001 .0005 .005 .05 .0075 .075 .75];
%p.alpha = [p.alpha 0.2 .01:.01:1 .15 .0075 .0025];
%p.alpha = [p.alpha .02 .03 .06 .0025 .75 .05 .01 .015 .0075 .005 .001 .2 .15 .04 .4 .75 1];
%p.alpha = [.0001 .00005 .0005 .001 .005 .01 .015 .02 .03 .04 .06 .1 .13 .175 .2];
p.alpha = exp([-3 -5 -7]);
%p.alpha = [exp(-3.5*.9) exp(-3.5*1.1) exp(-5.25) exp(-1.75)];
%p.enzdecayperhour = [1/24 1/12 1/6 .1]; 
%p.enzdecayperhour = [0 .1 .2 .05 .01];
p.enzdecayperhour = [0 1/2400 1/1200 1/600 1/4800 1/12000 1/240 1/180 1/120 1/360 1/480];
%p.enzdecayperhour = [1/6 1/12 1/18 1/24 1/36 1/48];
%p.enzdecayperhour = [1/18 1/24 1/36 1/48];
%p.enzdecayperhour = [0 1/1000 1/100]
%p.enzdecayperhour = [1/500 1/200 1/100 1/1000 0];
p.enzdecayperhour = 1/120;
%p.enzdecayperhour = [1/240 .9/240 1.1/240 1.5/240 0.5/240];

%p.deathrate = [0 .1 .2 .05 .01];
%p.deathrate = [0 0.01];
%p.deathrate = [0.00125 (0.00125 * .9) (0.00125 * 1.1) (0.00125 - 0.000625) (0.00125 + 0.000625)];
p.deathrate = 0.00125;
%p.deathrate = [.1 .2];
p.dilutionrate = 0; %fraction of media replaced per hour
p.mode = {'Fixed'};%{'Constitutive','Proportional','Fixed};
p.layouttype = {'1x1'}; %{'1x1','chemostat'};
p.diff_enz = 1e-6; %Default is 1e-6;
%p.costfactor = [1 .5 2];
p.costfactor = 1;
v.richmedia = false; %provide a bunch of everything

%%STRATEGY VARIABLES
v.defaultbound = 10; %default exchange upper bound
%v.maxEnzymeRate = 1; %maximum enzyme:biomass ratio by weight for fixed production. 1:1 w/ alpha=1 should mean half of all amino acids go into enzyme
v.alpha = 0.1; %fixed percentage of the maximum enzyme rate being expressed

%%ENZYME VARIABLES. See https://www.brenda-enzymes.org/enzyme.php?ecno=3.2.1.4
%values used: wild-type T.ressei on 4-methylumbelliferyl cellobioside. https://www.brenda-enzymes.org/literature.php?e=3.2.1.4&r=737187
v.kcat_cel = 5;%0.01; %mmol converted / mmol enzyme / second  %TODO: Confirm units
v.km_cel = 0.5;%0.125; %Units = mmoles/cm^3
v.enz_weight =  8.5853183e-20; %grams
v.enzByWeight = false;
v.atp_per_peptide = 2; %used to calculate energy costs in getEnzymeStoich
v.gtp_per_peptide = 2; %see https://www.ncbi.nlm.nih.gov/books/NBK224633/

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

%%WORLD/LAYOUT VARIABLES
v.xdim = 1;
v.ydim = 1;
v.diff_cel = 0; %Diffusion constants
%v.diff_enz = 1e-6; %Default is 1e-6;
v.diff_glc = 1e-6;
%v.enzdecayrate = 0; %degradation of enzyme over time
%v.deathrate = 0;
%v.enzymeperglc = 0.05; %how many units of enzyme are produced for every unit glc consumed?
v.initcellulose = 1;
v.initglc = 0; %needed to bootstrap production
v.initmedia = 1; %amount of all other media added
v.initialpop = 1e-3;
v.initO2 = 1000; %oxygen level in the media
v.initNH4 = 1000; %ammonium level in the media
v.initPi = 1000;    
v.initSO4 = 1000;
v.initK = 1000;
v.initCO2 = 1000;
v.initH2O = 1000;
v.format = p.layouttype{1};
v.continuous = startsWith(v.format,'c');

%%COMETS EXECUTION VARIABLES
v.maxcycles = 2000;
v.timestep = .1; %1/(60*60); %1/3600 hours = 1 second.
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
    model_default = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat'); %Saccaromyces cerevisiae iMM904 
    model_default = model_default.model;
    %S cerevisiae model ref: http://bigg.ucsd.edu/models/iMM904 https://www.ncbi.nlm.nih.gov/pubmed/19321003 
end

%build the list of new models
models = {};
for a = 1:length(p.alpha)
    for b = 1:length(p.enzdecayperhour)
        for c = 1:length(p.deathrate)
            for e = 1:length(p.mode)
                for f = 1:length(p.layouttype)
                    for g = 1:length(p.diff_enz)
                        for h = 1:length(p.costfactor)
                            v = v_default;
                            v.alpha = p.alpha(a);
                            v.mode = p.mode{e};
                            v.layouttype = p.layouttype{f};
                            v.diff_enz = p.diff_enz(g);
                            v.costfactor = p.costfactor(h);
                            if startsWith(p.layouttype{f},'c')
                                for d = 1:length(p.dilutionrate)
                                    v.enzdecayperhour = p.enzdecayperhour(b);
                                    v.deathrate = p.deathrate(c);
                                    v.dilutionrate = p.dilutionrate(d);
                                    desc = [v.mode '-A' num2str(v.alpha) '-D' num2str(v.dilutionrate) 'C' num2str(v.costfactor)];
                                    models{end+1} = buildScerevisiaeModel(v,model_default,desc,v.mode);
                                end
                            else
                                v.enzdecayperhour = p.enzdecayperhour(b);
                                v.deathrate = p.deathrate(c);
                                v.dilutionrate = 0;
                                desc = [v.mode '-A' num2str(v.alpha) '-DR' num2str(v.deathrate) '-ED' num2str(v.enzdecayperhour) '-C' num2str(v.costfactor)];
                                models{end+1} = buildScerevisiaeModel(v,model_default,desc,v.mode);
                            end
                            
                        end
                    end
                    
                end
            end
        end
    end
end

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
                    layout = buildCelOptLayout(model);
                    if strcmp(model.v.layouttype,'chemostat')
                        layout = makeChemostat(layout,deathrate,1);
                    end
                    layout = setMedia(layout,'enzyme[e]',model.v.alpha * model.v.initialpop);
                    dat.layout{currowid} = layout;
                end
            end
        else %if ~any(rowid)
            layout = buildCelOptLayout(model);
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
        layout = buildCelOptLayout(model);
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
        newlayout = buildCelOptLayout(newmodel);
        newlayout = setMedia(newlayout,'enzyme[e]',row.alpha * v.initialpop);
        dat.layout{i} = newlayout;       
    end
end

%execute and replace the layout for every model whose version is lower than
%the current version
curdir = pwd;
tdir = 'C:\sync\biomes\cellulose\optima\temp';
cd(tdir);
[nrows, ncols] = size(dat);
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
