%Replicate the conditions in the COMETS paper and observe results when the
%media has a small amount of insoluble cellulose

%Claim 1: colonies increase linearly in diameter over time

%Width is calculated by counting the number of cells that have more than
%10^-7g of biomass across a horizontal slice through the center of the
%colony

v.initialpop = 3e-7;

%initenzymeamts = v.initialpop * exp([-3 -5 -7]);
initenzymeamts = .01;
alphas = exp([-7 -6 -5 -4 -3 -2]);
maxcycles = 500;
includeProportionalInitEnzyme = true;

%temp shorter settings to confirm functionality
% initenzymeamts = [.001 .0001];
% alphas = exp(-4);
% maxcycles = 20;
% includeProportionalInitEnzyme = false;

% Layout conditions: (Cite COMETS paper 2014)
% uptake vmax 10mmol/g/h
% uptake km 10uM
% death rate 1%
% metabolite diffusion 5x10^-6 cm^2/s
% biomass diffusion 3x10^-10 cm2/s for colony growth sims, 3x10^-9 for
% others
% max colony height 200um
% oxygen concentration 250umol/cm^2
% space width 0.02 or 0.05 cm
% carbon concentration 0.0088 mmol glc/cm^2
% KM and kcat set according to Section 2 results
% time step 0.1 h
% initial biomass 3x10^-7 g
% plate size 50x50
v = getOptVars();

v.km_cel = .03;
v.kcat_cel = 27;

weightpercell = .027/(4e9); %6.75e-12
v.initialpop = 2e7 * weightpercell;
v.km_glc = v.km_glc / 10;
v.vmax_glc = v.vmax_glc * 1;
v.glcuptakerate = v.vmax_glc;
v.diff_glc = 5e-6;
v.deathrate = .01;
v.decayrate = 0; %?change?
v.initO2 = 0;
v.staticO2 = 0.25;
v.timestep = .1;
v.spacewidth =0.02;
v.xdim = 50;
v.ydim = 50;
v.diff_cel = 0;
v.diff_enz = 5e-6;
v.diff_glc = 5e-6;

v.alpha = exp(-4);

v.initenzyme = v.initialpop * v.alpha;

%calculate the initial cellulose
glc_con = 0.0088 * (v.spacewidth^2);
v.initcellulose = glc_con / v.glcpercellulose;

% %temp values to make sure it can grow at all
% % v.initenzyme = 3;
% v.initcellulose = 10;
%v.initcellulose = v.initcellulose * 10000000;

loadmodel = true;
%loadmodel = false;
if loadmodel
%     load('C:\sync\biomes\models\ecoli\ecoli_normalized.mat'); %'ecoli'
%     model = ecoli;
%     modelname = 'enzyme_ecoli';
%     model = load('C:\sync\biomes\models\Scerevisiae_iMM904_renamed.mat');
%     model = model.model;
%     modelname = 'iMM904_Scerevisiae_enzyme';
    model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
    model_default = model_default.model;
end
obj = model_default.rxns{find(model_default.c)};

nrows = length(alphas) * length(initenzymeamts);
firstInitGroup = 1;
if includeProportionalInitEnzyme
    firstInitGroup = 0;
    nrows = nrows + length(alphas);
end
biodat = table();
biodat.biomass = cell(nrows,1);
biodat.initenzyme = zeros(nrows,1);
biodat.initenzymegroup = zeros(nrows,1);
biodat.alpha = zeros(nrows,1);
biodat.logalpha = zeros(nrows,1);
currow = 1;
currdir = pwd;
cd 'C:\sync\biomes\cellulose\optima\temp';

for ieg = firstInitGroup:length(initenzymeamts) %ieg = init enzyme group
    for a = alphas
        v.alpha = a;
        
        if ieg == 0
            v.initenzyme = v.initialpop * v.alpha;
        else
            v.initenzyme = initenzymeamts(ieg);
        end
        disp(['Building Layout ' num2str(currow) ' of ' num2str(nrows)]);
        %model = buildSpeciesModel(v,model,modelname,'fixed',obj);
        model = prepareYeastGEMModel('model',model_default,'alpha',v.alpha);
        model.v = v;
        
        layout = buildCelOptLayout(model);
        layout.params.spaceWidth = v.spacewidth;
        layout.params.defaultDiffConst = 5e-6;
        layout.params.growthDiffRate = 3e-10;
        layout.params.flowDiffRate = 3e-10;
        layout = setGlobalStaticMedia(layout,'o2[e]',v.staticO2);
        layout.media_amt(stridx('cellulose',layout.mets)) = v.initcellulose;
        layout = setInitialMediaInCell(layout,25,25,'enzyme[e]',v.initenzyme);
        layout.params.writeMediaLog = false;
        layout.params.maxCycles = maxcycles;
%         
%         %TODO: remove
%         layout = setDiffusion(layout,'enzyme[e]',0);
%         layout = setDiffusion(layout,'glc-D[e]',0);
        
        disp('Executing Comets');
        runComets(layout);
        
        disp('Loading Biomass Log');
        biomass = parseBiomassLog(layout.params.biomassLogName);
        
        disp('Adding data to the biodat table'); 
        biodat.biomass{currow} = biomass;
        biodat.initenzyme(currow) = v.initenzyme;
        biodat.initenzymegroup(currow) = ieg;
        biodat.alpha(currow) = v.alpha;
        biodat.logalpha(currow) = log(v.alpha);
        
        currow = currow + 1;
    end
end
cd(currdir);