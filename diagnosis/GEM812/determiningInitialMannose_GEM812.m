%See also: mannose_initialization.m

%A script to figure out how much mannose should be added to the initial
%media for replicating the DenHaan experiment

%Created March 4 2020

%trying to determine an amount that is large enough to start cellulolysis,
%but small enough to be quickly depleted

%minimum mannose:biomass for growth is 100mmol:1g

%arbitrarily decided requirements:
t_mannose_depletion = 6; %max hours until mannose is almost entirely gone
mannose_depletion_threshold = 0.01; %how much remains when called depleted
t_lagphase = 48; %max hours until exponential growth should occur
lag_threshold = 0.02; %growth rate that we consider having moved out of stationary
yield_threshold = 2; %yield that must be achieved relative to without cellulose in order to consider growth a success

load('C:\sync\biomes\cellulose\optima\diagnosis\mannoselayout.mat');
model = layout.models{1};
anaerobicModel;
model.modifications = [model.modifications {'Convert to anaerobic mode with ''yeastGEM/ComplimentaryScripts/modelCuration/anaerobicModel'' script'}];
%as 'layout'

manidx = stridx('mannose exchange',layout.models{1}.rxnNames);
mannose_weight = 180.156; %g/mol

scales = [200 150 110 100 90 50 20 10 1 .5 .1 .05 .01 .005 .001];

t = table();
v = args.p2.v;
v.maxcycles = t_lagphase * 4;%1.5;
v.vmax_man = 1e10;
v.km_man = 1e-15;
layout.models{1}.vmax(manidx) = v.vmax_man;
layout.models{1}.km(manidx) = v.km_man;
layout.params.maxCycles = v.maxcycles;
layout_default = layout;
for i = 1:length(scales)
    temp = table;
    layout = layout_default;
    v.mannose_pct = scales(i);
    v.initman = v.initialpop * v.mannose_pct * 1000 / mannose_weight;
    layout = setMedia(layout,'man-D[e]',v.initman);
    model = layout.models{1};
    runComets(layout);
    bio = parseBiomassLog(layout.params.biomassLogName);
    media = parseMediaLog(layout.params.mediaLogName);
    temp.scale(1) = scales(i);
    temp.yield_g(1) = max(bio.biomass) - min(bio.biomass);
    
    
    layout = setMedia(layout,'cellulose',0);
    layout.params.writeMediaLog = false;
    runComets(layout);
    bio_poor = parseBiomassLog(layout.params.biomassLogName);
    temp.yield_g_fromManOnly(1) = max(bio_poor.biomass) - min(bio_poor.biomass);
    y = bio.biomass - min(bio.biomass);
    y = y / temp.yield_g_fromManOnly(1);
    lag_yield = bio.t(y >= yield_threshold);
    temp.lag_yield(1) = nan;
    if any(lag_yield)
        temp.lag_yield(1) = lag_yield(1);
    end
    
    %lag phase length
    temp.biodelta{1} = bio.biomass(2:end) - bio.biomass(1:end-1);
    temp.biodelta_relative{1} = temp.biodelta{1} ./ bio.biomass(1:end-1);
    lag = bio.t(temp.biodelta_relative{1} >= lag_threshold);
    temp.lag_rate(1) = nan;
    if any(lag)
        temp.lag_rate(1) = lag(1);
    end
    %--
    
    %time to mannose depletion
    media_man = media(strcmp('man-D[e]',media.metname),:);
    man_amt = media_man.amt;
    man_init = man_amt(1);
    man_relative = man_amt / man_init;
    depletion = media_man.t(man_relative <= mannose_depletion_threshold);
    temp.depletion(1) = nan;
    if any(depletion)
        temp.depletion(1) = depletion(1);
    end
    
    temp.bio{1} = bio;
    temp.media{1} = media;
    
    temp.pass_lag(1) = ~isnan(temp.lag_rate(1)) && temp.lag_rate(1) <= t_lagphase;
    temp.pass_depletion(1) = ~isnan(temp.depletion(1)) && temp.depletion(1) <= t_mannose_depletion;
    temp.pass_yield(1) = ~isnan(temp.lag_yield(1)) && temp.lag_yield(1) <= t_lagphase;    
    if ~any(size(t))
        t = temp;
    else
        t = vertcat(t,temp);
    end
end
mannosetable = t;
