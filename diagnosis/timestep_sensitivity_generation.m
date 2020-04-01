%Generates data for the live script timestep_sensitivity_glc.mlx
%Check if a homogenous layout is sensitive to differences in timestep.
%This is a simple layout, with unmodified S cerevisiae growing in a medium with 10mmol glucose.

timesteps = [10 1 .1 .05];
%timesteps = [10 1];
hours = 200;

mode = 1; %glucose in media
%mode = 2; %extracellular reactions for cellulolysis and enzyme decay

model = load('C:\sync\biomes\models\Scerevisiae_iMM904.mat');
model = model.model;

mets = {'nh4_e', 'o2_e', 'pi_e', 'so4_e'};
    
layout = createLayout(model);
layout = setMedia(layout,mets,1000);
layout = setInitialPop(layout,'1x1');
layout.params.writeBiomassLog = true;


if mode == 1 %glc
    layout = setMedia(layout,'glc__D_e',10);
else %cellulolysis
    %layout = setMedia(layout,'glc__D_e',0.001); %needs some movement in the first timestep so the cell doesn't get deemed dead by COMETS and removed
    layout = setMedia(layout,'glc__D_e',0);
    layout = addExternalReaction(layout,'degrade_cellulose',{'cellulose_e','glc__D_e'},[-1 200],'enzyme', 'enzyme_e', 'k', 32/200, 'km', 0.07);
    layout = addExternalReaction(layout,'enzyme_decay','enzyme_e',-1,'vmax',0.1/3600);
    layout = setMedia(layout,'cellulose_e',10/200);
    layout = setMedia(layout,'enzyme_e',1);
end

layouts_ts = cell(length(timesteps),1); 
for i = 1:length(timesteps)
    layout.params.timeStep = timesteps(i);
    layout.params.maxCycles = ceil(hours / timesteps(i));
    layouts_ts{i} = layout;
end

biodat = cell(length(timesteps),1);
curdir = pwd;
cd('C:\sync\biomes\cellulose\optima\temp');
for i = 1:length(layouts_ts)
    runComets(layouts_ts{i});
    bio = parseBiomassLog(layout.params.biomassLogName);
    biodat{i} = bio.biomass;
end

cd(curdir);