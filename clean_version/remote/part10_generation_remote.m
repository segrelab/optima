function layouts = part10_generation_remote(arg)
%PART6_GENERATION_REMOTE
%Creates a set of 2d layouts with variable enzyme diffusivity 

xdim = arg.v.xdim;
ydim = arg.v.ydim;
%distances = [0,2];
%distances = 5;

diffscales = arg.diffscales;
enzamts = arg.enzamts;

v = arg.v;
model_default = load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
model = model_default.model;
enzmodel = prepareYeastGEMModel('model',model,'alpha',v.alpha,'vmax_glc',v.vmax_glc,'maintRxnAsObjective',v.maintRxnAsObjective);
enzmodel.v = v;

layout = createLayout();
layout = setDims(layout,xdim,ydim);

drainAtEdges = false;%should enzyme and glc be removed at edges (as opposed to being allowed to bounce back)
if isfield(v,'drainatedges')
    drainAtEdges = v.drainatedges;
end
if drainAtEdges
    xdim = xdim + 2;
    ydim = ydim + 2;
    v.xdim = xdim;
    v.ydim = ydim;
    enzmodel.v.xdim = xdim;
    enzmodel.v.ydim = ydim;
    layout = setDims(layout,xdim,ydim);
end

%make the plate circular
layout = addCircularBarrier(layout);

if drainAtEdges
    mask = layout.barrier;
    mask(1:xdim,1) = 1;
    mask(1:xdim,ydim) = 1;
    mask(1,1:ydim) = 1;
    mask(xdim,1:ydim) = 1;
    
    for x = 2:xdim-1
        for y = 2:ydim-1
            edgeneighbors = sum(layout.barrier(sub2ind(size(layout.barrier),[(x-1) x x (x+1)],[y (y+1) (y-1) y])));
            if edgeneighbors > 0
                mask(x,y) = 1;
            end
        end
    end
    
    layout = applyStaticMediaMask(layout,mask,{'enzyme[e]','glc-D[e]','cellulose'},[0 0 0]);
    %layout.barrier = zeros(size(layout.barrier)); %allow diffusion into this space, don't bounce back
     % %disabled so other nutrients don't get lost in here
end
layout = addModel(layout,enzmodel);

enzposn = [ceil(xdim/2),ceil(xdim/2)]; %location of the enzyme producer
layout.initial_pop(1,enzposn(1),enzposn(2)) = v.initialpop;

layout = setMedia(layout,layout.mets,0);

layout = setMedia(layout,'cellulose',v.initcellulose);
layout = setInitialMediaInCell(layout,enzposn(1),enzposn(2),'enzyme[e]',v.initenzyme);
layout = setMedia(layout,'zymst[e]',v.initzymst);
layout = setMedia(layout,'ergst[e]',v.initergst);
% layout = setMedia(layout,'zymstest_SC[e]',v.initzymst);
% layout = setMedia(layout,'ergstest_SC[e]',v.initergst);
layout = setMedia(layout,'o2[e]',v.initO2);
layout = setMedia(layout,'h[e]',v.initH);
layout = setMedia(layout,'h2o[e]',v.initH2O);
layout = setMedia(layout,'co2[e]',v.initCO2);
layout = setMedia(layout,'pi[e]',v.initPi);
layout = setMedia(layout,'so4[e]',v.initSO4);
layout = setMedia(layout,'nh4[e]',v.initNH4);
layout = setMedia(layout,'k[e]',v.initK);
layout = setMedia(layout,'fe2[e]',v.initfe2);
layout = setMedia(layout,'glc-D[e]',v.initglc);
layout = setMedia(layout,'gthox[e]',v.initgthox);
layout = setMedia(layout,'oleate[e]',v.initoleate);
layout = setMedia(layout,'palmitoleate[e]',v.initpalmitoleate);

layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.maxSpaceBiomass = v.maxspacebiomass;
layout.params.writeBiomassLog = v.writebiomasslog;
layout.params.writeMediaLog = v.writemedialog;
layout.params.writeFluxLog = v.writefluxlog;
layout.params.biomassLogRate = v.lograte;
layout.params.mediaLogRate = v.lograte;
layout.params.fluxLogRate = v.lograte;
layout.params.biomassLogName = [v.filepath '\biomass.m'];
layout.params.mediaLogName = [v.filepath '\media.m'];
layout.params.fluxLogName = [v.filepath '\flux.m'];
layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.pixelScale = 10;
layout.params.spaceWidth = v.spaceWidth;
layout.params.deathRate = v.deathrate;


layout = addExternalReaction(layout,'degrade_cellulose',{'cellulose','glc-D[e]'},[-1 v.glcpercellulose],'enzyme', 'enzyme[e]', 'k', v.kcat_cel/v.glcpercellulose, 'km', v.km_cel);

if (isfield(v,'enzdecayperhour') && ~isfield(v,'enzdecaypersec'))
    v.enzdecaypersec = v.enzdecayperhour/3600;
end
%create the enzyme decay reaction
if isfield(v,'enzdecaypersec')
    if v.enzdecaypersec ~= 0
        layout = addExternalReaction(layout,'enzyme_decay','enzyme[e]',-1,'vmax',v.enzdecaypersec);
    end
end
if isfield(v,'numExRxnSubsteps')
    layout.params.numExRxnSubsteps = v.numExRxnSubsteps;
end

%modify diffusion rates
layout = setDiffusion(layout,'cellulose',v.diff_cel);
layout = setDiffusion(layout,'enzyme[e]',v.diff_enz);
layout = setDiffusion(layout,'glc-D[e]',v.diff_glc);

layouts = cell(length(diffscales) * length(enzamts),1);
idx = 0;
for i = 1:length(diffscales)
    for j = 1:length(enzamts)
        idx = idx + 1;
        layouts{idx} = setDiffusion(layout,'enzyme[e]',v.diff_enz * diffscales(i));
        layouts{idx}.models{1}.v.diff_enz = v.diff_enz * diffscales(i);
        layouts{idx}.models{1}.v.initenzyme = enzamts(j);
        layouts{idx} = setInitialMediaInCell(layouts{idx},ceil(v.xdim/2),ceil(v.ydim/2),'enzyme[e]',enzamts(j));
       
    end   
end
end


