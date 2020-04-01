function layout = buildCelOptLayout(model)
%Create a layout containing the given model, including the standard enzyme reaction
v = model.v;
layout = createLayout();
layout = addModel(layout,model);


% if isfield(v,'enzByWeight')
%     if v.enzByWeight
%         v.kcat_cel = v.kcat_cel / v.enz_weight;
%     end
% end
layout = addExternalReaction(layout,'degrade_cellulose',{'cellulose','glc-D[e]'},[-1 v.glcpercellulose],'enzyme', 'enzyme[e]', 'k', v.kcat_cel/v.glcpercellulose, 'km', v.km_cel);


layout = setDims(layout,v.xdim,v.ydim);
if ~isfield(v,'initialpop')
    v.initialpop = 1e-6;
end
layout = setInitialPop(layout,'colonies',v.initialpop,false);

media = {'ca2[e]', 'cbl1[e]', 'cl[e]', 'co2[e]', 'cobalt2[e]', 'cu2[e]', 'fe2[e]', ...
    'fe3[e]', 'h[e]', 'h2o[e]', 'k[e]', 'man-D[e]', 'mg2[e]', 'mn2[e]', 'mobd[e]', 'na1[e]', ...
    'nh4[e]', 'ni2[e]', 'o2[e]', 'pi[e]', 'sel[e]', 'slnt[e]', 'so4[e]', 'tungs[e]', ...
    'zn2[e]', 'zymst[e]', 'ergst[e]', 'oleate[e]', 'palmitoleate[e]'};

if isfield(v,'initmedia')
    initmedia = v.initmedia;
else
    initmedia = 1;
end

circularlayout = false;
if isfield(v,'circularlayout')
    circularlayout = v.circularlayout;
end
if circularlayout
    layout = addCircularBarrier(layout);
end

for i = 1 : length(media)
    if any(stridx(media{i},layout.mets,false))
        layout = setInitialMedia(layout,media{i},initmedia);
    end
end
layout = setMedia(layout,'cellulose',v.initcellulose);
if ~isfield(v,'initglc') && isfield(v,'initGlc')
    v.initglc = v.initGlc;
end

%set the intial enzyme and glucose to jump-start growth.
%only use the spatial intial_media if the layout is larger than 1x1
if max([v.xdim v.ydim]) > 1
    layout = setInitialMediaInCell(layout,ceil(v.xdim/2),ceil(v.ydim/2),'glc-D[e]',v.initglc);
    layout = setInitialMediaInCell(layout,ceil(v.xdim/2),ceil(v.ydim/2),'enzyme[e]',v.initenzyme);
else
    layout = setMedia(layout,'glc-D[e]',v.initglc);
    layout = setMedia(layout,'enzyme[e]',v.initenzyme);
end
    
layout = setDiffusion(layout,'cellulose', v.diff_cel);
layout = setDiffusion(layout,'enzyme[e]', v.diff_enz);
layout = setDiffusion(layout,'glc-D[e]', v.diff_glc);
if isfield(v,'diff_tre')
    layout = setDiffusion(layout,'tre[e]',v.diff_tre);
end
if isfield(v,'diff_gly')
    layout = setDiffusion(layout,'glycogen[c]',v.diff_gly);
end
if isfield(v,'diff_man')
    layout = setDiffusion(layout,'man-D[e]',v.diff_man);
end

if (~exist(v.filepath,'dir'))
    mkdir(v.filepath);
end

%apply death rate
if (isfield(v,'deathrateperhour') && ~isfield(v,'deathrate'))
    v.deathrate = v.deathrateperhour * v.timestep;
elseif (isfield(v,'deathrateperh') && ~isfield(v,'deathrate'))
    v.deathrate = v.deathrateperhour * v.timestep;
end

if ~isfield(v,'deathrate')
    v.deathrate = 0; %per timestep
end

% if (isfield(v,'decayrate') && ~isfield(v,'enzdecayperhour'))
%     v.enzdecayperhour = v.decayrate;
% end

if (isfield(v,'enzdecayperhour') && ~isfield(v,'enzdecaypersec'))
    v.enzdecaypersec = v.enzdecayperhour/3600;
end

%create the enzyme decay reaction
if isfield(v,'enzdecaypersec')
    if v.enzdecaypersec ~= 0
        layout = addExternalReaction(layout,'enzyme_decay','enzyme[e]',-1,'vmax',v.enzdecaypersec);
    end
end

%set certain metabolites;
if isfield(v,'initO2')
    layout = setMedia(layout,'o2[e]',model.v.initO2);
    %layout = applyStaticMediaMask(layout,mask,'o2[e]',v.initO2);
end
if isfield(v,'initCO2')
    layout = setMedia(layout,'co2[e]',model.v.initCO2);
    %layout = applyStaticMediaMask(layout,mask,'co2[e]',v.initCO2);
end
if isfield(v,'initH')
    layout = setMedia(layout,'h[e]',model.v.initH);
end
if isfield(v,'initH2O')
    layout = setMedia(layout,'h2o[e]',model.v.initH2O);
    %layout = applyStaticMediaMask(layout,mask,'h2o[e]',v.initH2O);
end
if isfield(v,'initNH4')
    layout = setMedia(layout,'nh4[e]',model.v.initNH4);
    %layout = applyStaticMediaMask(layout,mask,'hn4[e]',v.initNH4);
end
if isfield(v,'initPi')
    layout = setMedia(layout,'pi[e]',model.v.initPi);
    %layout = applyStaticMediaMask(layout,mask,'pi[e]',v.initPi);
end
if isfield(v,'initSO4')
    layout = setMedia(layout,'so4[e]',model.v.initSO4);
    %layout = applyStaticMediaMask(layout,mask,'so4[e]',v.initSO4);
end
if isfield(v,'initK')
    layout = setMedia(layout,'k[e]',model.v.initK);
    %layout = applyStaticMediaMask(layout,mask,'k[e]',v.initK);
end
if isfield(v,'initergst')
    layout = setMedia(layout,'ergst[e]',model.v.initergst);
    %layout = applyStaticMediaMask(layout,mask,'k[e]',v.initK);
end
if isfield(v,'initzymst')
    layout = setMedia(layout,'zymst[e]',model.v.initzymst);
    %layout = applyStaticMediaMask(layout,mask,'k[e]',v.initK);
end
if isfield(v,'initoleate')
    layout = setMedia(layout,'oleate[e]',model.v.initoleate);
    %layout = applyStaticMediaMask(layout,mask,'k[e]',v.initK);
end
if isfield(v,'initpalmitoleate')
    layout = setMedia(layout,'palmitoleate[e]',model.v.initoleate);
    %layout = applyStaticMediaMask(layout,mask,'k[e]',v.initK);
end
if isfield(v,'initgthox')
    layout = setMedia(layout,'gthox[e]',v.initgthox);
end
if isfield(v,'initman')
    layout = setMedia(layout,'man-D[e]',v.initman);
end
%if isfield(v,'inittre')
%    layout = setMedia(layout,'tre[e]',v.inittre);
%end
if isfield(v,'trehalose_per_g')
    layout = setMedia(layout,'tre[e]',v.initialpop * v.trehalose_per_g);
end

if isfield('initgly_pct',v)
    if v.initgly_pct > 0
        layout = applyInitGlycogen(layout,v);
    end
end

layout.params.maxSpaceBiomass = v.maxspacebiomass;
layout.params.writeBiomassLog = v.writebiomasslog;
layout.params.writeMediaLog = v.writemedialog;
layout.params.writeFluxLog = v.writefluxlog;
layout.params.biomassLogRate = v.lograte;
layout.params.mediaLogRate = v.lograte;
layout.params.fluxLogRate = v.lograte;
layout.params.biomassLogName = [v.filepath '\log_biomass.m'];
layout.params.mediaLogName = [v.filepath '\log_media.m'];
layout.params.fluxLogName = [v.filepath '\log_flux.m'];
layout.params.maxCycles = v.maxcycles;
layout.params.timeStep = v.timestep;
layout.params.pixelScale = 10;
layout.params.spaceWidth = v.spaceWidth;
layout.params.deathRate = v.deathrate;

if isfield(v,'continuous')
    if v.continuous
        %layout is a bioreactor for continuous culture
        
        %add dilution for all metabolites
        %units are h^-1
        for i = 1:length(layout.mets)
            layout = addExternalReaction(layout,['DILUTE_' layout.mets{i}],layout.mets{i},-1,'vmax',v.dilutionrate);
        end
        %add dilution for organisms
        layout.params.deathRate = v.dilutionrate + (layout.params.deathRate * (1-v.dilutionrate));
    end
end

end