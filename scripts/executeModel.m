%run the model and attach its biomass and media log to it
function model = executeModel(model)
layout = buildCelOptLayout(model);
createCometsFiles(layout,model.v.filepath,'layout_temp.txt');
runCometsOnDirectory(model.v.filepath);

tab = parseBiomassLog([model.v.filepath '\log_biomass.m']);
model.biomass = tab.biomass; %no need to sum if there's no biomass diffusion
model.timestep = tab.t;

if model.v.writemedialog
    media = parseMediaLog([model.v.filepath '\log_media.m'],{'enzyme_e', 'cellulose', 'glc__D_e'});
    model.enzyme_amt = zeros(length(model.timestep),1);
    model.cellulose_amt = zeros(length(model.timestep),1);
    model.glc_amt = zeros(length(model.timestep),1);
    for i = 1:length(model.timestep)
        time = model.timestep(i);
        enzrows = (media.t == time) & strcmp(media.metname,'enzyme_e');
        celrows = (media.t == time) & strcmp(media.metname,'cellulose');
        glcrows = (media.t == time) & strcmp(media.metname,'glc__D_e');
        subtab = media(enzrows,:);
        model.enzyme_amt(i) = sum(subtab.amt);
        subtab = media(celrows,:);
        model.cellulose_amt(i) = sum(subtab.amt);
        subtab = media(glcrows,:);
        model.glc_amt(i) = sum(subtab.amt);
    end
end

end