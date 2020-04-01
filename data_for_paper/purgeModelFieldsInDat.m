function dat = purgeModelFieldsInDat(dat)
%PURGEMODELFIELDSINDAT Clears a bunch of fields in the models in the given
%data table to reduce its size
%   Assumes the table contains a column "layout" and each layout is a
%   struct containing a cell array "models" containing one model

[nrows, ncols] = size(dat);
for i = 1:nrows
    layout = dat.layout{i};
    model = layout.models{1};
    
    %make a new struct "m" with just what we want to keep
    m.v = model.v;
    m.description = model.description;
    m.enzExpressionMode = model.enzExpressionMode;
    m.biomass = model.biomass;
    m.timestep = model.timestep;
    m.enzyme_amt = model.enzyme_amt;
    m.cellulose_amt = model.cellulose_amt;
    m.glc_amt = model.glc_amt;
    m.acetate_amt = model.acetate_amt;
    m.etoh_amt = model.etoh_amt;
    
    layout.models{1} = m;
    dat.layout{i} = layout;
end

end

