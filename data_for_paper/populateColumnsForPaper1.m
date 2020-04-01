function dat = populateColumnsForPaper1(dat)
%POPULATECOLUMNSFORPAPER1 Adds columns needed for plots_enzymepaper1b.mlx
%   Adds max growth rate, max enzyme concentration, max cellulose digestion
%   rate, lag phase length

[rows, cols] = size(dat);

%add growth rate and max growth rate
mu = zeros(rows,1);
for i = 1:rows
    biomass = dat.layout{i}.models{1}.biomass;
    deltas = biomass(2:end) - biomass(1:end-1);
    deltas = deltas ./ biomass(1:end-1); %scale to the amount of biomass present
    timestep = dat.layout{i}.models{1}.v.timestep;
    deltas = deltas / timestep;
    mu(i) = max(deltas);
    dat.layout{i}.models{1}.growthrate = deltas; 
end
dat.growthrate_max = mu;


%add max enzyme concentration
ec = zeros(rows,1);
for i = 1:rows
    enzyme = dat.layout{i}.models{1}.enzyme_amt;
    ec(i) = max(enzyme);
end
dat.enzyme_max = ec;

%add max cellulose digestion rate
vc = zeros(rows,1);
for i = 1:rows
    cel = dat.layout{i}.models{1}.cellulose_amt;
    deltas = cel(2:end) - cel(1:end-1);
    deltas = -deltas;
    timestep = dat.layout{i}.models{1}.v.timestep;
    deltas = deltas / timestep;
    vc(i) = max(deltas);
end
dat.cellulose_negvmax = vc;

%add lag phase length
%   Lag phase ends when the second derivative of biomass goes from convex
%   to concave
lag_steps = zeros(rows,1);
lag_hours = zeros(rows,1);
for i = 1:rows
    biomass = dat.layout{i}.models{1}.biomass;
    deltas = biomass(2:end) - biomass(1:end-1);
    d2 = deltas(2:end) - deltas(1:end-1);
    d2(1) = 0; %this point's a bit messed up...
    pos = find(d2 < 0);
    if length(pos) < 1
        pos = length(d2);
    end
    timestep = dat.layout{i}.models{1}.v.timestep;
    lag_steps(i) = pos(1);
    lag_hours(i) = pos(1) * timestep;
end
dat.lag_steps = lag_steps;
dat.lag_hours = lag_hours;

%did all carbon get consumed?
finished = zeros(rows,1);
for i = 1:rows
    finished(i) = (dat.layout{i}.models{1}.cellulose_amt(end) + dat.layout{i}.models{1}.glc_amt(end)) < 1e-8;
end
dat.sim_finished = finished; 

%update model.timestep to be the time in hours
for i = 1:rows
    t = 0:length(dat.layout{i}.models{1}.timestep)-1;
    t = t * dat.layout{i}.models{1}.v.timestep;
    dat.layout{i}.models{1}.timestep = t;
end

end

