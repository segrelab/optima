function [figure4, models] = plotRates(ramodels)
%the instantaneous growth rate of each model
models = executeModels(ramodels);

biodat = zeros(length(models{1}.timestep),length(models)); %y axis
xdat = zeros(length(models{1}.timestep),length(models));
for i = 1:length(models)
    xdat(1:end,i) = models{i}.timestep * models{i}.v.timestep;
    model = models{i};
    deltas = log(model.biomass(2:end)) - log(model.biomass(1:end-1));
    biodat(2:end,i) = deltas; %first delta is left as 0.
end
figure4 = figure;

% Create axes
axes4 = axes('Parent',figure4);
hold(axes4,'on');
% Create plot
plot4 = plot(xdat,biodat,'Parent',axes4);
for i = 1:length(models)
    m = models{i};
    s = '-'; %solid line
    w = 0.75; %line width. Default is 0.5
    if m.v.alpha < 1
        s = '--'; %dashed line
        w = 0.5;
    end
    set(plot4(i),'DisplayName',m.description,'LineStyle',s,'LineWidth',w);
end
xlabel('Timestep');
ylabel('Growth rate');
title('Growth rate over time with varying Alpha');
box(axes4,'on');
legend4 = legend(axes4,'show');
set(legend4,'Location','eastoutside');

% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 0.0001]);
box(axes4,'on');

end