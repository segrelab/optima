function [figure2, models] = plotYAModels(yamodels)
% Plot yield over alpha
models = executeModels(yamodels);
    
biodat = zeros(length(models),1);
xdat = zeros(length(models),1);
for i = 1:length(models)
    biodat(i) = log10(models{i}.biomass(end));
    xdat(i) = models{i}.v.alpha;
end

figure2 = figure;

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');
% Create plot
plot(xdat,biodat,'DisplayName','Final Biomass','Marker','diamond');

% Create labels
xlabel('Alpha (Enzyme:Biomass ratio)');
ylabel('log10(Final Biomass');
title('Final biomass with varying Alpha');

%ylim(axes2,[0 0.005]);
box(axes2,'on');
end