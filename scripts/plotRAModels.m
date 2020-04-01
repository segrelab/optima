function [figure3, models] = plotRAModels(ramodels)
%plot k as in the exponential growth formula x(t) = x(0)e^kt, for a set of
%models with varying alpha
models = executeModels(ramodels);

biodat = zeros(length(models),1);
xdat = zeros(length(models),1);
cel_amts = zeros(length(models),1);
for i = 1:length(models)
    xdat(i) = models{i}.v.alpha;
    %get the rate a few steps before the curve levels off
    model = models{i};
    deltas = log(model.biomass(2:end)) - log(model.biomass(1:end-1));
    %steps = round(finish/10); %step back 1/10th of the way
    %drange = 1; %use this to check a few different deltas
    %a = max([finish-steps-drange 1]);
    %b = min ([finish-steps+drange length(deltas)]);
    %d = mean(deltas(a:b)); 
    %biodat(i) = d;
    biodat(i) = max(deltas);
    cel_amts(i) = models{i}.v.initcellulose;
end
figure3 = figure;

% Create axes
axes3 = axes('Parent',figure3);
hold(axes3,'on');
% Create plot
ucels = unique(cel_amts);
for i = 1:length(ucels)
    c = ucels(i);
    idx = cel_amts == c;
    plot(xdat(idx),biodat(idx),'DisplayName','Rate','Marker','diamond');
end
    
% Create labels
xlabel('Alpha (Enzyme:Biomass ratio)');
ylabel('Log-phase growth rate');
title('Exponential growth rate with varying Alpha');
lgd = legend(arrayfun(@(x) string(x),ucels),'Location','best');
title(lgd,'Initial Cellulose (mmol)');
% Uncomment the following line to preserve the Y-limits of the axes
%ylim(axes3,[0 3.5]);
box(axes3,'on');

hold(axes3,'off');
end