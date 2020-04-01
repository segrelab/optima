function [figure10, models] = plotRateOverYield(btmodels)
%plot k as in the exponential growth formula x(t) = x(0)e^kt
models = executeModels(btmodels);

biodat = zeros(length(models),1);
xdat = zeros(length(models),1);
alphas = zeros(length(models),1);
for i = 1:length(models)
    alphas(i) = models{i}.v.alpha;
    %get the rate a few steps before the curve levels off
    model = models{i};
    deltas = log(model.biomass(2:end)) - log(model.biomass(1:end-1));
    finish = find(deltas > 0, 1, 'last');
    if finish == length(deltas)
        finish = length(deltas) - round(length(deltas)/10);
    end
    %steps = round(finish/10); %step back 1/10th of the way
    %drange = 1; %use this to check a few different deltas
    %a = max([finish-steps-drange 1]);
    %b = min ([finish-steps+drange length(deltas)]);
    %d = mean(deltas(a:b)); 
    %biodat(i) = d;
    biodat(i) = log10(max(deltas));
    xdat(i) = log10(model.biomass(end));
end
figure10 = figure;

% Create axes
axes10 = axes('Parent',figure10);
hold(axes10,'on');
% Create plot
plot(xdat,biodat,'DisplayName','Rate','Marker','diamond');

datalabels = cellstr( num2str(alphas) );
text(xdat,biodat, datalabels, 'VerticalAlignment','bottom','HorizontalAlignment','left');

% Create labels
xlabel('Final Biomass ( Log10 GDW )');
ylabel('Log10 (Max Growth Rate');
title('Growth Rate Over Final Yield [Labels = alpha]');

% Uncomment the following line to preserve the Y-limits of the axes
%ylim(axes3,[0 3.5]);
box(axes10,'on');

end