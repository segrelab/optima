function [parent, models] = plotBiomassOverTime(models,parent)
models = executeModels(models);

if nargin < 2
    fig = figure;
    parent = gca;
end

logmode = false;

%assemble data
nrecords = 1 + (models{1}.v.maxcycles / models{1}.v.lograte);
biodat = zeros(nrecords,length(models));
tdat = zeros(nrecords,length(models));
for i = 1:length(models)
    biodat(1:length(models{i}.biomass),i) = models{i}.biomass;
    tt = models{i}.timestep * models{i}.v.timestep;
    tdat(1:length(tt),i) = tt;
end

colors = jet(length(models));

% Create axes
%axes1 = axes('Parent',fig);%,'Position',[0.128692810457516 0.11 0.775 0.815]);
hold(parent,'on');

% Create the log(biomass)/time plot-
if logmode
    biodat = log10(biodat);
end
plot1 = plot(tdat,biodat,'Parent',parent,'lineWidth',2);
for i = 1:length(models)
    m = models{i};
    s = '-'; %solid line
    w = 0.75; %line width. Default is 0.5
%     if m.v.alpha < 1
%         s = '--'; %dashed line
%         w = 0.5;
%     end
    set(plot1(i),'DisplayName',m.description,'LineStyle',s,'LineWidth',w,'Color',colors(i,:));
end
xlabel('Timestep (Hours)','HorizontalAlignment','center');
ylabel('log10(Biomass)','HorizontalAlignment','center');
title('Total Biomass over Time');
box(parent,'on');
legend1 = legend(parent,'show');
set(legend1,'Location','eastoutside');

end