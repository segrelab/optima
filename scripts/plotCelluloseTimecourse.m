function [fig, models] = plotCelluloseTimecourse(models)
models = executeModels(models);

%assemble data
nrecords = length(models{1}.timestep);
ydat = zeros(nrecords,length(models));
tdat = zeros(nrecords,length(models));
for i = 1:length(models)
    npoints = length(models{i}.cellulose_amt);
    ydat(1:npoints,i) = models{i}.cellulose_amt;
    tdat(1:npoints,i) = models{i}.timestep * models{i}.v.timestep;
end

fig = figure;
% Create axes
axes1 = axes('Parent',fig,...
    'Position',[0.128692810457516 0.11 0.775 0.815]);
hold(axes1,'on');

% Create the log(biomass)/time plot
plot1 = plot(tdat,log10(ydat),'Parent',axes1);
for i = 1:length(models)
    m = models{i};
    s = '-'; %solid line
    w = 0.75; %line width. Default is 0.5
%     if m.v.alpha < 1
%         s = '--'; %dashed line
%         w = 0.5;
%     end
    set(plot1(i),'DisplayName',m.description,'LineStyle',s,'LineWidth',w);
end
xlabel('Timestep (Hours)','HorizontalAlignment','center');
ylabel('log10(Cellulose)','HorizontalAlignment','center');
title('Cellulose Amount with varying Alpha');

xlim(axes1,[min(tdat(:)) max(tdat(:))]);
ylim(axes1,[-6 1]);
box(axes1,'on');
legend1 = legend(axes1,'show');
set(legend1,'Location','eastoutside');

%add a line at the half-saturation concentration
kmcon = models{1}.v.km_cel;
width = models{1}.v.spaceWidth;
vol = width ^ 3;
km = kmcon * vol;
rl = refline(0,log10(km));
rl.LineStyle = ':';
rl.DisplayName = 'KM';
rl.LineWidth = 1.5;

end