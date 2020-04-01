function [fig, models] = plotGlcRates(models)
models = executeModels(models);

%assemble data
nrecords = length(models{1}.timestep);
ydat_scaled = zeros(nrecords,length(models));
ydat = zeros(nrecords,length(models));
tdat = zeros(nrecords,length(models));
for i = 1:length(models)
    npoints = length(models{i}.glc_amt);
    deltas = models{i}.glc_amt(2:npoints) - models{i}.glc_amt(1:npoints-1); 
    b = models{i}.biomass(1:npoints-1);
    tt = models{i}.timestep * models{i}.v.timestep;
    tdat(1:npoints,i) = tt;
    tdat(length(tt)+1:end,i) = max(tt);
    tdeltas = tt(2:npoints) - tt(1:npoints-1);
    deltas = deltas./tdeltas;
    ydat(2:npoints,i) = deltas;
    deltas_scaled = deltas./b;
    ydat_scaled(2:npoints,i) = deltas_scaled;
end

fig = figure;

% Create the scaled log(biomass)/time plot
subplot(2,1,1);
% Create axes
%axes1 = axes('Parent',fig,...
%    'Position',[0.128692810457516 0.11 0.775 0.815]);
%axes1 = axes('Parent',fig);
axes1 = gca;
hold(axes1,'on');
plot1 = plot(tdat,ydat_scaled,'Parent',axes1);
%plot1 = plot(tdat,ydat_scaled);
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
ylabel('mmol Glucose / GDW * h','HorizontalAlignment','center');
title('Glucose Production/Consumption Rate per GDW');
box(axes1,'on');
legend1 = legend(axes1,'show');
set(legend1,'Location','eastoutside');
hold(axes1,'off');
% Create the total log(biomass)/time plot
subplot(2,1,2);
% Create axes
%axes2 = axes('Parent',fig,...
%    'Position',[0.128692810457516 0.11 0.775 0.815]);
%axes2 = axes('Parent',fig);
axes2 = gca;
hold(axes2,'on');
plot2 = plot(tdat,ydat,'Parent',axes2);
%plot2 = plot(tdat,ydat);
for i = 1:length(models)
    m = models{i};
    s = '-'; %solid line
    w = 0.75; %line width. Default is 0.5
%     if m.v.alpha < 1
%         s = '--'; %dashed line
%         w = 0.5;
%     end
    set(plot2(i),'DisplayName',m.description,'LineStyle',s,'LineWidth',w);
end
xlabel('Timestep (Hours)','HorizontalAlignment','center');
ylabel('mmol Glucose / h','HorizontalAlignment','center');
title('Net Glucose Production/Consumption Rate');
box(axes2,'on');
legend2 = legend(axes2,'show');
set(legend2,'Location','eastoutside');
hold(axes2,'off');

end