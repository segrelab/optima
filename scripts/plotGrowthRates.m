function [fig, models] = plotGrowthRates(models)
models = executeModels(models);

%assemble data
nrecords = length(models{1}.timestep);
ydat = zeros(nrecords,length(models));
ydat_scaled = zeros(nrecords,length(models));
tdat = zeros(nrecords,length(models));
for i = 1:length(models)
    deltas = models{i}.biomass(2:end) - models{i}.biomass(1:end-1); 
    b = models{i}.biomass(1:end-1);
    %deltas = [0.0; deltas];
    tt = models{i}.timestep * models{i}.v.timestep;
    tdeltas = tt(2:end) - tt(1:end-1);
    deltas = deltas./tdeltas;
    tdat(1:length(tt),i) = tt;
    tdat(length(tt)+1:end,i) = max(tt);
    ydat(2:length(deltas)+1,i) = deltas;
    deltas_scaled = deltas./b;
    ydat_scaled(2:length(deltas_scaled)+1,i) = deltas_scaled;
end

fig = figure;

% Create the scaled log(biomass)/time plot
subplot(2,1,1);
axes1 = gca;
hold(axes1,'on');
plot1 = plot(tdat,ydat_scaled,'Parent',axes1);
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
ylabel('GDW/GDW/h','HorizontalAlignment','center');
title('Biomass Production per hour per GDW');
box(axes1,'on');
legend1 = legend(axes1,'show');
set(legend1,'Location','eastoutside');
hold(axes1,'off');

% Create the total log(biomass)/time plot
subplot(2,1,2);
axes2 = gca;
hold(axes2,'on');
plot2 = plot(tdat,ydat,'Parent',axes2);
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
ylabel('GDW/h','HorizontalAlignment','center');
title('Total Biomass Production Rate Per Hour');
box(axes2,'on');
legend2 = legend(axes2,'show');
set(legend2,'Location','eastoutside');
hold(axes2,'off');

end