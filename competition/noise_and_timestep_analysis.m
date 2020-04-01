
%TODO: add cell death

cdat_temp = sortrows(cdat,{'maint_replaced','decayrate','alpha'});
cdat_temp = cdat_temp(cdat_temp.alpha == exp(-6),:);
glc = cdat_temp.glc_amt;
t = cdat_temp.t;
[nrows, ncols] = size(glc);
prod = cell(nrows,1);
cheat = cell(nrows,1);
for i = 1:nrows
    prod{i} = cdat_temp.biomass_producer{i};
    delta_prod{i} = prod{i}(2:end) - prod{i}(1:end-1);
    cheat{i} = cdat_temp.biomass_cheater{i};
    delta_cheat{i} = cheat{i}(2:end) - cheat{i}(1:end-1);
end

colorset = jet(nrows);
colorset = vertcat(colorset,colorset);
legendstrs = cell(1,nrows);
linetypes = {'-' '-' '-' '-' '-' '-' '-' '--' '--' '--' '--' '--' '--' '--' };

fig1 = figure();
hold on;
for i = 1:nrows
    plot(t{i},glc{i},linetypes{i},'color',colorset(i,:));
    lam = 'na';
    if cdat_temp.decayrate(i) > 0
        lam = num2str(1/cdat_temp.decayrate(i));
    end
    legendstrs{i} = ['logA:' num2str(log(cdat_temp.alpha(i))) ';Lambda:' lam];
end
title('Glc Amount');
hold off;

legendstrs = {'No Decay, Orig Maint' 'Decay, Orig Maint' 'No Decay, Maint Replaced' 'Decay, Maint Replaced'};

leg1 = legend(legendstrs,'location','eastoutside');

fig2 = figure();
hold on;
for i = 1:nrows
    plot(t{i}(2:end),delta_prod{i},linetypes{i},'color',colorset(i,:));
end
title('Enzyme Producer Growth Rate');
leg2 = legend(legendstrs,'location','eastoutside');
hold off;

fig3 = figure();
hold on;
for i = 1:nrows
    plot(t{i}(2:end),delta_cheat{i},linetypes{i},'color',colorset(i,:));
end
title('Cheater Growth Rate');
leg3 = legend(legendstrs,'location','eastoutside');
hold off;

fig4 = figure();
hold on;
for i = 1:nrows
    plot(t{i},prod{i},linetypes{i},'color',colorset(i,:));
end
title('Enzyme Producer Biomass');
leg4 = legend(legendstrs,'location','eastoutside');
hold off;

fig5 = figure();
hold on;
for i = 1:nrows
    plot(t{i},cheat{i},linetypes{i},'color',colorset(i,:));
end
title('Cheater Biomass');
leg5 = legend(legendstrs,'location','eastoutside');
hold off;

fig6 = figure();
hold on;
for i = 1:nrows
    plot(t{i},cheat{i}./prod{i},'color',colorset(i,:));
end
title('Cheater/Producer Biomass Ratio')
leg6 = legend(legendstrs,'location','eastoutside');
hold off;

%% compare yield when timestep is 1h vs 30 min
cdat_half = sortrows(cdat,{'maint_replaced','decayrate','alpha'});
fig7 = figure();
hold on;
for i = 1:4
    plot(cdat1h.biomass_producer{i}(1:end-1) ./ cdathalf.biomass_producer{i}(1:2:end-1));
end
leg7 = legend(legendstrs,'location','eastoutside');
title('Producer biomass @ ts = 1h / biomass @ ts = 0.5h)');
hold off;

fig8 = figure();
hold on;
for i = 1:4
    plot(cdat1h.biomass_cheater{i}(1:end-1) ./ cdathalf.biomass_cheater{i}(1:2:end-1));
end
leg8 = legend(legendstrs,'location','eastoutside');
title('Cheater biomass (ts = 1h / ts = 0.5h)');
hold off;

fig9 = figure();
hold on;
for i = 1:4
    plot(cdat1h.glc_amt{i}(1:end-1) ./ cdathalf.glc_amt{i}(1:2:end-1));
end
leg9 = legend(legendstrs,'location','eastoutside');
title('Glucose concentration (ts = 1h / ts = 0.5h)');
hold off;
%%
fig10 = figure();
hold on;
x = 3;
plot(cdat1h.t{x},cdat1h.biomass_producer{x},cdat1h.t{x},cdat1h.biomass_cheater{x},...
    cdathalf.t{x},cdathalf.biomass_producer{x},cdathalf.t{x},cdathalf.biomass_cheater{x})
title('Biomass at timestep = 1h and 1/2h');
leg10 = legend({'Producer 1h' 'Cheater 1h' 'Producer 1/2h' 'Cheater 1/2h'});
hold off;