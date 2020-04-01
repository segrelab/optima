function [parent,model] = plotAllTimecourses(model,scalemode,parent)
%PLOTALLTIMECOURSES plot the timecourse for glc, cellulose, enzyme, and
%biomass in one figure
%   If a second argument is given, it will be used as the parent of the
%   figure. Otherwise a new figure is returned
if nargin < 2
    scalemode = ''; %blank for no scaling or one of {'pctmax' 'log10'}
end
ylabtext = 'Concentration (mmol)';

if ~isfield(model,'biomass')
    model = executeModel(model)
end

if nargin < 3
    fig = figure;
    parent = gca;
end

%assemble data
biodat = model.biomass;
tdat = model.timestep * model.v.timestep;
edat = model.enzyme_amt;
cdat = model.cellulose_amt;
gdat = model.glc_amt;
ethdat = model.etoh_amt;

%scale data to percent of max
if strcmp(scalemode,'pctmax')
    edat = edat / max(edat);
    cdat = cdat / max(cdat);
    gdat = gdat / max(gdat);
    ethdat = ethdat / max(ethdat);
    ylabtext = '% Max Concentration';
end

if strcmp(scalemode,'log10')
    edat = log10(edat);
    cdat = log10(cdat);
    gdat = log10(gdat);
    ethdat = log10(ethdat);
    ylabtext = 'Concentration (mmol) (log10)';
end

%plot the metabolites
hold on;
plot(tdat,cdat,'Parent',parent);
plot(tdat,edat,'Parent',parent);
plot(tdat,gdat,'Parent',parent);
plot(tdat,ethdat,'Parent',parent);


xlabel('Timestep (Hours)','HorizontalAlignment','center');
ylabel(ylabtext,'HorizontalAlignment','center');
%title('');

box(parent,'on');
%legend1 = legend(parent,'show');
legend1 = legend('Cellulose', 'Enzyme', 'Glucose', 'Ethanol');
set(legend1,'Location','eastoutside');


hold off;
end

