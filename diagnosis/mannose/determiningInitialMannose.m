%See also: mannose_initialization.m

%A script to figure out how much mannose should be added to the initial
%media for replicating the DenHaan experiment

%Requires execution of pipeline0.mlx before runnning

%Created March 4 2020

%trying to determine an amount that is large enough to start cellulolysis,
%but small enough to be quickly depleted

%minimum mannose:biomass for growth is 100mmol:1g

%output saved as determiningInitialMannose_table.mat

%arbitrarily decided requirements:
t_mannose_depletion = 6; %max hours until mannose is almost entirely gone
mannose_depletion_threshold = 0.01; %how much remains when called depleted
t_lagphase = 48; %max hours until exponential growth should occur
lag_threshold = 0.02; %growth rate that we consider having moved out of stationary
yield_threshold = 2; %yield that must be achieved relative to without cellulose in order to consider growth a success

%load('C:\sync\biomes\cellulose\optima\diagnosis\mannose\mannoselayout.mat');

%model = layout.models{1};
%model = anaerobicModelWithNGAM(model,0.7,23.49);
%as 'layout'

mannose_weight = 180.156; %g/mol

% ions = {'ca2[e]' 'cl[e]' 'cu2[e]' 'mg2[e]' 'mn2[e]' 'na1[e]' 'zn2[e]'};
% for i = 1:length(ions)
%     layout = setMedia(layout,ions{i},v.initmedia);
% end

scales = [1000 200 150 100 75 20  5 2 1 .5 .1 .05 .01];
scales = [1000 200 100 25 5 1];
scales = [1000 750 500 250 100];
scales = 1000:-50:100;
% scales = 1000;
scales = [25 20 15 10 5 2 1];

t = table();
v = args.p2.v;
v_default = v;
v.maxcycles = t_lagphase * 4;%1.5;
v.maxcycles = 1000;
% v.maxcycles = 300;

% %temp: turn down mannose uptake to reasonable levels
% v.vmax_man = 1e10;
% v.km_man = 1e-15;
% v.vmax_man = v_default.vmax_glc * 10;
% v.km_man = v_default.km_glc / 10;
v.vmax_man = v_default.vmax_glc;
v.km_man = v_default.km_glc;

% % %temp: scale oleates by 100
%  v.initoleate = v.initoleate * 100;
%  v.initpalmitoleate = v.initpalmitoleate * 100;
% 
% % %temp: scale cellulose
%  v.initcellulose = v_default.initcellulose * 100;
% 
% %temp: scale ergosterol
% v.initergst = v_default.initergst * 100;

%temp:
% args.p2.v.NGAM = 0;
% args.p2.v.GAM = 30.49;
% v.NGAM = 0;
% v.GAM = 30.49;


%alphas = v.alpha * [10 3 1 0.3 0.1];
alphas = v.alpha;

for k = 1:length(alphas)
    v.alpha = alphas(k);
    
    model = prepareYeastGEMModel('v',v);
    model.v = v;
    
    layout = buildCelOptLayout(model);
      
    layout_default = layout;
    
    for i = 1:length(scales)
        temp = table;
        layout = layout_default;
        v.mannose_pct = scales(i);
        v.initman = v.initialpop * v.mannose_pct * 1000 / mannose_weight;
        layout = setMedia(layout,'man-D[e]',v.initman);
        model = layout.models{1};
        runComets(layout);
        bio = parseBiomassLog(layout.params.biomassLogName);
        media = parseMediaLog(layout.params.mediaLogName);
        temp.scale(1) = scales(i);
        temp.alpha(1) = alphas(k);
        temp.yield_g(1) = max(bio.biomass) - min(bio.biomass);
        
        %celamt = layout_amt(stridx('cellulose',layout.mets));
        layout = setMedia(layout,'cellulose',0);
        layout.params.writeMediaLog = false;
        runComets(layout);
        bio_poor = parseBiomassLog(layout.params.biomassLogName);
        temp.yield_g_fromManOnly(1) = max(bio_poor.biomass) - min(bio_poor.biomass);
        temp.yield_g_fromGlc(1) = temp.yield_g(1) - temp.yield_g_fromManOnly(1);
        y = bio.biomass - min(bio.biomass);
        y = y / temp.yield_g_fromManOnly(1);
        lag_yield = bio.t(y >= yield_threshold);
        temp.lag_yield(1) = nan;
        if any(lag_yield)
            temp.lag_yield(1) = lag_yield(1);
        end
        
        %lag phase length
        temp.biodelta{1} = bio.biomass(2:end) - bio.biomass(1:end-1);
        temp.biodelta_relative{1} = temp.biodelta{1} ./ bio.biomass(1:end-1);
        lag = bio.t(temp.biodelta_relative{1} >= lag_threshold);
        temp.lag_rate(1) = nan;
        if any(lag)
            temp.lag_rate(1) = lag(1);
        end
        %--
        
        %time to mannose depletion
        media_man = media(strcmp('man-D[e]',media.metname),:);
        man_amt = media_man.amt;
        man_init = man_amt(1);
        man_relative = man_amt / man_init;
        depletion = media_man.t(man_relative <= mannose_depletion_threshold);
        temp.depletion(1) = nan;
        if any(depletion)
            temp.depletion(1) = depletion(1);
        end
        
        temp.bio{1} = bio;
        temp.media{1} = media;
        
        temp.pass_lag(1) = ~isnan(temp.lag_rate(1)) && temp.lag_rate(1) <= t_lagphase;
        temp.pass_depletion(1) = ~isnan(temp.depletion(1)) && temp.depletion(1) <= t_mannose_depletion;
        temp.pass_yield(1) = ~isnan(temp.lag_yield(1)) && temp.lag_yield(1) <= t_lagphase;
        if ~any(size(t))
            t = temp;
        else
            t = vertcat(t,temp);
        end
    end
end
mannosetable = t;

%% 
% 
% xmax = 500;
% rnum = 4;
% 
% figure();
% plotBiomassTimecourse(mannosetable.bio{rnum},{'yeastGEMxml_noMaintRxn.txt'},'l');
% xlim([0 xmax]);
% figure();
% plotMediaTimecourse(mannosetable.media{rnum},{'man-D[e]','glc-D[e]','cellulose'},'l');
% xlim([0 xmax]);
