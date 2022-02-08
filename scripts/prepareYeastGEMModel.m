function model = prepareYeastGEMModel(varargin)
%PREPAREYEAST7MODEL Loads the Yeast 8.3.2 model and modifies it so that... 
%   * it can run anaerobically
%   * metabolites are renamed to interact with existing scripts
%   * glucose uptake rate is set 
%
%Inputs:
%   'model' : COBRA model of yeast 8.3.2. It will be loaded from file if
%   not provided
%   'alpha' : mmol enzyme produced per gram biomass growth
%   'stoich': number of each amino acid in each molecule of enzyme
%   'seqlength': total length of the enzyme protein
%   'objname': RxnName of the objective reaction if you'd like to change it
%   'atp': number of atp consumed per peptide bond
%   'gtp': number of gtp consumed per peptide bond
%   'addEnzymeToMaint': Should the maintenance reaction create enzyme?
%   'maintRxnAsObjective': reframe the maintenance reaction as the
%           highest-priority objective instead of a rxn with lb > 0
%   'applyCorrections': Apply corrections according to yeastGEM/ComplimentaryScripts/modelCuration/modelCorrections
%   'relaxExchangeLbs': set lower bounds for any default exchange reactions
%           from 0->defaultLB

v.lb = -1000;
v.alpha = 1e-4;
modelloaded = false;

%T reesei endoglucanase EGI, Cel7b
%Source: https://genome.jp/dbget-bin/www_bget?tre:TRIREDRAFT_123989
%We use Cel7b because that's the one that has mass data on BRENDA
aaseq_EGI = 'MAPSVTLPLTTAILAIARLVAAQQPGTSTPEVHPKLTTYKCTKSGGCVAQDTSVVLDWNYRWMHDANYNSCTVNGGVNTTLCPDEATCGKNCFIEGVDYAASGVTTSGSSLTMNQYMPSSSGGYSSVSPRLYLLDSDGEYVMLKLNGQELSFDVDLSALPCGENGSLYLSQMDENGGANQYNTAGANYGSGYCDAQCPVQTWRNGTLNTSHQGFCCNEMDILEGNSRANALTPHSCTATACDSAGCGFNPYGSGYKSYYGPGDTVDTSKTFTIITQFNTDNGSPSGNLVSITRKYQQNGVDIPSAQPGGDTISSCPSASAYGGLATMGKALSSGMVLVFSIWNDNSQYMNWLDSGNAGPCSSTEGNPSNILANNPNTHVVFSNIRWGDIGSTTNSTAPPPPPASSTTFSTTRRSSTTSSSPSCTQTHWGQCGGIGYSGCKTCTSGTTCQYSNDYYSQCL';
%Saccharomycopsis fibuligeria beta-glucosidase
%Ref: Machida M.,Ohtsuki I.,Fukui S.,Yamashita I.	Nucleotide sequences of Saccharomycopsis fibuligera genes for extracellular beta-glucosidases as expressed in Saccharomyces cerevisiae.	Appl. Environ. Microbiol.	54	3147-3155	1988	3146949
aaseq_BGL1 = 'MLMIVQLLVFALGLAVAVPIQNYTQSPSQRDESSQWVSPHYYPTPQGGRLQDVWQEAYARAKAIVGQMTIVEKVNLTTGTGWQLDPCVGNTGSVPRFGIPNLCLQDGPLGVRFADFVTGYPSGLATGATFNKDLFLQRGQALGHEFNSKGVHIALGPAVGPLGVKARGGRNFEAFGSDPYLQGTAAAATIKGLQENNVMACVKHFIGNEQEKYRQPDDINPATNQTTKEAISANIPDRAMHALYLWPFADSVRAGVGSVMCSYNRVNNTYACENSYMMNHLLKEELGFQGFVVSDWGAQLSGVYSAISGLDMSMPGEVYGGWNTGTSFWGQNLTKAIYNETVPIERLDDMATRILAALYATNSFPTEDHLPNFSSWTTKEYGNKYYADNTTEIVKVNYNVDPSNDFTEDTALKVAEESIVLLKNENNTLPISPEKAKRLLLSGIAAGPDPIGYQCEDQSCTNGALFQGWGSGSVGSPKYQVTPFEEISYLARKNKMQFDYIRESYDLAQVTKVASDAHLSIVVVSAASGEGYITVDGNQGDRKNLTLWNNGDKLIETVAENCANTVVVVTSTGQINFEGFADHPNVTAIVWAGPLGDRSGTAIANILFGKANPSGHLPFTIAKTDDDYIPIETYSPSSGEPEDNHLVENDLLVDYRYFEEKNIEPRYAFGYGLSYNEYEVSNAKVSAAKKVDEELPEPATYLSEFSYQNAKDSKNPSDAFAPADLNRVNEYLYPYLDSNVTLKDGNYEYPDGYSTEQRTTPNQPGGGLGGNDALWEVAYNSTDKFVPQGNSTDKFVPQLYLKHPEDGKFETPIQLRGFEKVELSPGEKKTVDLRLLRRDLSVWDTTRQSWIVESGTYEALIGVAVNDIKTSVLFTI';
aaseq = [aaseq_EGI aaseq_BGL1];
v.stoich = zeros(1,20);
v.stoich(1) = count(aaseq,'A');
v.stoich(2) = count(aaseq,'R');
v.stoich(3) = count(aaseq,'N');
v.stoich(4) = count(aaseq,'D');
v.stoich(5) = count(aaseq,'C');
v.stoich(6) = count(aaseq,'Q');
v.stoich(7) = count(aaseq,'E');
v.stoich(8) = count(aaseq,'G');
v.stoich(9) = count(aaseq,'H');
v.stoich(10) = count(aaseq,'I');
v.stoich(11) = count(aaseq,'L');
v.stoich(12) = count(aaseq,'K');
v.stoich(13) = count(aaseq,'M');
v.stoich(14) = count(aaseq,'F');
v.stoich(15) = count(aaseq,'P');
v.stoich(16) = count(aaseq,'S');
v.stoich(17) = count(aaseq,'T');
v.stoich(18) = count(aaseq,'W');
v.stoich(19) = count(aaseq,'Y');
v.stoich(20) = count(aaseq,'V');

v.seqlength = sum(v.stoich);

v.NGAM = 0;
v.GAM = 30.49;

v.objname = 'growth';
v.changeobj = false;
v.atpperbond = 8;
v.gtpperbond = 4;
v.costfactor = 1;
v.vmax_glc = 0;
v.vmax_man = 0;
v.vmax_tre = 0;
v.vmax_gly = 0;
v.km_glc = nan;
v.km_man = nan;
v.km_tre = nan;
v.km_gly = nan;
v.maintRxnAsObjective = false;
v.addEnzymeToMaint = false;
v.anaerobic = false;
v.applyCorrections = false; %deprecated in some version sub 8.3.4
v.relaxExchangeLbs = false; %set lower bounds for any default exchange reactions from 0->lb
for i = 1:2:length(varargin)
    label = varargin{i};
    val = varargin{i+1};
    if strcmpi('model',label)
        model = val;
        modelloaded = true;
    end
    if strcmpi('alpha',label)
        v.alpha = val;
    end
    if strcmpi('stoich',label)
        v.stoich = val;
        v.seqlength = sum(val);
    end
    if strcmpi('seqlength',label)
        v.seqlength = val;
    end
    if strcmpi('objname',label)
        v.changeobj = true;
        v.objname = val;
    end
    if strcmpi('atp',label)
        v.atpperbond = val;
    end
    if strcmpi('gtp',label)
        v.gtpperbond = val;
    end
    if strcmpi('costfactor',label)
        v.costfactor = val;
    end
    if strcmpi('vmax_glc',label)
        v.vmax_glc = abs(val);
    end
    if strcmpi('maintRxnAsObjective',label)
        v.maintRxnAsObjective = val;
    end
    if strcmpi('addEnzymeToMaint',label)
        v.addEnzymeToMaint = val;
    end
    if strcmpi('makeAnaerobic',label)
        v.anaerobic = val;
    end
    if strcmpi('anaerobic',label)
        v.anaerobic = val;
    end
    if strcmpi('applyCorrections',label)
        v.applyCorrections = val;
    end
    if strcmpi('relaxExchangeLbs',label)
        v.relaxExchangeLbs = val;
    end
    if strcmpi('vmax_man',label)
        v.vmax_man =  abs(val);
    end
    if strcmpi('vmax_tre',label)
        v.vmax_tre =  abs(val);
    end
    if strcmpi('vmax_gly',label)
        v.vmax_gly = abs(val);
    end
    if strcmpi('v',label)
        %         allfields = {'GAM' 'NGAM' 'alpha' 'seqlength' 'objname' 'atp' 'gtp' 'costfactor' 'vmax_glc' 'maintRxnAsObjective' 'addEnzymeToMaint' 'anaerobic' 'applyCorrections' 'relaxExchangeLbs' 'vmax_man' 'vmax_tre' 'vmax_gly' 'km_glc' 'km_man' 'km_tre' 'km_gly'};
        %         allfields = {'alpha' 'enzdecayperhour' 'deathrate' 'dilutionrate'...
        %             'mode' 'layouttype' 'diff_enz' 'costfactor' 'richmedia'...
        %             'defaultbound' 'kcat_cel' 'km_cel' 'n_enzymes' 'enz_weight'...
        %             'atp_per_peptide' 'gtp_per_peptide' 'vmax_glc' 'km_glc'...
        %             'glcpercellulose' 'glcuptakerate' 'xdim' 'ydim' 'diff_cel'...
        %             'diff_glc' 'initcellulose' 'initglc' 'initmedia'...
        %             'initialpop' 'initO2' 'initNH4' 'initPi' 'initSO4' 'initK'...
        %             'initCO2' 'initH2O' 'initergst' 'initzymst' 'initoleate'...
        %             'initpalmitoleate' 'maxcycles' 'timestep' 'filepath'...
        %             'writemedialog' 'writebiomasslog' 'writefluxlog' 'lograte'...
        %             'maxspacebiomass' 'spaceWidth' 'maintRxnAsObjective' 'anaerobic'...
        %             'initH' 'initgthox' 'initfe2' 'vmax_gly' 'gly_molarmass'...
        %             'initgly_pct' 'mass_glcmonomer' 'initenzyme' 'vmax_man'...
        %             'km_man' 'diff_man' 'mannose_pct' 'initman' 'NGAM' 'GAM'};
        allfields = fields(val);
        
        for i = 1:length(allfields)
            %            if isfield(val,allfields{i})
            v = setfield(v,allfields{i},getfield(val,allfields{i}));
            %            end
        end
        if isfield(val,'stoich')
            v.stoich = val.stoich;
            v.seqlength = sum(val.stoich);
        end
    end
end


v.atpstoich = v.atpperbond * v.seqlength;
v.gtpstoich = v.gtpperbond * v.seqlength;

if ~modelloaded
    %load('C:\sync\biomes\models\Scerevisiae_Yeast7.6_norm.mat');
    %load('C:\sync\biomes\models\yeastGEM_8.3.2.mat');
    %load 'C:\sync\biomes\models\yeast-GEM-8.3.5-dev\ModelFiles\mat\yeastGEM_8-3-4-dev' model;
    load 'C:\sync\biomes\models\yeast-GEM-8.3.5-dev\ModelFiles\mat\yeastGEM_8-3-5-dev' model;
    %loads as "model"
end

if ~isfield(model,'modifications')
    model.modifications = cell(0);
end
model.mets_original = model.mets;
model.modifications = {'Duplicated original names of mets into field mets_original.'};

if v.anaerobic
    model = anaerobicModelWithNGAM(model,v.NGAM,v.GAM);
    model.modifications = [model.modifications {['Convert to anaerobic mode with anaerobicModelWithGAM script. NGAM:' num2str(v.NGAM) ', GAM:' num2str(v.GAM)]}];
end

% % Deprecated for model version 8.3.4 
if v.applyCorrections
    modelCorrections;
    model.modifications = [model.modifications {'Apply corrections according to ''yeastGEM/ComplimentaryScripts/modelCuration/modelCorrections'' script'}];
end

if v.changeobj
    model.c(1:end) = 0;
    model.c(stridx(objname,model.rxnNames)) = 1;
    model.modifications = [model.modifications {['Changed objective reaction to ' objname]}];
end

if v.relaxExchangeLbs
    nlb = 0;
    exIdxs = find(findExchRxns(model));
    for i = 1:length(exIdxs)
        if model.lb(exIdxs(i)) == 0
            model.lb(exIdxs(i)) = lb;
            nlb = nlb + 1;
        end
    end
    if nlb > 0
    model.modifications = [model.modifications {['Changed lb of ' num2str(nlb) ' reactions from 0 to ' num2str(v.lb)]}];
    end
end

if ~isfield(model,'vmax')
    model.vmax = nan(size(model.rxns));
end

if abs(v.vmax_glc) > 0
    model.lb(stridx('D-glucose exchange' ,model.rxnNames)) = -abs(v.vmax_glc);
    model.vmax(stridx('D-glucose exchange' ,model.rxnNames)) = abs(v.vmax_glc);
    model.modifications = [model.modifications {['Changed max glucose uptake bound to ' num2str(-abs(v.vmax_glc)) ' mmol/g/h']}];
end

if abs(v.vmax_man) > 0
    model.lb(stridx('D-mannose exchange' ,model.rxnNames)) = -abs(v.vmax_man);
    model.vmax(stridx('D-mannose exchange' ,model.rxnNames)) = abs(v.vmax_man);
    model.modifications = [model.modifications {['Changed max mannose uptake bound to ' num2str(-abs(v.vmax_man)) ' mmol/g/h']}];
end

if abs(v.vmax_tre) > 0
    model.lb(stridx('trehalose exchange' ,model.rxnNames)) = -v.vmax_tre;
    model.vmax(stridx('trehalose exchange' ,model.rxnNames)) = abs(v.vmax_tre);
    model.modifications = [model.modifications {['Changed max trehalose uptake bound to ' num2str(-abs(v.vmax_tre)) ' mmol/g/h']}];
end

if abs(v.vmax_gly) > 0
    model = addExchangeRxn(model,{'s_0773[c]'},-abs(v.vmax_gly));
    model.rxnNames{end} = 'glycogen exchange';
    model.rxns{end} = 'EX_glycogen[c]';
    model.vmax(stridx('glycogen exchange' ,model.rxnNames)) = abs(v.vmax_gly);
    model.modifications = [model.modifications {['Added glycogen exchange with bound ' num2str(-abs(v.vmax_gly)) ' mmol/g/h']}];
end

if ~isfield(model,'km')
    model.km = nan(size(model.rxns));
end

if v.km_glc > 0
    model.km(stridx('D-glucose exchange' ,model.rxnNames)) = abs(v.km_glc);
    model.modifications = [model.modifications {['Changed glucose KM to ' num2str(abs(v.km_glc)) ' mmol']}];
end
if v.km_man > 0
    model.km(stridx('D-mannose exchange' ,model.rxnNames)) = abs(v.km_man);
    model.modifications = [model.modifications {['Changed mannose KM to ' num2str(abs(v.km_man)) ' mmol']}];
end
if v.km_gly > 0
    model.km(stridx('glycogen exchange' ,model.rxnNames)) = abs(v.km_gly);
    model.modifications = [model.modifications {['Changed glycogen KM to ' num2str(abs(v.km_gly)) ' mmol']}];
end
if v.km_tre > 0
    model.km(stridx('trehalose exchange' ,model.rxnNames)) = abs(v.km_tre);
    model.modifications = [model.modifications {['Changed trehalose KM to ' num2str(abs(v.km_tre)) ' mmol']}];
end

model.S = sparse(model.S);

%add a fied to track changes
model.modifications = [model.modifications {'Made S matrix sparse'}];

%% add exchange of metabolites necessary for anaerobic growth
%the exchange reactions for episterol, ergosterol, fecosterol, lanosterol, zymosterol and 
%ergosta-5,7,22,24(28)-tetraen-3beta-ol must have nonzero lower bounds
%cite Yeast 6 paper https://doi.org/10.1093/database/bat059
% ex.episterol = stridx('episterol exchange',model.rxnNames,false);
% ex.ergosterol = stridx('ergosterol exchange',model.rxnNames,false);
% ex.fecosterol = stridx('fecosterol exchange',model.rxnNames,false);
% ex.lanosterol = stridx('lanosterol exchange',model.rxnNames,false);
% ex.zymosterol = stridx('zymosterol exchange',model.rxnNames,false);
% ex.ergosta = stridx('ergosta-5,7,22,24(28)-tetraen-3beta-ol exchange',model.rxnNames,false);
% sterolrxns = [ex.episterol ex.ergosterol ex.fecosterol ex.lanosterol ex.zymosterol ex.ergosta];

%model.modifications = [model.modifications {['Change ''episterol exchange'' lb from ' num2str(model.lb(ex.episterol)) ' to ' defaultlb]}];
%model.modifications = [model.modifications {['Change ''ergosterol exchange'' lb from ' num2str(model.lb(ex.ergosterol)) ' to ' lb]}];
%model.modifications = [model.modifications {['Change ''fecosterol exchange'' lb from ' num2str(model.lb(ex.fecosterol)) ' to ' defaultlb]}];
%model.modifications = [model.modifications {['Change ''lanosterol exchange'' lb from ' num2str(model.lb(ex.lanosterol)) ' to ' defaultlb]}];
%model.modifications = [model.modifications {['Change ''zymosterol exchange'' lb from ' num2str(model.lb(ex.zymosterol)) ' to ' defaultlb]}];
%model.modifications = [model.modifications {['Change ''ergosta-5,7,22,24(28)-tetraen-3beta-ol exchange'' lb from ' num2str(model.lb(ex.ergosta)) ' to ' defaultlb]}];
%model.lb(sterolrxns) = defaultlb;
%model.lb(ex.ergosterol) = lb;
% model.lb(ex.zymosterol) = lb;
% model.lb(ex.episterol) = lb;
% model.lb(ex.fecosterol) = lb;
% model.lb(ex.lanosterol) = lb;
%model.lb(ex.ergosta) = lb;
% 
% model.lb(sterolrxns(1)) = 0;
% model.lb(sterolrxns(2)) = 0;
% model.lb(sterolrxns(3)) = 0;
% model.lb(sterolrxns(4)) = 0;
% model.lb(sterolrxns(5)) = 0;
% model.lb(sterolrxns(6)) = 0;
%adding ergosterol + zymosterol + fecosterol gets very low growth


%model.modifications = [model.modifications {['Change ''iron(2+) exchange'' lb from ' num2str(model.lb(stridx('iron(2+) exchange',model.rxnNames))) ' to ' defaultlb]}];
% model.modifications = [model.modifications {['Change ''glutathione disulfide exchange'' lb from ' num2str(model.lb(stridx('glutathione disulfide exchange',model.rxnNames))) ' to ' num2str(lb)]}];
% %model.lb(stridx('iron(2+) exchange',model.rxnNames)) = defaultlb;
% model.lb(stridx('glutathione disulfide exchange',model.rxnNames)) = lb;

% oxy = stridx('oxygen exchange',model.rxnNames);
% model.lb(oxy) = 0;

%% Change metabolite names to be consistent with exchanges and enzyme production requirements

model.mets{268} = 'ac[e]';
model.mets{288} = 'adp[c]';
model.mets{311} = 'nh4[e]';
model.mets{324} = 'atp[c]';
model.mets{345} = 'co2[e]';
model.mets{431} = 'glc-D[e]';
model.mets{437} = 'man-D[e]';
model.mets{444} = 'xyl-D[e]';
model.mets{512} = 'ergst[e]';
model.mets{520} = 'etoh[e]';
model.mets{556} = 'gdp[c]';
model.mets{572} = 'gthox[e]';
model.mets{582} = 'glycogen[c]';
model.mets{592} = 'gtp[c]';
model.mets{603} = 'h[e]';
model.mets{612} = 'h2o[e]';
model.mets{720} = 'fe2[e]';
model.mets{748} = 'ala-L[c]';
model.mets{758} = 'arg-L[c]';
model.mets{762} = 'asn-L[c]';
model.mets{766} = 'asp-L[c]';
model.mets{774} = 'cys-L[c]';
model.mets{784} = 'glu-L[c]';
model.mets{791} = 'gln-L[c]';
model.mets{795} = 'gly[c]';
model.mets{798} = 'his-L[c]';
model.mets{807} = 'ile-L[c]';
model.mets{812} = 'leu-L[c]';
model.mets{816} = 'lys-L[c]';
model.mets{820} = 'met-L[c]';
model.mets{823} = 'phe-L[c]';
model.mets{826} = 'pro-L[c]';
model.mets{830} = 'ser-L[c]';
model.mets{836} = 'thr-L[c]';
model.mets{839} = 'trp-L[c]';
model.mets{842} = 'tyr-L[c]';
model.mets{847} = 'val-L[c]';
model.mets{1005} = 'o2[e]';
model.mets{1037} = 'pi[e]';
model.mets{1055} = 'k[e]';
model.mets{1135} = 'so4[e]';
model.mets{1168} = 'tre[e]'; %trehalose is used as an internal carbon store
model.mets{1215} = 'zymst[e]';
model.mets{stridx('Cu2(+) [extracellular]',model.metNames,false)} = 'cu2[e]';
model.mets{stridx('Mn(2+) [extracellular]',model.metNames,false)} = 'mn2[e]';
model.mets{stridx('Zn(2+) [extracellular]',model.metNames,false)} = 'zn2[e]';
model.mets{stridx('Mg(2+) [extracellular]',model.metNames,false)} = 'mg2[e]';
model.mets{stridx('Ca(2+) [extracellular]',model.metNames,false)} = 'ca2[e]';
model.mets{stridx('chloride [extracellular]',model.metNames,false)} = 'cl[e]';
model.mets{stridx('sodium [extracellular]',model.metNames,false)} = 'na1[e]';
model.mets{stridx('oleate [extracellular]',model.metNames,false)} = 'oleate[e]';
model.mets{stridx('palmitoleate [extracellular]',model.metNames,false)} = 'palmitoleate[e]';

model.modifications = [model.modifications {'Replaced the following ids in mets field: ac[e] adp[c] nh4[e] atp[c] co2[e] glc-D[e] ergst[e] etoh[e] gdp[c] gthox[e] gtp[c] h[e] h2o[e] fe2[e] man-D[e] ala-L[c] arg-L[c] asn-L[c] asp-L[c] cys-L[c] glu-L[c] gln-L[c] gly[c] his-L[c] ile-L[c] leu-L[c] lys-L[c] met-L[c] phe-L[c] pro-L[c] ser-L[c] thr-L[c] trp-L[c] tyr-L[c] val-L[c] o2[e] pi[e] k[e] so4[e] tre[e] zymst[e] oleate[e] palmitoleate[e] mn2[e] zn2[e] mg2[e] ca2[e] cl[e] na1[e]'}];

%% Add the enzyme reactions
%NOTE: ATP/GTP costs not included here!

trna_bound_names = {'Ala-tRNA(Ala) [cytoplasm]' 'Arg-tRNA(Arg) [cytoplasm]' ...
    'Asn-tRNA(Asn) [cytoplasm]' 'Asp-tRNA(Asp) [cytoplasm]' 'Cys-tRNA(Cys) [cytoplasm]' 'Gln-tRNA(Gln) [cytoplasm]' ...
    'Glu-tRNA(Glu) [cytoplasm]' 'Gly-tRNA(Gly) [cytoplasm]' 'His-tRNA(His) [cytoplasm]' 'Ile-tRNA(Ile) [cytoplasm]' ...
    'Leu-tRNA(Leu) [cytoplasm]' 'Lys-tRNA(Lys) [cytoplasm]' 'Met-tRNA(Met) [cytoplasm]' 'Phe-tRNA(Phe) [cytoplasm]' ...
    'Pro-tRNA(Pro) [cytoplasm]' 'Ser-tRNA(Ser) [cytoplasm]' 'Thr-tRNA(Thr) [cytoplasm]' 'Trp-tRNA(Trp) [cytoplasm]' ...
    'Tyr-tRNA(Tyr) [cytoplasm]' 'Val-tRNA(Val) [cytoplasm]'};
trna_bound_idx = zeros(20,1);
for i = 1:length(trna_bound_names)
    trna_bound_idx(i) = stridx(trna_bound_names{i},model.metNames,false);
    model.mets(trna_bound_idx(i)) = trna_bound_names(i);
end
trna_free_names = {'tRNA(Ala) [cytoplasm]' 'tRNA(Arg) [cytoplasm]' 'tRNA(Asn) [cytoplasm]' ...
    'tRNA(Asp) [cytoplasm]' 'tRNA(Cys) [cytoplasm]' 'tRNA(Gln) [cytoplasm]' 'tRNA(Glu) [cytoplasm]' ...
    'tRNA(Gly) [cytoplasm]' 'tRNA(His) [cytoplasm]' 'tRNA(Ile) [cytoplasm]' 'tRNA(Leu) [cytoplasm]' ...
    'tRNA(Lys) [cytoplasm]' 'tRNA(Met) [cytoplasm]' 'tRNA(Phe) [cytoplasm]' 'tRNA(Pro) [cytoplasm]' ...
    'tRNA(Ser) [cytoplasm]' 'tRNA(Thr) [cytoplasm]' 'tRNA(Trp) [cytoplasm]' 'tRNA(Tyr) [cytoplasm]' ...
    'tRNA(Val) [cytoplasm]'};
trna_free_idx = zeros(20,1);
for i = 1:length(trna_free_names)
    trna_free_idx(i) = stridx(trna_free_names{i},model.metNames,false);
    model.mets(trna_free_idx(i)) = trna_free_names(i);
end

%stoich = stoich / seqlength;
fullstoich = [-v.stoich v.stoich -v.atpstoich v.atpstoich -v.gtpstoich v.gtpstoich 1/v.costfactor];
mets = [model.mets(trna_bound_idx); model.mets(trna_free_idx); 'atp[c]'; 'adp[c]'; 'gtp[c]'; 'gdp[c]'; 'enzyme[c]'];
model = addReaction(model,'produce_enzyme','metaboliteList',mets,'stoichCoeffList',fullstoich,'reversible',false);
model.modifications = [model.modifications {'Added reaction produce_enzyme'}];

%% Add the enzyme flow reactions to set alpha
model = addExchangeRxn(model,{'enzyme[e]'},0,1000);
model.modifications = [model.modifications {'Added exchange reaction for enzyme[e]'}];
objidx = stridx(v.objname,model.rxnNames,false);
objabbrevname = model.rxns(objidx);
objS = model.S(:,objidx);
objmets = find(objS);
objstoich = objS(objmets);
model = changeRxnMets(model,model.mets(objmets),[model.mets(objmets) 'enzyme[c]' 'enzyme[e]'],objabbrevname,[objstoich; -v.alpha; v.alpha]);
model.modifications = [model.modifications {['Added enzyme secretion to objective reaction ' v.objname]}];
model = refreshCsense(model);

%% Change the way maintenance flux is handled
if v.addEnzymeToMaint
    model = addEnzymeToMaintRxn(model,v.alpha);
end
if v.maintRxnAsObjective
    model = maintenanceAsObjective(model);
end

%%
model.v = v;