function [stoich,metnames] = getEnzymeStoich(v)
%GETENZYMECOST Get the stoichiometry and metabolite names for enzyme
%production. Units of enzyme produced are grams/GDW
%   v must contain the following fields (defaults in parens): costfactor(1)

%NOTE: Don't use this! The resulting tiny fluxes get rounded off to 0 
% if isfield(v,'enzByWeight')
%     unitsAreWeight = v.enzByWeight;
% else
%     unitsAreWeight = false;
% end

if ~isfield(v,'costfactor')
    v.costfactor = 1;
end


% aaseq = 'MFRTAALVSFVFLATASAQQVGTSAAETHPALSYQTCTSSGCTNVAGSVTLDSNWRWVHNVGGYTNCYTGNQWDATLCPDPATCAKNCALDGADYSGTYGITSSGSALTLAFKTGSNVGSRVYLMKDDSTYETFNLLNKEFTMDVDVSKLPCGLNGAVYFSAMSANGDASSTNAAGAKYGTGYCDSQCPHDLKFIQGQANIIGWNGTSANSGTGSIGACCSEMDIWEANSISEAVTPHPCTGTGITACTGSTCTSTCDQAGCDFNSYRMGNTTFYGPGMTVDTTKPITVVTQFLTSDNTTTGSLSEIRRFYVQNGVVIPNSQSTIAGVSGNSITDSWCAAQKTAFGDTNTFQSLGGLKQLGAAFAKGMVLVLSIWDDYAVNMLWLDSDYPTNKDASAPGVKRGTCATTSGVPADVEANGGSISVTYSNIKVGDIGTTFKAGQSSGSGSSGSSGSSTSTASGSGSTGSSGSVAKYGQCGGTGYAGATTCASGSTCTVLNPYYSQCL';
% %Punctularia strigosozonata CBHI(1) http://www.uniprot.org/uniprot/R7S4L5
% %Note this is not the exact sequence Melisa's using

%T reesei endoglucanase EGI, Cel7b
%Source: https://genome.jp/dbget-bin/www_bget?tre:TRIREDRAFT_123989
%We use Cel7b because that's the one that has mass data on BRENDA
aaseq_EGI = 'MAPSVTLPLTTAILAIARLVAAQQPGTSTPEVHPKLTTYKCTKSGGCVAQDTSVVLDWNYRWMHDANYNSCTVNGGVNTTLCPDEATCGKNCFIEGVDYAASGVTTSGSSLTMNQYMPSSSGGYSSVSPRLYLLDSDGEYVMLKLNGQELSFDVDLSALPCGENGSLYLSQMDENGGANQYNTAGANYGSGYCDAQCPVQTWRNGTLNTSHQGFCCNEMDILEGNSRANALTPHSCTATACDSAGCGFNPYGSGYKSYYGPGDTVDTSKTFTIITQFNTDNGSPSGNLVSITRKYQQNGVDIPSAQPGGDTISSCPSASAYGGLATMGKALSSGMVLVFSIWNDNSQYMNWLDSGNAGPCSSTEGNPSNILANNPNTHVVFSNIRWGDIGSTTNSTAPPPPPASSTTFSTTRRSSTTSSSPSCTQTHWGQCGGIGYSGCKTCTSGTTCQYSNDYYSQCL';

%Saccharomycopsis fibuligeria beta-glucosidase
%Ref: Machida M.,Ohtsuki I.,Fukui S.,Yamashita I.	Nucleotide sequences of Saccharomycopsis fibuligera genes for extracellular beta-glucosidases as expressed in Saccharomyces cerevisiae.	Appl. Environ. Microbiol.	54	3147-3155	1988	3146949
aaseq_BGL1 = 'MLMIVQLLVFALGLAVAVPIQNYTQSPSQRDESSQWVSPHYYPTPQGGRLQDVWQEAYARAKAIVGQMTIVEKVNLTTGTGWQLDPCVGNTGSVPRFGIPNLCLQDGPLGVRFADFVTGYPSGLATGATFNKDLFLQRGQALGHEFNSKGVHIALGPAVGPLGVKARGGRNFEAFGSDPYLQGTAAAATIKGLQENNVMACVKHFIGNEQEKYRQPDDINPATNQTTKEAISANIPDRAMHALYLWPFADSVRAGVGSVMCSYNRVNNTYACENSYMMNHLLKEELGFQGFVVSDWGAQLSGVYSAISGLDMSMPGEVYGGWNTGTSFWGQNLTKAIYNETVPIERLDDMATRILAALYATNSFPTEDHLPNFSSWTTKEYGNKYYADNTTEIVKVNYNVDPSNDFTEDTALKVAEESIVLLKNENNTLPISPEKAKRLLLSGIAAGPDPIGYQCEDQSCTNGALFQGWGSGSVGSPKYQVTPFEEISYLARKNKMQFDYIRESYDLAQVTKVASDAHLSIVVVSAASGEGYITVDGNQGDRKNLTLWNNGDKLIETVAENCANTVVVVTSTGQINFEGFADHPNVTAIVWAGPLGDRSGTAIANILFGKANPSGHLPFTIAKTDDDYIPIETYSPSSGEPEDNHLVENDLLVDYRYFEEKNIEPRYAFGYGLSYNEYEVSNAKVSAAKKVDEELPEPATYLSEFSYQNAKDSKNPSDAFAPADLNRVNEYLYPYLDSNVTLKDGNYEYPDGYSTEQRTTPNQPGGGLGGNDALWEVAYNSTDKFVPQGNSTDKFVPQLYLKHPEDGKFETPIQLRGFEKVELSPGEKKTVDLRLLRRDLSVWDTTRQSWIVESGTYEALIGVAVNDIKTSVLFTI';
aaseq = [aaseq_EGI aaseq_BGL1];
seqlength = length(aaseq);

%weight = v.enz_weight; %8.5853183e-20; %grams

% aminos = {'ala__L_c' 'arg__L_c' 'asn__L_c' 'asp__L_c' 'cys__L_c' 'gln__L_c' 'glu__L_c'...
%     'gly_c' 'his__L_c' 'ile__L_c' 'leu__L_c' 'lys__L_c' 'met__L_c' 'phe__L_c' 'pro__L_c' ...
%     'ser__L_c' 'thr__L_c' 'trp__L_c' 'tyr__L_c' 'val__L_c'};
aminos = {'ala-L[c]' 'arg-L[c]' 'asn-L[c]' 'asp-L[c]' 'cys-L[c]' 'gln-L[c]' 'glu-L[c]'...
    'gly[c]' 'his-L[c]' 'ile-L[c]' 'leu-L[c]' 'lys-L[c]' 'met-L[c]' 'phe-L[c]' 'pro-L[c]' ...
    'ser-L[c]' 'thr-L[c]' 'trp-L[c]' 'tyr-L[c]' 'val-L[c]'};
% stoich = zeros(length(aminos)+1,1);
% 
stoich = zeros(26,1);
stoich(1) = count(aaseq,'A');
stoich(2) = count(aaseq,'R');
stoich(3) = count(aaseq,'N');
stoich(4) = count(aaseq,'D');
stoich(5) = count(aaseq,'C');
stoich(6) = count(aaseq,'Q');
stoich(7) = count(aaseq,'E');
stoich(8) = count(aaseq,'G');
stoich(9) = count(aaseq,'H');
stoich(10) = count(aaseq,'I');
stoich(11) = count(aaseq,'L');
stoich(12) = count(aaseq,'K');
stoich(13) = count(aaseq,'M');
stoich(14) = count(aaseq,'F');
stoich(15) = count(aaseq,'P');
stoich(16) = count(aaseq,'S');
stoich(17) = count(aaseq,'T');
stoich(18) = count(aaseq,'W');
stoich(19) = count(aaseq,'Y');
stoich(20) = count(aaseq,'V');

%Composition here generated by looking at the average enzymes across three
%species. See the file "s4.2_AA_comp_species" for more details. This is the
%"totalcounts" variable.



if v.n_enzymes == 3
    %This array takes all three enzymes into account, so it should be used with
    %costfactor=1.
    stoich = [157.4903; 49.5298; 91.4945; 97.2040; 34.4946; 60.4280; 69.3807; 163.1554; 24.0718; 69.1210; 118.2207; 53.2529; 28.5820; 59.3903; 82.0254; 151.1506; 138.9847; 35.7296; 71.2610; 102.3696];
    seqlength = 1657.4;
end
metnames = aminos;

stoich = stoich * -1;

%if (unitsAreWeight)
%    stoich(21) = 1 * weight / v.costfactor; %units for enzyme amount are in milligrams
%else
stoich(21) = 1/v.costfactor; %units for enzyme are number of molecules
%end
metnames = aminos;
metnames{21} = 'enzyme[c]';

%include energy costs
metnames{22} = 'atp[c]';
metnames{23} = 'adp[c]';
metnames{24} = 'gtp[c]';
metnames{25} = 'gdp[c]';
metnames{26} = 'pi[c]';
atpcost = seqlength * v.atp_per_peptide;
gtpcost = seqlength * v.gtp_per_peptide;
stoich = [stoich; -atpcost; atpcost; -gtpcost; gtpcost; (atpcost+gtpcost)];

end

