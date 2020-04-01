function [rss,r2,V,K,rmse] = findFitQuality(timepoints_exp,expdata,timepoints_sim,simdata,varargin)
%FINDFITQUALITY Compares the fit of the given experimental and simulated data
%Optional arguments:
%   order: max power to which to raise the predictor
%   mode: if 'log10' scale the data to log scale before finding RSS
order = 3;
%mode = 'normal';
mode = 'log10';

timepoints_exp = timepoints_exp(1:8);
expdata = expdata(1:8);

if nargin > 4
    order = varargin{1};
end
if nargin >5
    mode = varargin{2};
end

if strcmp(mode,'log10') %change the data to log scale
    expdata = log10(expdata);
    simdata = log10(simdata);
end

%Step 1: Find the equation that describes the simulated data
X = ones(size(timepoints_exp'));
for i = 1:order
    X = horzcat(X, power(timepoints_exp',i));
end
timepoints = ismember(timepoints_sim,timepoints_exp);
timepoints = find(timepoints);
sdt = simdata(timepoints);
%X = X
try 
[b,bint,r,rint,stats] = regress(sdt,X);
catch ME
    errorobj = ME;
    rethrow(ME);
end

%stats contains the R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
yfit = b(1);
for i = 1:order
    yfit = yfit + (power(timepoints_exp,i) * b(i+1));
end


%Step 2: Find the likelihood of the experimental data arising from the sim
%equation
rss = sum(power(expdata-yfit,2));

%Root mean square error. RMSE = 0 for a perfect fit
rmse = sqrt(rss/length(expdata));

%R-squared = 1 - residual sum of squares / total sum of squares
sstot = sum(power(expdata-mean(expdata),2));
r2 = 1 - (rss/sstot);

%Kuiper's test with the null hypothesis that the errors are
%normally distributed
%[phat] = mle(mean(r)); %maximum likelihood estimation for a normal distribution
%phat(1) is the MLE for the mean, phat(2) is the MLE for the standard
%deviation
mad = abs(mean(r)); %Mean Absolute Difference from the assumed mean, 0
%Kuiper's test: Define: Samples X are drawn from cumulative distribution function F
%       z_i = F(x_i)
%       D+ = max[i/n - z_i]
%       D- = max[z_i - (i - 1)/n]
%       V = D+ + D-
%   Smaller V = better fit
n = length(r);
rsort = sort(r);
dplus = 0;
dminus = 0;
for i = 1:n
    zi = normcdf(rsort(i),0,sqrt(mad));
    p = (i/n) - zi;
    if p > dplus
        dplus = p;
    end
    m = zi - ((i-1)/n);
    if m > dminus
        dminus = m;
    end
end
V = dplus + dminus;

%Kolmogorov-Smirnov test
K = 0;
for i = 1:n
    s = 0;
    for j = 1:i
        s = s + rsort(j);
    end
    fn = (s/n) - normcdf(rsort(i),0,sqrt(mad));
    K = max(K,fn);
end

