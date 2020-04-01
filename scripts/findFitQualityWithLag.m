function [rmse,shift] = findFitQualityWithLag(timepoints_exp,expdata,timepoints_sim,simdata,varargin)
%FINDFITQUALITY Compares the fit of the given experimental and simulated data
%Optional arguments:
%   order: max power to which to raise the predictor
%   mode: if 'log10' scale the data to log scale before finding RSS

%slide the sim plot left and right, and report the best fit and the
%corresponding shift
exp_end = max(timepoints_exp);
sim_end = max(timepoints_sim);
maxshift = abs(sim_end - exp_end);

%prepend a bunch of negative timepoints T<0 to the sim timepoints
pre = -maxshift:-1;
ts = [pre'; timepoints_sim];
%prepend pre-growth data to the sim data when T<0
predat = repmat(simdata(1),length(pre),1);
tsd = [predat; simdata];

rmse = 99999;
shift = 99999;
for s = -maxshift:maxshift
    %shift the sim timeline
    ts_shifted = ts + s;
    [rss,r2,V,K,temprmse] = findFitQuality(timepoints_exp,expdata,ts_shifted,tsd,varargin{:});
    if temprmse < rmse
        rmse = temprmse;
        shift = s;
    end
end