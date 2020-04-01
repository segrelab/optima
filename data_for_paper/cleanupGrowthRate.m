function dat = cleanupGrowthRate(dat)
%Set all grwoth rates prior to the end of the lag phase to 0.
%This is intended to make sure the reported growth rate isn't due to the
%initial glucose inclusion
%The input dat should be a table that has been through the function
%populateColumnsForPaper1()
[nrows,ncols] = size(dat);
for i = 1:nrows
    row = dat(i,:);
    lag = row.lag_steps;
    for j = 1:length(lag)
        r = row.layout{1}.models{j}.growthrate;
        t_oldmax = find(r == max(r));
        t_oldmax = t_oldmax(1);
        r(lag(j):-1:1) = nan;
        dat.layout{i}.models{j}.growthrate;
        dat.growthrate_max(i,j) = max(r);
        t_newmax = find(r == max(r));
        t_newmax = t_newmax(1);
        
%         %diagnostic
%         if (t_oldmax == t_newmax)
%             warning(['Max rate timing matched: ' num2str(t_newmax)]);
%         else
%             warning(['New max rate timing: ' num2str(t_oldmax) ' -> ' num2str(t_newmax)]);
%         end
        
    end
end