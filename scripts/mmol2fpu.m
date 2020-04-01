function fpu = mmol2fpu(mmol,kcat,km)
%MMOL2FPU Calculate the Filter Paper Units from the given cellulase
%concentration.
%   See https://www.nrel.gov/docs/gen/fy08/42628.pdf
%   The value of 2.0mg of reducing sugar as glucose from 50mg of filter
%   paper (4% conversion) in 60 minutes is designated as the intercept for
%   calculating FPU by IUPAC.

%find two points close to 2.0
dist = 100;
mindist = .01;
left = [0 0]; %coordinates are glucose released x enzyme dilution
right = [99 99];
con = mmol/2; %start with a 1:1 dilution
dilution = .5; %ml of original stock
while dist > mindist
    
    %g = enz * kcat * 3600 * (49/KM+49)
    g = con * kcat * 3600 * 49 / (km + 49);
    
    if g < 2
        if g > left(1)
            left = [g dilution];
        end
        dilution = 1.5 * dilution;
        con = 1.5 * con;
    else % g > 2
        if g < right(1)
            right = [g dilution];
        end
        dilution = dilution / 2;
        con = con / 2;        
    end 
    dist = sqrt((right(1)-left(1))^2 + (right(2)-left(2))^2);
end

%The point [2,y] lies on a straight line on a semilog plot between the two
%points found above
x = [left(1) right(1)];
y = [log10(left(2)) log10(right(2))];
p = polyfit(x,y,1);
d = p(1)*2 + p(2);
d = 10^d;
fpu = 0.37/d;
