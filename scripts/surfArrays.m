function [p,X,Y,Z] = surfArrays(x,y,z)
%Convert the three given vectors to a format that's compatable with the
%surf function
%Requires that |z| = |x| = |y|


ux = unique(x);
uy = unique(y);
[X,Y] = meshgrid(ux,uy);
Z = nan(length(ux),length(uy));
for i = 1:length(z)
    row = find(x(i) == ux);
    col = find(y(i) == uy);
    Z(row,col) = z(i);
end
Z = Z';
p = surf(X,Y,Z);     
end