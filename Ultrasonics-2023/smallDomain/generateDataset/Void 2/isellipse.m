% ISELLIPSE    This program checks whether a point (x,y) lies inside,outside
%            or on a ellipse defined by centre at (h,k) and major axis a,
%            minor axis b. and with an angle A with x axis.
%
%   Syntax: iscircle(h,k,a,b,A)
%
%           Program checks whether point (x,y) lies inside,outside or on the ellipse.
%           ans = 1  ==> lie on the ellipse.
%           ans > 1  ==> lie outside the ellipse.
%           ans < 1 ==> lie inside the ellipse.

function result=isellipse(H,K,x,y,a,b,A)



val=((x-H)*cosd(A)+ (y-K)*sind(A))^2/a^2 + ((x-H)*sind(A)- (y-K)*cosd(A))^2/b^2;

result=(val);
