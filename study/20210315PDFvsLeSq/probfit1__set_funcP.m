numPts = size(xVals,2)
s1 = numPts
sxx = sum(xVals.*xVals)
syy = sum(yVals.*yVals)
sxt2 = sum(xVals) * 2.0
syt2 = sum(yVals) * 2.0
sxyt2 = sum(xVals.*yVals) * 2.0
%
funch_p = @(c0,c1,sr) (abs(sr).^numPts) .* exp(-(sr.^2).*( ...
   c0.*( s1*c0 + sxt2*c1 - syt2 ) ...
 + c1.*( sxx*c1 - sxyt2 )  + syy ));

