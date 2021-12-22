clear;
%
% integral -inf +inf du u * exp(-(a*u^2+b*u+c))
%  = (-0.5*b/a) * exp( b^2/(4*a) - c ) * sqrt(pi/a)
numPts = 1 + 1E6

a = abs(randn)
b = randn
c = randn
%
uMax = 20.0*(1.0+sqrt(abs(b/a))+abs(c/a));
tempVals0 = linspace(-1.0,1.0,numPts);
tempVals1 = abs(tempVals0).^4;
tempVals2 = tempVals1.*sign(tempVals0);
uEdges = uMax*tempVals2;
%
duVals = uEdges(2:end)-uEdges(1:end-1);
uVals = (uEdges(2:end)+uEdges(1:end-1))/2.0;
%
eVals = exp(-(  a*(uVals.^2) + b*uVals + c ));
%
i0_numerical = sum(eVals.*duVals)
i0_theoretical = exp( b^2/(4*a) - c ) * sqrt(pi/a)
%
i1_numerical = sum(uVals.*eVals.*duVals)
i1_theoretical = exp( b^2/(4*a) - c ) * sqrt(pi/a) * (-0.5*b/a)
%
i2_numerical = sum(uVals.*uVals.*eVals.*duVals)
i2_theoretical = exp( b^2/(4*a) - c ) * sqrt(pi/a) ...
 * ( (0.5*b/a)^2 + 0.5/a )
%
epsval = sqrt(eps);
i0_res = (epsval+abs(i0_numerical))/(epsval+abs(i0_theoretical)) - 1.0
i1_res = (epsval+abs(i1_numerical))/(epsval+abs(i1_theoretical)) - 1.0
i2_res = (epsval+abs(i2_numerical))/(epsval+abs(i2_theoretical)) - 1.0
