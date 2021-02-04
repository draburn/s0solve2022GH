clear;
setprngstates(0);
%
numMeas = 100;
y0 = sqrt(3.0);
"Looks like this calculation isn't quite correct."
%bigI_expected = pi/2.0;
%%%f0Norm = 1.0
%%%f0Norm = 8.32437496310950e-47
%%%f0Norm = 8.32434618918507e-47;
f0Norm = 8.32437496310918e-47
%
%funch_f = @(x)( sqrt(1.0 - (x.*x)) );
funch_f0 = @(x)( exp( -0.5*numMeas*(y0^2)./(x.^2) ) ./ (x.^(numMeas)) );
funch_f1 = @(x)( funch_f0(x) / f0Norm );
%
%funch_f = @(x)( funch_f1(x) );
%bigI_expected = 1.0
%
funch_f = @(x)( (x.^2).*funch_f1(x) );
bigI_expected = y0^2
%
x0 = 0.01;
x1 = 5.0;
numIntervals = 100000
%
vec_x = linspace( x0, x1, numIntervals );
vec_f = funch_f( vec_x );
%
vec2_fAvg = ( vec_f(2:end) + vec_f(1:end-1) )/2.0;
vec2_xDif = ( vec_x(2:end) - vec_x(1:end-1) );
bigI = sum( vec2_fAvg .* vec2_xDif );
%
res = bigI - bigI_expected;
%
plot( vec_x, vec_f/max(vec_f), 'o-' );
echo__bigI = bigI
echo__bigIExpected = bigI_expected
echo__res = res
echo__resRel = res / bigI_expected
