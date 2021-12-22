clear;
setprngstates(0);
%
bigI_expected = pi/2.0;
funch_f = @(x)( sqrt(1.0 - (x.*x)) );
x0 = -1.0;
x1 = 1.0;
bigN = 1000;
%
vec_x = linspace( x0, x1, bigN );
vec_f = funch_f( vec_x );
%
vec2_fAvg = ( vec_f(2:end) + vec_f(1:end-1) )/2.0;
vec2_xDif = ( vec_x(2:end) - vec_x(1:end-1) );
bigI = sum( vec2_fAvg .* vec2_xDif );
%
res = bigI - bigI_expected;
%
echo__bigI = bigI
echo__bigIExpected = bigI_expected
eco__res = res