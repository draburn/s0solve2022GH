clear;
setprngstates(0);
%
y0 = 1.0;
%
if (0)
	numMeas = 4;
	x0 = 1.0E-2
	x1 = 1.0E3
	%f0Norm = 1.0
	%f0Norm = 0.156663933871106
	f0Norm = 0.156664266831862
elseif (0)
	numMeas = 20;
	x0 = 1.0E-2
	x1 = 1.0E3
	%f0Norm = 1.0
	f0Norm = 1.88617943795974e-05
elseif (0)
	numMeas = 100;
	x0 = 1.0E-2
	x1 = 1.0E3
	%f0Norm = 1.0
	f0Norm = 3.45027525098795e-23;
else
	numMeas = 1000;
	x0 = 5.0E-1
	x1 = 1.0E2
	%f0Norm = 1.0
	%f0Norm = 2.05792070439881e-219
	f0Norm = 3.99698303844724e-219
end
%
funch_f0 = @(x)( exp( -0.5*numMeas*(y0^2)./(x.^2) ) ./ (x.^(numMeas)) );
funch_f1 = @(x)( funch_f0(x) / f0Norm );
%
if (0)
	funch_f = @(x)( funch_f1(x) );
	bigI_expected = 1.0
else
	funch_f = @(x)( (x.^2).*funch_f1(x) );
	bigI_expected = y0^2
end
%
numIntervals = 10000000
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
%%%plot( vec_x, vec_f, 'o-' );
%%%%semilogy( vec_x, vec_f, 'o-' );
%%%grid on;
echo__bigI = bigI
echo__bigIExpected = bigI_expected
echo__res = res
echo__resRel = res / bigI_expected
