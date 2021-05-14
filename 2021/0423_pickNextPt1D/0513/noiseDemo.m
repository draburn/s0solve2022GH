clear
%
%setprngstates(0);
%setprngstates
%setprngstates(92361344)
setprngstates(99983264)
%
numTerms = 2+round( 10.0 * abs(randn) )
xValsA = linspace(-1e2,1e2,10001);
xValsB = linspace(-1e0,1e0,10001);
xValsC = linspace(-1e-2,1e-2,10001);
fValsA = zeros(size(xValsA));
fValsB = zeros(size(xValsB));
fValsC = zeros(size(xValsC));
%
for n=1:numTerms
	c0 = randn;
	phi0 = 2*pi*rand;
	omega0 = randn*exp(3.0*randn)
	c1 = randn;
	phi1 = 2*pi*rand
	omega1 = 100.0*randn*exp(3.0*abs(randn))
	%
	omega0 /= c0;
	c1 /= (c0*omega1);
	funchy_term = @(x)( c0 * cos( phi0 + omega0*x + c1*cos( phi1 + omega1*x ) ) );
	%
	fValsA += funchy_term(xValsA);
	fValsB += funchy_term(xValsB);
	fValsC += funchy_term(xValsC);
end
%
numFigs = 0;
%
numFigs++; figure(numFigs);
plot( xValsA, fValsA, 'o-' );
grid on;
%
numFigs++; figure(numFigs);
plot( xValsB, fValsB, 'o-' );
grid on;
%
numFigs++; figure(numFigs);
plot( xValsC, fValsC, 'o-' );
grid on;