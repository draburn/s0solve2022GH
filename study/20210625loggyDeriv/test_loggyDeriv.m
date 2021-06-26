clear;
commondefs;
thisFile = "test_groot1d";
numFigs = 0;
tic();

setprngstates();
verbLev = VERBLEV__COPIOUS;
caseNum = 1;

prm = []; % Unless overwritten by case.
dat = []; % Unless overwritten by case.

switch (caseNum)
case (0)
	setprngstates(0);
	msg_main( verbLev, thisFile, __LINE__, ...
	  sprintf("Case: Simple test case.") );
	bigA_secret = 1.0;
	bigB_secret = 1.0;
	bigX_secret = 0.0;
	bigP_secret = 5.1;
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ) .^ bigP_secret );
	numPts = 21;
	xVals = bigX_secret + linspace( -2.0, 2.0, numPts );
case (1)
	setprngstates(0);
	msg_main( verbLev, thisFile, __LINE__, ...
	  sprintf("Case: Basic test case.") );
	bigA_secret = randn()
	bigB_secret = randn()
	bigX_secret = randn()
	bigP_secret = 2.0 + abs(randn())
	funchF = @(x)( bigA_secret + bigB_secret * abs( x - bigX_secret ) .^ bigP_secret );
	numPts = round( 1000 + abs(randn()*exp(randn())) );
	xVals = sort([ ...
	  bigX_secret-abs(randn(1,2)), ...
	  bigX_secret+abs(randn(1,2)), ...
	  bigX_secret+randn(1,numPts-4) ]);
otherwise
	error(sprintf("Unsupported value of caseNum (%d).",caseNum) );
end
fVals = funchF(xVals);

numFigs++; figure(numFigs);
plot( xVals, fVals, 'o-' );
grid on;

polyOrder = 2;
for nLo = 1:numPts-polyOrder
	nHi = nLo+polyOrder;
	vecX = xVals(nLo:nHi)';
	matX = ones(polyOrder+1,1);
	for m=1:polyOrder
		matX = [ matX, vecX.^m ];
	end
	vecF = fVals(nLo:nHi)';
	vecC = matX\vecF;
	%
	thisIndexTimesTwo = nLo+nHi;
	if ( 0==mod(thisIndexTimesTwo,2) )
		thisX = xVals(round(thisIndexTimesTwo/2.0));
	else
		thisX = 0.5*( xVals(floor(thisIndexTimesTwo/2.0)) + xVals(ceil(thisIndexTimesTwo/2.0)) );
	end
	%
	thisF = vecC(1);
	for m=1:polyOrder
		thisF += vecC(m+1) * (thisX.^m);
	end
	%
	thisFP = vecC(2);
	for m=2:polyOrder
		thisFP += vecC(m+1) * (m.*thisX.^(m-1));
	end
	%
	thisFPP = 2.0*vecC(3);
	for m=3:polyOrder
		thisFPP += vecC(m+1) * (m.*(m-1.0).*thisX.^(m-2));
	end
	%
	calcd_index = nLo;
	calcd_x(calcd_index) = thisX;
	calcd_f(calcd_index) = thisF;
	calcd_fp(calcd_index) = thisFP;
	calcd_fpp(calcd_index) = thisFPP;
end
calcd_fofp = calcd_f./calcd_fp;
calcd_fpof = calcd_fp./calcd_f;
calcd_fpofpp = calcd_fp./calcd_fpp;

numFigs++; figure(numFigs);
plot( calcd_x, calcd_fofp, 'o-' );
grid on;

numFigs++; figure(numFigs);
plot( calcd_x, calcd_fpof, 'o-' );
grid on;

numFigs++; figure(numFigs);
plot( calcd_x, calcd_fpofpp, 'o-' );
grid on;
