function [ omega, bigA, bigB, bigC ] = extFinder_getOmega( bigS, bigP, xVals, gVals, nOfPtWiseMin, prm=[] );
	thisFile = "extFinder_getOmega";
	%
	if ( mygetfield( prm, "beFast", false ) )
		n = nOfPtWiseMin;
		%
		xOfPtWiseMin = xVals(n);
		vecX = xVals(n-1:n+1)';
		vecG = gVals(n-1:n+1)';
		bigDelta = vecX(3) - vecX(1);
		bigG0 = vecG(2);
		bigG1 = vecG(1) + vecG(3) - 2.0*vecG(2);
		%
		vecD = ( vecX - bigS ) / bigDelta;
		vecY = ( vecG - bigG0 ) / bigG1;
		matM = [ abs(vecD).^bigP, vecD, ones(3,1) ];
		vecCoeff = matM \ vecY;
		%
		bigA = vecCoeff(1);
		bigB = vecCoeff(2);
		bigC = vecCoeff(3);
		%
		funchGModel = @(x)( bigG0 + bigG1*( ...
		  bigA*abs((x-bigS)/bigDelta).^bigP + bigB*(x-bigS)/bigDelta + bigC ) );
		xExt = bigS - bigDelta*mypow( bigB/(bigA*bigP), 1.0/(bigP-1.0) );
		gExt = funchGModel( xExt );
		%
		maskVals = logical(ones(size(xVals)));
		maskVals(nOfPtWiseMin-1:nOfPtWiseMin+1) = 0;
		wVals = mygetfield( prm, "wVals", (1.0/sum(1.0./sqrt(gVals(maskVals))))./sqrt(gVals) );
		rhoWVals = wVals.*(funchGModel(xVals)-gVals).^2;
		omegaFit = 0.5 * sum( rhoWVals(maskVals) );
		%
		wSus = mygetfield( prm, "wSus", 1.0 );
		omegaSus = ( gExt < 0.0 ) * 0.5 * wSus * ( gExt / bigG1 )^2;
		%
		epsA = mygetfield( prm, "epsA", eps^0.75 );
		epsB = mygetfield( prm, "epsB", eps^0.50 );
		epsC = mygetfield( prm, "epsC", eps^0.75 );
		epsS = mygetfield( prm, "epsS", eps^0.75 );
		epsP = mygetfield( prm, "epsP", eps^0.75 );
		omegaRegu = 0.5 * ( epsA*bigA^2 + epsB*bigB^2 + epsC*bigC^2 + ...
		  epsS*(bigS-xOfPtWiseMin)^2 + epsP*(bigP-2.0)^2 );
		%
		omega = omegaFit * ( 1.0 + omegaSus ) + omegaRegu;
		%
		return;
	end
	%
	assert( isrealscalar(bigS) );
	assert( isrealscalar(bigP) );
	numPts = size(xVals,2);
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(gVals,[1,numPts]) );
	assert( isrealscalar(nOfPtWiseMin) );
	assert( fleq(nOfPtWiseMin,round(nOfPtWiseMin)) );
	assert( 1 <= nOfPtWiseMin );
	assert( nOfPtWiseMin <= numPts );
	n = nOfPtWiseMin;
	allGValsArePositive = (0==sum( 0.0 > gVals ) );
	assert( allGValsArePositive );
	isPtWiseMin = ( (gVals(n+1)>=gVals(n)) && (gVals(n-1)>=gVals(n)) );
	assert( isPtWiseMin );
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	[ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, gVals, nOfPtWiseMin, prm );
	assert( isrealscalar(bigA) );
	assert( isrealscalar(bigB) );
	assert( isrealscalar(bigC) );
	%
	n = nOfPtWiseMin;
	xOfPtWiseMin = xVals(n);
	bigDelta = xVals(n+1) - xVals(n-1);
	bigG0 = gVals(n);
	bigG1 = gVals(n+1) + gVals(n-1) - 2.0*gVals(n);
	funchGModel = @(x)( bigG0 + bigG1*( ...
	  bigA*abs((x-bigS)/bigDelta).^bigP + bigB*(x-bigS)/bigDelta + bigC ) );
	xExt = bigS - bigDelta*mypow( bigB/(bigA*bigP), 1.0/(bigP-1.0) );
	gExt = funchGModel( xExt );
	%
	maskVals = logical(ones(size(xVals)));
	maskVals(nOfPtWiseMin-1:nOfPtWiseMin+1) = 0;
	wVals = mygetfield( prm, "wVals", (1.0/sum(1.0./sqrt(gVals(maskVals))))./sqrt(gVals) );
	assert( isrealarray(wVals,[1,numPts]) );
	noWValuesAreNegative = ( 0 == sum(wVals<0.0) );
	assert( noWValuesAreNegative );
	rhoWVals = wVals.*(funchGModel(xVals)-gVals).^2;
	omegaFit = 0.5 * sum( rhoWVals(maskVals) );
	%
	wSus = mygetfield( prm, "wSus", 1.0 );
	assert( isrealscalar(wSus) );
	assert( wSus >= 0.0 );
	omegaSus = ( gExt < 0.0 ) * 0.5 * wSus * ( gExt / bigG1 )^2;
	%
	epsA = mygetfield( prm, "epsA", eps^0.75 );
	assert( isrealscalar( epsA ) );
	assert( 0.0 <= epsA );
	epsB = mygetfield( prm, "epsB", eps^0.50 );
	assert( isrealscalar( epsB ) );
	assert( 0.0 < epsB );
	epsC = mygetfield( prm, "epsC", eps^0.75 );
	assert( isrealscalar( epsC ) );
	assert( 0.0 <= epsC );
	epsS = mygetfield( prm, "epsS", eps^0.75 );
	assert( isrealscalar( epsS ) );
	assert( 0.0 <= epsS );
	epsP = mygetfield( prm, "epsP", eps^0.75 );
	assert( isrealscalar( epsP ) );
	assert( 0.0 <= epsP );
	omegaRegu = 0.5 * ( epsA*bigA^2 + epsB*bigB^2 + epsC*bigC^2 + ...
	  epsS*(bigS-xOfPtWiseMin)^2 + epsP*(bigP-2.0)^2 );
	%
	omega = omegaFit * ( 1.0 + omegaSus ) + omegaRegu;
	%
	%
	%
return;
end
