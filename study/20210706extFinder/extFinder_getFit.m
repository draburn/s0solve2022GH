function [ bigA, bigB, bigC ] = extFinder_getFit( bigS, bigP, xVals, gVals, nOfPtWiseMin, prm=[] );
	thisFile = "extFinder_getFit";
	if ( mygetfield( prm, "beFast", false ) )
		n = nOfPtWiseMin;
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
		return;
	end
	%
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
	%
return;
end
