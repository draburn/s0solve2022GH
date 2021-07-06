function vecCoeff = quadFit( xVals, fVals, prm=[], datIn=[] );
	thisFile = "quadFit";
	%
	%
	numPts = size(xVals,2);
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	%
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	wVals = mygetfield( prm, "wVals", ones(size(xVals)) );
	assert( isrealarray(wVals,[1,numPts]) );
	wValsAreAllNonNegative = (0==sum( 0.0 > wVals ));
	assert( wValsAreAllNonNegative );
	atLeastThreeWValsArePositive = (3<=sum( 0.0 < wVals ));
	assert( atLeastThreeWValsArePositive );
	%
	vecX = xVals';
	vecF = fVals';
	matW = diag(wVals);
	matX = [ vecX.^2, vecX, ones(numPts,1) ];
	vecCoeff = (matW*matX)\(matW*vecF);
	%
	%
return;
end
