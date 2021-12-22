function [ bigX, bigP, datOut ] = extFinder_dlogdf( xVals, fVals, prm=[], datIn=[] );
	commondefs;
	thisFile = "extFinder_dlogdf";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	datOut = [];
	%
	%
	%
	numPts = size(xVals,2);
	assert( numPts >= 5 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	%
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	fValsAreAllNonzero = (0==sum( 0.0 == fVals ));
	assert( fValsAreAllNonzero );
	%
	fValsAllHaveSameSign = (0==sum( 0.0 >= sign(fVals(1)) * fVals ));
	assert( fValsAllHaveSameSign );
	gVals = abs(fVals);
	%
	%
	%
	numPts_h = 0;
	xVals_h = [];
	hVals_h = [];
	for n=2:numPts-1
		x = xVals(n);
		p = polyfit( xVals(n-1:n+1), gVals(n-1:n+1), 2 );
		df = 2.0*p(1)*x + p(2);
		ddf = 2.0*p(1);
		if ( abs(ddf) > sqrt(eps)*abs(df) )
			numPts_h++;
			xVals_h(numPts_h) = x;
			hVals_h(numPts_h) = df/abs(ddf);
		end
	end
	assert( isrealarray(xVals_h,[1,numPts_h]) );
	assert( isrealarray(hVals_h,[1,numPts_h]) );
	hScale = sqrt(sum(hVals_h.^2));
	%
	% Linear fit, with weight...
	vecX = xVals_h';
	vecH = hVals_h';
	vecW = 1.0./( sqrt(abs(vecH)) + sqrt(eps*hScale) );
	matX = [ vecX, ones(numPts_h,1) ];
	matW = diag(vecW);
	vecCoeff = (matW*matX)\(matW*vecH);
	assert( 0.0 < vecCoeff(1) );
	%
	bigX = -vecCoeff(2)/vecCoeff(1);
	bigP = 1.0 + 1.0/vecCoeff(1);
	if ( mygetfield(prm,"makePlot",false) )
		plot( ...
		  xVals, gVals, 'o-', ...
		  bigX*[1,1], [min(gVals),max(gVals)], 'x-' );
		grid on;
		figure();
		plot( ...
		  xVals_h, hVals_h, 'o-', ...
		  xVals_h, vecCoeff(1)*xVals_h + vecCoeff(2), 'x-' );
		grid on;
	end
	%
thisFile = [ "RETURN from " thisFile ];	
return;
end
