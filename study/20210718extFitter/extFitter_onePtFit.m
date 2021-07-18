function [ omegaMesh, bigF0Mesh, bigF1Mesh, datOut ] = extFitter_onePtFit( ...
 sMesh, pMesh, xVals, fVals, nFit, prm=[] );
	thisFile = "extFitter_onePtFit";
	datOut = [];
	%
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	assert( isrealarray(sMesh,[size1,size2]) );
	assert( isrealarray(pMesh,[size1,size2]) );
	numPts = size(xVals,2);
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	assert( isrealscalar(nFit) );
	assert( fleq(nFit,round(nFit)) );
	assert( 1 <= nFit );
	assert( nFit <= numPts );
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	xFit = xVals(nFit);
	fFit = fVals(nFit);
	%bigX0 = xFit;
	%if ( nFit > (numPts+1)/2 )
	%	% Do I really want to do this?
	%	bigX1 = xVals(end)-xVals(1);
	%else
	%	bigX1 = xVals(1)-xVals(end);
	%end
	bigX0 = 0.0;
	bigX1 = 1.0;
	%
	datOut.bigX0 = bigX0;
	datOut.bigX1 = bigX1;
	%
	yVals = (xVals-bigX0)/bigX1;
	yFit = (xFit-bigX0)/bigX1;
	wVals_default = 1.0./( abs(fVals) + eps*max(abs(fVals)) );
	wVals_default /= sum(wVals_default);
	wVals = mygetfield( prm, "wVals", wVals_default );
	assert( isrealarray(wVals,[1,numPts]) );
	%
	parfor i1=1:size1
	parfor i2=1:size2
		s = sMesh(i1,i2);
		p = pMesh(i1,i2);
		%
		gFit = abs( yFit -s ).^p;
		gVals = abs( yVals - s).^p;
		cVals = 1.0 - gVals/gFit;
		dVals = fVals - fFit*gVals/gFit;
		bigF0 = sum( wVals .* cVals .* dVals ) / sum( wVals .* cVals .* cVals );
		bigF1 = ( fFit - bigF0 ) / gFit;
		omega = 0.5 * sum( wVals .* ( bigF0*cVals - dVals ).^2 );
		%
		bigF0Mesh(i1,i2) = bigF0;
		bigF1Mesh(i1,i2) = bigF1;
		omegaMesh(i1,i2) = omega;
	end
	end
return;
end
