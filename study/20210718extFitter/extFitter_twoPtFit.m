function [ omegaMesh, bigF0Mesh, bigFLMesh, bigFRMesh, datOut ] = extFitter_twoPtFit( ...
 sMesh, pLMesh, pRMesh, xVals, fVals, nL, nR, prm=[] );
	thisFile = "extFitter_twoPtFit";
	datOut = [];
	%
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	assert( isrealarray(sMesh,[size1,size2]) );
	assert( isrealarray(pLMesh,[size1,size2]) );
	assert( isrealarray(pRMesh,[size1,size2]) );
	numPts = size(xVals,2);
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	assert( isrealscalar(nL) );
	assert( isrealscalar(nR) );
	assert( fleq(nL,round(nL)) );
	assert( fleq(nR,round(nR)) );
	assert( 1 <= nL );
	assert( nL < nR );
	assert( nR <= numPts );
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	xL = xVals(nL);
	fL = fVals(nL);
	xR = xVals(nR);
	fR = fVals(nR);
	%bigX0 = xL;
	%bigX1 = xR - xL; %Yeah???
	bigX0 = 0.0;
	bigX1 = 1.0;
	%
	datOut.bigX0 = bigX0;
	datOut.bigX1 = bigX1;
	%
	yVals = (xVals-bigX0)/bigX1;
	yL = (xL-bigX0)/bigX1;
	yR = (xR-bigX0)/bigX1;
	wVals_default = 1.0./( abs(fVals) + eps*max(abs(fVals)) );
	wVals_default /= sum(wVals_default);
	wVals = mygetfield( prm, "wVals", wVals_default );
	assert( isrealarray(wVals,[1,numPts]) );
	%
	parfor i1=1:size1
	parfor i2=1:size2
		s = sMesh(i1,i2);
		pL = pRMesh(i1,i2);
		pR = pRMesh(i1,i2);
		%
		rR = abs( yR - s ).^pR;
		lL = abs( yL - s ).^pL;
		lVals = (yVals<=yL).*abs( yVals - s ).^pL;
		rVals = (yVals>=yR).*abs( yVals - s ).^pR;
		cVals = 1.0 - lVals/lL - rVals/rR;
		dVals = fVals - fL*lVals/lL - fR*rVals/rR;
		bigF0 = sum( wVals .* cVals .* dVals ) / sum( wVals .* cVals .* cVals );
		bigFL = ( fL - bigF0 ) / lL;
		bigFR = ( fR - bigF0 ) / rR;
		omega = 0.5 * sum( wVals .* ( bigF0*cVals - dVals ).^2 );
		%
		bigF0Mesh(i1,i2) = bigF0;
		bigFLMesh(i1,i2) = bigFL;
		bigFRMesh(i1,i2) = bigFR;
		omegaMesh(i1,i2) = omega;
	end
	end
return;
end
