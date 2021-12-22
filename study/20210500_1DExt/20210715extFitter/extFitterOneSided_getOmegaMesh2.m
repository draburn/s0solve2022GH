function [ omegaMesh, bigF0Mesh, bigF1Mesh, datOut ] = extFitterOneSided_getOmegaMesh2( ...
  sMesh, pMesh, xVals, fVals, nOfClosestPt, prm=[] );
	thisFile = "extFitterOneSided_getOmegaMesh2";
	datOut = [];
	%
	%
	%
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	assert( isrealarray(sMesh,[size1,size2]) );
	assert( isrealarray(pMesh,[size1,size2]) );
	numPts = size(xVals,2);
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	assert( isrealscalar(nOfClosestPt) );
	assert( fleq(nOfClosestPt,round(nOfClosestPt)) );
	assert( 1 <= nOfClosestPt );
	assert( nOfClosestPt <= numPts );
	xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
	assert( xValsAreStrictlyIncreasing );
	%
	xc = xVals(nOfClosestPt);
	fc = fVals(nOfClosestPt);
	delta = max(xVals)-min(xVals);
	datOut.bigX0 = xc;
	datOut.bigX1 = delta;
	%
	yVals = (xVals-xc)/delta;
	%
	parfor i1=1:size1
	parfor i2=1:size2
		s = sMesh(i1,i2);
		p = pMesh(i1,i2);
		%
		hVals = abs( yVals - s ).^p- abs( s ).^p;
		fmfcVals = fVals - fc;
		%bigF1Mesh(i1,i2) = sum( fmfcVals .* hVals ) ./ sum( hVals .* hVals );
		%omegaMesh(i1,i2) = sum(( bigF1Mesh(i1,i2)*hVals - fmfcVals ).^2);
		%
		%wVals = 1.0./(abs(hVals)+sqrt(eps)*max(abs(hVals)));
		%wVals = 1.0./(abs(hVals)+eps*max(abs(hVals)));
		%wVals = wVals.^2.0;
		%wVals = 1.0./( abs(fVals)+sqrt(eps)*max(abs(fVals)) );
		wVals = 1.0./abs(fVals);
		%wVals = wVals.^-2.0;
		wVals /= sum(wVals);
		bigF1 = sum( wVals .* fmfcVals .* hVals ) ./ sum( wVals .* hVals .* hVals );
		bigF0 = fc - bigF1*( abs(s).^p );
		omega = 0.5*sum( wVals .* ( bigF1*hVals - fmfcVals ).^2);
		%
		bigF0Mesh(i1,i2) = bigF0;
		bigF1Mesh(i1,i2) = bigF1;
		omegaMesh(i1,i2) = omega;
	end
	end
return;
end
