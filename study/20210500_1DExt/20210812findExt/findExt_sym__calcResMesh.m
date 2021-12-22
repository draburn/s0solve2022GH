function [ omegaMesh, fExtMesh, f1Mesh ] = findExt_sym__calcResMesh( xExtMesh, pMesh, xVals, fVals, wVals )
	size1 = size(xExtMesh,1);
	size2 = size(xExtMesh,2);
	numPts = size(xVals,2);
	assert( isrealarray(xExtMesh,[size1,size2]) );
	assert( isrealarray(pMesh,[size1,size2]) );
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
	assert( isrealarray(wVals,[1,numPts]) );
	noWValIsNegative = (0==sum(wVals<0.0));
	assert( noWValIsNegative );
	assert( sum(wVals) > 0.0 );
	%
	sigma1 = sum(wVals);
	sigmaF = sum(wVals.*fVals);
	sigmaY = zeros(size1,size2);
	sigmaYSq = zeros(size1,size2);
	sigmaYF = zeros(size1,size2);
	for n=1:numPts
		yNMesh = abs( xVals(n) - xExtMesh ).^pMesh;
		sigmaY += wVals(n)*yNMesh;
		sigmaYSq += wVals(n)*(yNMesh.^2);
		sigmaYF += wVals(n)*yNMesh*fVals(n);
	end
	denomMesh = sigma1*sigmaYSq - (sigmaY.^2);
	fExtMesh = ( sigmaYSq*sigmaF - sigmaY.*sigmaYF ) ./ denomMesh;
	f1Mesh = ( sigma1*sigmaYF - sigmaY*sigmaF ) ./ denomMesh;
	%
	omegaMesh = zeros(size1,size2);
	for n=1:numPts
		yNMesh = abs( xVals(n) - xExtMesh ).^pMesh;
		rhoNMesh = fExtMesh + f1Mesh.*yNMesh - fVals(n);
		omegaMesh += wVals(n)*(rhoNMesh.^2);
	end
	omegaMesh *= 0.5;
return;
end
