function [ bigF0Mesh, bigF1Mesh, omegaMesh ] = extFit__calcMesh( ...
  sMesh, pMesh, xVals, fVals, nExactFit, wVals=[], prm=[] )
	thisFile = "extFit__calcMesh";
	doChecks = mygetfield( prm, "doChecks", true );
	%
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	size1 = size(sMesh,1);
	size2 = size(sMesh,2);
	numPts = size(xVals,2);
	if ( doChecks )
		assert( isrealarray(sMesh,[size1,size2]) );
		assert( isrealarray(pMesh,[size1,size2]) );
		allPMeshValsArePositive = ( 0 == sum(sum(pMesh<=0.0)) );
		assert( allPMeshValsArePositive );
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isrealscalar(nExactFit) );
		assert( fleq(nExactFit,round(nExactFit)) );
		assert( 1 <= nExactFit );
		assert( nExactFit <= numPts );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	% Least-squares fit to f = F0 + F1 * | x - s |^p,
	% subject to an exact match at nExactFit.
	%
	% For comparison, use:
	%prm.doChecks = false;
	%for i1=1:size1
	%for i2=1:size2
	%	[ rhoVals, bigF0, bigF1, omega ] = extFit__calcAtPt( ...
	%	  sMesh(i1,i2), pMesh(i1,i2), xVals, fVals, nExactFit, wVals, prm );
	%	bigF0Mesh(i1,i2) = bigF0;
	%	bigF1Mesh(i1,i2) = bigF1;
	%	omegaMesh(i1,i2) = omega;
	%end
	%end
	%return;
	%
	% Optimize for meshes, loop over vals...
	dfVals = fVals - fVals(nExactFit);
	numerMesh = zeros(size1,size2);
	denomMesh = zeros(size1,size2);
	for n=1:numPts
		dyMesh = abs( xVals(n) - sMesh ).^pMesh - abs( xVals(nExactFit) - sMesh).^pMesh;
		numerMesh += wVals(n) * dyMesh * dfVals(n);
		denomMesh += wVals(n) * dyMesh.^2;
	end
	%
	bigF1Mesh = numerMesh ./ denomMesh;
	bigF0Mesh = fVals(nExactFit) - bigF1Mesh .* abs( xVals(nExactFit) - sMesh).^pMesh;
	%
	omegaMesh = zeros( size1, size2 );
	for n=1:numPts
		omegaMesh += wVals(n) * ( ...
		  bigF0Mesh + bigF1Mesh .* abs( xVals(n) - sMesh ).^pMesh - fVals(n) ).^2;
	end
	omegaMesh *= 0.5;
return;
end
