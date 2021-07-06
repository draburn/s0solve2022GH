function [ omegaMesh ] = extFinder_getOmegaMesh( ...
  bigSMesh, bigPMesh, xVals, gVals, nOfPtWiseMin, prm=[] );
	thisFile = "extFinder_getOmegaMesh";
	%
	size1 = size(bigSMesh,1);
	size2 = size(bigSMesh,2);
	assert( isrealarray(bigSMesh,[size1,size2]) );
	assert( isrealarray(bigPMesh,[size1,size2]) );
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
	n = nOfPtWiseMin;
	xOfPtWiseMin = xVals(n);
	bigDelta = xVals(n+1) - xVals(n-1);
	bigG0 = gVals(n);
	bigG1 = gVals(n+1) + gVals(n-1) - 2.0*gVals(n);
	vecX = xVals(n-1:n+1)';
	vecG = gVals(n-1:n+1)';
	vecY = ( vecG - bigG0 ) / bigG1;
	%
	maskVals = logical(ones(size(xVals)));
	maskVals(nOfPtWiseMin-1:nOfPtWiseMin+1) = 0;
	%wVals = mygetfield( prm, "wVals", (1.0/sum(1.0./sqrt(gVals(maskVals))))./sqrt(gVals) );
	wVals = mygetfield( prm, "wVals", [] );
	if ( isempty(wVals) )
		wVals = 1.0./gVals;
		wVals /= sum(wVals);
	end
	assert( isrealarray(wVals,[1,numPts]) );
	noWValuesAreNegative = ( 0 == sum(wVals<0.0) );
	assert( noWValuesAreNegative );
	%
	wSus = mygetfield( prm, "wSus", 1.0 );
	assert( isrealscalar(wSus) );
	assert( wSus >= 0.0 );
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
	%
	%
	%
	myExpMesh = 1.0./(bigPMesh-1.0);
	for i1=1:size1
	for i2=1:size2
		vecD = ( vecX - bigSMesh(i1,i2) ) / bigDelta;
		vecCoeff = [ abs(vecD).^bigPMesh(i1,i2), vecD, ones(3,1) ] \ vecY;
		funchGModel = @(x)( bigG0 + bigG1*( ...
		   vecCoeff(1)*abs((x-bigSMesh(i1,i2))/bigDelta).^bigPMesh(i1,i2) ...
		 + vecCoeff(2)*(x-bigSMesh(i1,i2))/bigDelta ...
		 + vecCoeff(3) ) );
		foo = vecCoeff(2)/(vecCoeff(1)*bigPMesh(i1,i2));
		xExt = bigSMesh(i1,i2) - bigDelta*sign(foo)*abs(foo)^myExpMesh(i1,i2);
		gExt = funchGModel( xExt );
		%
		rhoWVals = wVals.*(funchGModel(xVals)-gVals).^2;
		omegaFit = 0.5 * sum( rhoWVals(maskVals) );
		omegaSus = ( gExt < 0.0 ) * 0.5 * wSus * ( gExt / bigG1 )^2;
		omegaRegu = 0.5 * ( ...
		 + epsA*vecCoeff(1)^2 ...
		 + epsB*vecCoeff(2)^2 ...
		 + epsC*vecCoeff(3)^2 ...
		 + epsS*(bigSMesh(i1,i2)-xOfPtWiseMin)^2 ...
		 + epsP*(bigPMesh(i1,i2)-2.0)^2 );
		%omegaMesh(i1,i2) = omegaFit * ( 1.0 + omegaSus ) + omegaRegu;
		omegaMesh(i1,i2) = omegaFit;% + omegaRegu/sqrt(eps);
	end
	end
	%
	%
	%
return;
end
