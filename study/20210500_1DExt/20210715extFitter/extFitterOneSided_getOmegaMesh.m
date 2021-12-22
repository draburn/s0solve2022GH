function [ omegaMesh, bigG0Mesh, bigG1Mesh, datOut ] = extFitterOneSided_getOmegaMesh( ...
  sMesh, pMesh, xVals, fVals, nOfClosestPt, prm=[] );
	thisFile = "extFitterOneSided_getOmegaMesh";
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
	bigX1 = (xVals(end)-xVals(1))/2.0;
	datOut.bigX1 = bigX1;
	%
	if (0)
		bigX0 = xVals(nOfClosestPt);
		bigF1 = (max(fVals)-min(fVals))/2.0;
		bigF0 = fVals(nOfClosestPt);
	else
	n = median([ nOfClosestPt, 2, numPts-1 ]);
	vecX = xVals(n-1:n+1)';
	vecF = fVals(n-1:n+1)';
	matX = [ ones(3,1), vecX, vecX.^2 ];
	quadFitMatrixIsNonSingular = ( rcond(matX) > eps );
	assert( quadFitMatrixIsNonSingular );
	vecC = matX \ vecF;
	quadFitIsNotLinear = ( abs(vecC(3)) > sqrt(eps)*abs(vecC(2))) ;
	assert( quadFitIsNotLinear );
	bigX0 = -vecC(2)/ (2.0*vecC(3));
	bigF1 = vecC(3)*(bigX1^2);
	bigF0 = vecC(1) - vecC(3)*(bigX0^2);
	clear vecC;
	clear matX;
	clear vecF;
	clear vecX;
	end
	datOut.bigX0 = bigX0;
	datOut.bigF1 = bigF1;
	datOut.bigF0 = bigF0;
	%
	gVals = ( fVals - bigF0 ) / bigF1;
	yVals = ( xVals - bigX0 ) / bigX1;
	%
	n = median([ nOfClosestPt, 2, numPts-1 ]);
	% Model: g = G0 + G1 * | y - Y_1 * s |^p.
	% For exact match at "nOfClosestPt"...
	yc = yVals(nOfClosestPt);
	gc = gVals(nOfClosestPt);
	% Model: g = gc + G1 * ( |y-Y1*s|^p - |yc-Y1*s|^p ).
	%  = gc + g1 * h( y; s, p; Y1, yc ).
	%
	parfor i1=1:size1
	parfor i2=1:size2
		hVals = ...
		   abs( yVals - sMesh(i1,i2) ).^pMesh(i1,i2) ...
		 - abs( yc    - sMesh(i1,i2) ).^pMesh(i1,i2);
		gmgcVals = gVals - gc;
		%bigG1Mesh(i1,i2) = sum( gmgcVals .* hVals ) ./ sum( hVals .* hVals );
		%omegaMesh(i1,i2) = sum(( bigG1Mesh(i1,i2)*hVals - gmgcVals ).^2);
		wVals = 1.0./(abs(hVals)+sqrt(eps)*max(abs(hVals)));
		wVals /= sum(wVals);
		bigG1Mesh(i1,i2) = sum( wVals .* gmgcVals .* hVals ) ./ sum( wVals .* hVals .* hVals );
		omegaMesh(i1,i2) = sum( wVals .* ( bigG1Mesh(i1,i2)*hVals - gmgcVals ).^2);
		
	end
	end
	omegaMesh /= 2.0;
	bigG0Mesh = gc - bigG1Mesh.*( abs(yc-sMesh).^pMesh );
	
	
	
	return;
	
	mesh_bigY1TimeS = bigY1*meshS;
	ary3_h = zeros( size1, size2, numPts );
	parfor n=1:numPts
		ary3_h(:,:,n) = yVals(n) - 	mesh_bigY1TimeS;
	end
	ary3_h = abs(ary3_h).^meshP;
	mesh_term2 = abs(yc - mesh_bigY1TimeS).^meshP;
	parfor n=1:numPts
		ary3_h(:,:,n) -= mesh_term2;
	end
	
	return;
	
	n = nOfPtWiseMin;
	allFValsArePositive = (0==sum( 0.0 > gVals ) );
	assert( allGValsArePositive );
	isPtWiseMin = ( (gVals(n+1)>=gVals(n)) && (gVals(n-1)>=gVals(n)) );
	assert( isPtWiseMin );
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
	matX = [ vecX.^2, vecX, ones(3,1) ];
	vecCoeff = matX\vecG;
	assert( vecCoeff(1) > 0.0 );
	xOfQuadMin = -vecCoeff(2)/(2.0*vecCoeff(1));
	%
	maskVals = logical(ones(size(xVals)));
	maskVals(nOfPtWiseMin-1:nOfPtWiseMin+1) = 0;
	%wVals = mygetfield( prm, "wVals", (1.0/sum(1.0./sqrt(gVals(maskVals))))./sqrt(gVals) );
	wVals = mygetfield( prm, "wVals", [] );
	if ( isempty(wVals) )
		wVals = 1.0./(gVals.^0.0);
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
	epsA = mygetfield( prm, "epsA", eps^1.5 );
	assert( isrealscalar( epsA ) );
	assert( 0.0 <= epsA );
	epsB = mygetfield( prm, "epsB", eps^1.5 );
	assert( isrealscalar( epsB ) );
	assert( 0.0 < epsB );
	epsC = mygetfield( prm, "epsC", eps^1.5 );
	assert( isrealscalar( epsC ) );
	assert( 0.0 <= epsC );
	epsS = mygetfield( prm, "epsS", eps^1.5 );
	assert( isrealscalar( epsS ) );
	assert( 0.0 <= epsS );
	epsP = mygetfield( prm, "epsP", eps^1.5 );
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
		 + epsS*((bigSMesh(i1,i2)-xOfQuadMin)/bigDelta)^2 ...
		 + epsP*(bigPMesh(i1,i2)-2.0)^2 );
		omegaMesh(i1,i2) = omegaFit + omegaRegu;
		%omegaMesh(i1,i2) = omegaFit * ( 1.0 + omegaSus ) + omegaRegu;
		%omegaMesh(i1,i2) = omegaFit;% + omegaRegu/sqrt(eps);
		bigAMesh(i1,i2) = vecCoeff(1);
		bigBMesh(i1,i2) = vecCoeff(2);
		bigCMesh(i1,i2) = vecCoeff(3);
	end
	end
	omegaMesh = sqrt( omegaMesh );
	%
	%
	%
return;
end
