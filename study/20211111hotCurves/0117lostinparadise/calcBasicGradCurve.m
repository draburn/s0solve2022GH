function [ vecXVals, datOut ] = calcBasicGradCurve( vecX0, omega0, vecG0, matH, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(omega0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	%
	datOut = [];
	%
	[ matPsi, matLambda ] = eig( matH );
	vecLambda = diag(matLambda);
	%
	vecY0 = matPsi'*vecX0;
	vecGamma = (matPsi'*vecG0) - (matLambda*vecY0);
	%
	numVals = 1001;
	sVals = linspace( 0.0, 1.0, numVals );
	%
	sizeY = sizeX;
	vecGOL = vecGamma./vecLambda;
	vecYVals = zeros(sizeY,numVals);
	for n=1:sizeY
	if ( 0.0 == vecLambda(n) )
		vecYVals(n,:) = vecY0(n) + vecGamma(n)*log(sVals);
	else
		vecYVals(n,:) = ( vecY0(n) + vecGOL(n) ) * (sVals.^vecLambda(n)) - vecGOL(n);
	end
	end
	%
	ySumSqVals = sumsq( vecYVals, 1 );
	nanVals = ( (isnan(ySumSqVals)) | isinf(ySumSqVals) );
	ySumSqMax = 100.0;
	skipVals = nanVals | (ySumSqVals>ySumSqMax);
	vecYVals_trim = vecYVals(:,~skipVals);
	vecXVals = matPsi*vecYVals_trim;
return;
end
