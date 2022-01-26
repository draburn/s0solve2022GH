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
	lambdaAbsMin = min(abs(vecLambda))+sqrt(eps)*max(abs(vecLambda));
	numVals = 101;
	sVals = (linspace( 0.0, 1.0, numVals )).^(1.0/lambdaAbsMin);
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
	%echo__vecYVals = vecYVals
	vecXVals = matPsi*vecYVals;
	vecDVals = vecXVals - vecX0;
	omegaVals = omega0 + (vecG0'*vecDVals) + 0.5 * sum( vecDVals.*(matH*vecDVals), 1 );
	%
	ySumSqVals = sumsq( vecYVals, 1 );
	nanVals = ( (isnan(ySumSqVals)) | isinf(ySumSqVals) );
	skipVals = nanVals | (omegaVals<0);
	if ( sum(double(skipVals))>0 )
		msg( __FILE__, __LINE__, "Discarding invalid points." );
	end
	vecXVals = vecXVals(:,~skipVals);
return;
end
