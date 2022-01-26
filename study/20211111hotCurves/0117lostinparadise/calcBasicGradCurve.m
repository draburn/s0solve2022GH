function [ vecXVals, datOut ] = calcBasicGradCurve( vecX0, omega0, vecG0, matH, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(omega0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	%
	matS = mygetfield( prm, "matS", [] );
	if ( ~isempty(matS) )
		msg( __FILE__, __LINE__, "Calling with GUESS scaled values." );
		prmMod = prm;
		prmMod.matS = [];
		matSInv = inv(matS);
		vecGMod = matSInv*vecG0;
		matHMod = matSInv'*matH*matSInv;
		vecXVals = calcBasicGradCurve( vecX0, omega0, vecGMod, matHMod, prmMod );
		vecXVals = vecX0 + matSInv*(vecXVals-vecX0);
		return;
	end
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
	numVals = 401;
	sVals = (linspace( 1.0, 0.0, numVals )).^(1.0/lambdaAbsMin);
	omegaMin = 0.0;
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
	skipVals = nanVals | (omegaVals<omegaMin);
	if ( 0 == sum(double(skipVals)) )
		return;
	end
	clear vecDVals;
	clear omegaVals;
	%
	msg( __FILE__, __LINE__, "Discarding invalid points." );
	vecXVals = vecXVals(:,~skipVals);
	%
	numVals = size(vecXVals,2);
	if ( 1 == numVals )
		msg( __FILE__, __LINE__, "Only one point was accepted." );
		return;
	end
	%
	msg( __FILE__, __LINE__, "Extrapolating to omegaMin." );
	vecXA = vecXVals(:,end-1);
	vecXB = vecXVals(:,end);
	vecZ = vecXB-vecXA; % We'll take a step in this direction...
	%
	vecDB = vecXB - vecX0;
	omegaB = omega0 + (vecG0'*vecDB) + 0.5*(vecDB'*matH*vecDB);
	vecGB = vecG0 + matH*vecDB;
	%
	% DRaburn 2022.01.25:
	%  I'm not 100%, given the above, that the s here will always hit omegaMin;
	%  but, when it doesn't, this should still be a meaningful point... Right?
	s = calcLinishRootOfQuad( 0.5*(vecZ'*matH*vecZ), vecGB'*vecZ, omegaB-omegaMin );
	assert( isrealscalar(s) );
	if ( s < 0.0 )
		msg( __FILE__, __LINE__, "Extrapolation goes in wrong direction!" );
	elseif ( s < sqrt(eps) )
		msg( __FILE__, __LINE__, "Extrapolation would be too small!" );
	elseif ( s > 1.0/sqrt(eps) )
		msg( __FILE__, __LINE__, "Extrapolation would be too large!" );
	end
	vecXVals(:,end+1) = vecXB + s*vecZ;
return;
end
