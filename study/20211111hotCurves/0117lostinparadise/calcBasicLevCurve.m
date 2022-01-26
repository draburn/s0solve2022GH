function [ vecXVals, datOut ] = calcBasicLevCurve( vecX0, omega0, vecG0, matH, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(omega0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	%
	matS = mygetfield( prm, "matS", eye(sizeX,sizeX) );
	matD = matS'*matS;
	%
	datOut = [];
	eps075 = eps^0.75;
	omegaMin = 0.0;
	%
	numVals = 401;
	for n=1:numVals
		s = (n-1.0)/(numVals-1.0);
		matM = (s*matH) + ((1.0-s)*matD);
		[ matR, cholFlag ] = chol( matM );
		if (0~=cholFlag)
			msg( __FILE__, __LINE__, sprintf( "chol() failed for s = %f.", s ) );
			break;
		end
		vecDelta = matR \ ( matR' \ ( -s*vecG0) );
		%
		omega = omega0 + vecG0'*vecDelta + 0.5*vecDelta'*matH*vecDelta;
		if ( omega < omegaMin )
			msg( __FILE__, __LINE__, sprintf( "omega = %f.", omega ) );
			break;
		end
		%
		vecXVals(:,n) = vecX0 + vecDelta;
	end
	%
	if ( 1+numVals == n )
		% All points were accepted.
		return;
	elseif ( 2 == n )
		msg( __FILE__, __LINE__, "Only one point was accepted." );
		return;
	end
	vecXVals = vecXVals(:,1:n-1);
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
