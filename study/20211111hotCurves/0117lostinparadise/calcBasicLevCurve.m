function [ vecXVals, datOut ] = calcBasicLevCurve( vecX0, omega0, vecG0, matH, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(omega0) );
	assert( isrealarray(vecG0,[sizeX,1]) );
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	%
	datOut = [];
	eps075 = eps^0.75;
	%
	matI = eye(sizeX,sizeX);
	numVals = 1001;
	for n=1:numVals
		s = (n-1.0)/(numVals-1.0);
		matM = (s*matH) + ((1.0-s)*matI);
		[ matR, cholFlag ] = chol( matM );
		if (0~=cholFlag)
			msg( __FILE__, __LINE__, sprintf( "chol() failed for s = %f.", s ) );
			return;
		end
		vecDelta = matR \ ( matR' \ ( -s*vecG0) );
		%
		omega = omega0 + vecG0'*vecDelta + 0.5*vecDelta'*matH*vecDelta;
		if ( omega < -abs(eps075*omega0) )
			msg( __FILE__, __LINE__, sprintf( "omega = %f.", omega ) );
			return;
		end
		%
		vecXVals(:,n) = vecX0 + vecDelta;
	end
return;
end
