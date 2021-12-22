function vecGBallSurf = calcBallSurfGrad( funchG, vecX0, bigR, vecX, t=0.0, prm=[] )
	% t is an argument of funchG.
	thisFile = "calcBallSurfGrad";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(bigR) );
	assert( 0.0 < bigR );
	assert( isrealarray(vecX,[sizeX,1]) );
	%
	normScale = mygetfield( prm, "normScale", 1.0 );
	assert( isrealscalar(normScale) );
	assert( 0.0 < normScale );
	%
	vecDelta = vecX - vecX0;
	deltaNorm = norm(vecDelta);
	%if ( deltaNorm < bigR )
	%	vecGBallSurf = funchG(vecX,t);
	%	return;
	%end
	if ( 0.0 == deltaNorm )
		vecGTemp = funchG(vecX0,t);
		gTempNorm = norm(vecGTemp);
		assert( 0.0 ~= gTempNorm );
		vecDeltaHat = vecGTemp/gTempNorm;
	else
		vecDeltaHat = vecDelta/deltaNorm;
	end
	vecXBallSurf = vecX0 + bigR*vecDeltaHat;
	vecG = funchG(vecXBallSurf,t);
	vecGBallSurf = vecG - vecDeltaHat * ( vecDeltaHat'*vecG ) + normScale * vecDeltaHat * ( deltaNorm - bigR );
	%
return;
end
