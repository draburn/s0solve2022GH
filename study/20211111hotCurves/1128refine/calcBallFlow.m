function vecFlow = calcBallFlow( funchFlow, vecX0, bigR, vecX, t=0.0, prm=[] )
	% Flow is the negative of the gradient.
	thisFile = "calcBallFlow";
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
	vecFlow = funchFlow(vecX,t);
	flowNorm = norm(vecFlow);
	assert( 0.0 < flowNorm );
	vecDelta = vecX - vecX0;
	deltaNorm = norm(vecDelta);
	if ( 0.0 == deltaNorm )
		return;
	end
	vecDeltaHat = vecDelta / deltaNorm;
	z = deltaNorm - bigR;
	vecSurfaceAttractionFlow = -vecDeltaHat * z;
	%
	% Only modify the outward flow if it exceeds the surface attraction...
	if ( vecDeltaHat'*vecFlow > vecDeltaHat'*vecSurfaceAttractionFlow )
		vecFlow = vecFlow - vecDeltaHat*(vecDeltaHat'*vecFlow) + vecSurfaceAttractionFlow;
	end
	%
return;
end
