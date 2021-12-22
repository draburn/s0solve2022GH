function matX = NEO_calcLevCurve_general( funchG, vecX0, prm=[] )
	thisFile = "NEO_calcLevCurve_general";
	%
	msg( thisFile, __LINE__, "*THIS* is promising, but needs refinement!" );
	msg( thisFile, __LINE__, "Make sure we don't stick to wrong side of ball." );
	msg( thisFile, __LINE__, "Use RK4/5 or Heun/Euler." );
	msg( thisFile, __LINE__, "Use adaptive params and step sizes." );
	msg( thisFile, __LINE__, "Stop when gradient is no longer inwards." );
	%
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	numDVals = 51;
	dHi = 10.0;
	dVals = linspace( 0.0, dHi, numDVals );
	%
	matX(:,1) = vecX0;
	for n=2:numDVals
		d = dVals(n);
		funchFlow_original = @(x,t)( -funchG(x,t) );
		funchFlow = @(x,t)( calcBallFlow( funchFlow_original, vecX0, d, x, t ) );
		%
		vecX = matX(:,n-1);
		for m=1:10
			vecFlow = funchFlow(vecX,0);
			flowNorm = norm(vecFlow);
			if ( 0.0 == flowNorm )
				break;
			end
			vecX += 0.1*vecFlow/flowNorm;
		end
		matX(:,n) = vecX;
	end
	%
return;
end
