function matX = NEO_calcLevCurve_lsodeOnSurface( funchG, vecX0, prm=[] )
	thisFile = "NEO_calcLevCurve_lsodeOnSurface";
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
	%barrierThickness = dHi/(5.0*(numDVals-1.0));
	barrierThickness = 0.01;
	%
	matX(:,1) = vecX0;
	for n=2:numDVals
		d = dVals(n);
		funchAntiG = @(x,t)( -calcBallSurfGrad( funchG, vecX0, d, x, t ) );
		%
		%matX_temp = lsode( funchAntiG, matX(:,n-1), [0.0,1.0] )';
		%matX(:,n) = matX_temp(:,2);
		%
		vecX = matX(:,n-1);
		for m=1:30
			vecG = funchAntiG(vecX,0);
			gNorm = norm(vecG);
			if ( 0.0 == gNorm )
				break;
			end
			vecX += 0.05*vecG/gNorm;
		end
		matX(:,n) = vecX;
	end
	%
return;
end
