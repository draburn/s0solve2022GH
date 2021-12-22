function matX = NEO_calcLevCurve_lsodePlusCost( funchG, vecX0, prm=[] )
	thisFile = "NEO_calcLevCurve_lsodePlusCost";
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	%
	numSteps = mygetfield( prm, "numSteps", 100 );
	assert( isposintscalar(numSteps) );
	finalStepEps = mygetfield( prm, "finalStepEps", 0.5 );
	assert( isrealscalar(finalStepEps) );
	assert( 0.0 <= finalStepEps );
	%
	matI = eye(sizeX,sizeX);
	n = 1;
	vecX = vecX0;
	while (1)
		matX(:,n) = vecX;
		n++;
		if ( n > numSteps )
			break;
		end
		s = n/(numSteps+finalStepEps);
		funchAntiG = @(x,t)( -(s*funchG(x,t)+(1.0-s)*(x-vecX0)) );
		matX_temp = lsode( funchAntiG, matX(:,n-1), [0.0,5.0] )';
		vecX = matX_temp(:,2);
	end
	%funchAntiG = @(x)( -funchG(x) );
	%matX = lsode( funchAntiG, vecX0, 1000.0*linspace(0.0,1.0,1001) )';
	%
return;
end
