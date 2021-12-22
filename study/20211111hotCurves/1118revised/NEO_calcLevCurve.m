function matX = NEO_calcLevCurve( funcPrm, vecX0, prm=[] )
	thisFile = "NEO_calcLevCurve";
	%
	% Using "barrier" method...
	numDVals = 1001;
	dVals = 10.0*linspace( 0.0, 1.0, numDVals );
	if (0)
		assert( 0 );
		funchAntiG = @(x)( -testFunc_gForLSODE(x,funcPrm) );
		matX = lsode( funchAntiG, vecX0, 1000.0*linspace(0.0,1.0,1001) );
	elseif (1)
		matX(1,:) = vecX0;
		for n=2:numDVals
			%funchAntiG = @(x)( -testFunc_gForLSODE(x,funcPrm) );
			funchAntiG = @(x)( -testFunc_gLevForLSODE(x,funcPrm,vecX0,dVals(n)) );
			matX_temp = lsode( funchAntiG, matX(n-1,:), 0.1*linspace(0.0,1.0,10) );
			matX(n,:) = matX_temp(end,:);
		end
	end
return;
end
