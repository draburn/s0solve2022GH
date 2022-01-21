function [ vecXVals, datOut ] = calcLevCurve( vecX0, funchOmega, prm=[] )
	%
	n = 1;
	vecX = vecX0;
	vecXVals(:,n) = vecX0;
	numSVals = 51;
	vecS = linspace( 0.0, 1.0, numSVals );
	s = vecS(n);
	%printf( "%5d, %0.3f, %9.3f, %9.3f, %10.3e, %10.3e, %10.3e\n", ...
	%  n, s, vecX(1), vecX(2), funchOmega(vecX), 0.5*(vecX-vecX0)'*(vecX-vecX0), ...
	%  s*funchOmega(vecX) + 0.5*(1.0-s)*(vecX-vecX0)'*(vecX-vecX0) );
	%
	function [ omegaVals, vecNablaOmegaVals ] = funchOmegaLev_sansS( vecXVals, funchOmegaBase, s, vecXCent )
		if ( 1 == nargout )
			omegaBaseVals = funchOmegaBase( vecXVals );
			omegaVals = s*omegaBaseVals + 0.5*(1.0-s)*sumsq(vecXVals-vecXCent,1);
		else
			[ omegaBaseVals, vecNablaOmegaBaseVals ] = funchOmegaBase( vecXVals );
			omegaVals = s*omegaBaseVals + 0.5*(1.0-s)*sumsq(vecXVals-vecXCent,1);
			vecNablaOmegaVals = s*vecNablaOmegaBaseVals + (1.0-s)*(vecXVals-vecX0);
		end
	return;
	end
	%
	for n=2:numSVals
		s = vecS(n);
		fcn = @(dummyX)( funchOmegaLev_sansS( dummyX, funchOmega, s, vecX0 ) );
		opts = optimset( 'GradObj', 'on' );
		vecX = fminunc( fcn, vecX, opts );
		vecXVals(:,n) = vecX;
	end
	%
return;
end
