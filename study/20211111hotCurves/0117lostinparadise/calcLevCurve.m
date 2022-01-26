function [ vecXVals, datOut ] = calcLevCurve( vecX0, funchOmega, prm=[] )
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	numSVals = 101;
	vecS = linspace( 0.0, 1.0, numSVals );
	%
	matS = mygetfield( prm, "matS", [] );
	if (isempty(matS))
		matD = [];
	else
		matD = matS'*matS;
	end
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
	function [ omegaVals, vecNablaOmegaVals ] = funchOmegaLev_withS( vecXVals, funchOmegaBase, s, vecXCent, matD )
		vecDVals = vecXVals - vecXCent;
		if ( 1 == nargout )
			omegaBaseVals = funchOmegaBase( vecXVals );
			omegaVals = s*omegaBaseVals + 0.5*(1.0-s)*sum( vecDVals'.*matD*vecDVals, 1 );
		else
			[ omegaBaseVals, vecNablaOmegaBaseVals ] = funchOmegaBase( vecXVals );
			omegaVals = s*omegaBaseVals + 0.5*(1.0-s)*sum( vecDVals'.*matD*vecDVals, 1 );
			vecNablaOmegaVals = s*vecNablaOmegaBaseVals + (1.0-s)*matD*vecDVals;
		end
	return;
	end
	%
	vecX = vecX0;
	n = 1;
	s = vecS(n);
	vecXVals(:,n) = vecX0;
	%
	for n=2:numSVals
		s = vecS(n);
		if ( isempty(matD) )
			fcn = @(dummyX)( funchOmegaLev_sansS( dummyX, funchOmega, s, vecX0 ) );
		else
			fcn = @(dummyX)( funchOmegaLev_withS( dummyX, funchOmega, s, vecX0, matD ) );
		end
		opts = optimset( 'GradObj', 'on' );
		vecX = fminunc( fcn, vecX, opts );
		vecXVals(:,n) = vecX;
	end
	%
return;
end
