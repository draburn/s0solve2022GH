function [ vecX, retCode, datOut ] = findOmegaMinWithinSurf( vecX0, funchSurf, funchOmega, prm=[] )
	debugMode = mygetfield( prm, "debugMode", true );
	%
	% Validate main input.
	sizeX = size(vecX0,1);
	if ( debugMode )
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
		[ vecS0, vecU0, vecV0, matNablaST0 ] = funchSurf( vecX0 );
		[ omega0, vecNablaOmega0 ] = funchOmega( vecX0 );
		assert( isrealarray(vecS0,[sizeX,1]) );
		assert( isrealarray(vecU0,[sizeX,1]) );
		assert( isrealarray(vecV0,[sizeX,1]) );
		assert( isrealarray(matNablaST0,[sizeX,sizeX]) );
		assert( abs(norm(vecU0)-1.0) < sizeX*(eps^0.75) );
		assert( abs(norm(vecV0)-1.0) < sizeX*(eps^0.75) );
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecNablaOmega0,[sizeX,1]) );
	end
	datOut = [];
	%
	%
	% Set up the objective function that combines the input funchSurf and funchOmega.
	tauX = mygetfield( prm, "tauX", 1e-3 );
	h0 = mygetfield( prm, "h0", norm(vecNablaOmega0)/tauX );
	epsX = mygetfield( prm, "epsX", 1e-5 );
	useProvidedGradients = mygetfield( prm, "useProvidedGradients", true );
	if ( debugMode )
		assert( isrealscalar(tauX) );
		assert( isrealscalar(h0) );
		assert( isrealscalar(epsX) );
		assert( 0.0 < epsX );
		assert( isscalar(useProvidedGradients) );
		assert( isbool(useProvidedGradients) );
	end
	prm.h0 = h0;
	prm.tau = tauX;
	prm.epsX = epsX;
	funchBigF = @(dummyX)( funcOmegaWithinSurf( dummyX, funchOmega, funchSurf, prm ) );
	%
	%
	% Do work.
	if (useProvidedGradients)
		fminunc_opts = optimset( 'GradObj', 'on' );
		vecX = fminunc( funchBigF, vecX0, fminunc_opts );
	else
		vecX = fminunc( funchBigF, vecX0 );
	end
	if ( debugMode )
		assert( isrealarray(vecX,[sizeX,1]) );
	end
	retCode = 0;
return;
end
