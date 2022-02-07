% Function...
%  [ vecX, retCode, datOut ] = findOmegaMinWithinSurf( vecX0, funchSurf, funchOmegaG, prm=[] )
% Unless prm.useProvidedGradients is set to FALSE, funchOmegaG should support the following interface.
%  [ omega, vecG ] = funchOmegaG( vecX )

function [ vecX, retCode, datOut ] = findOmegaMinWithinSurf( vecX0, funchSurf, funchOmegaG, prm=[] )
	debugMode = mygetfield( prm, "debugMode", false );
	%
	% Validate main input.
	sizeX = size(vecX0,1);
	if ( debugMode )
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
	end
	[ omega0, vecNablaOmega0 ] = funchOmegaG( vecX0 );
	if ( debugMode )
		[ vecS0, vecU0, vecV0, matNablaST0 ] = funchSurf( vecX0 );
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
	% Set up the objective function that combines the input funchSurf and funchOmegaG.
	tauX = mygetfield( prm, "tauX", 1e-2 );
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
	prm_fowis = [];
	prm_fowis.h0 = h0;
	prm_fowis.tau = tauX;
	prm_fowis.epsX = epsX;
	funchBigF = @(dummyX)( funcOmegaWithinSurf( dummyX, funchOmegaG, funchSurf, prm_fowis ) );
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
