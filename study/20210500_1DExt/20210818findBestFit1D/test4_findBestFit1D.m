	clear;
	commondefs;
	thisFile = "test4_findBestFit1D";
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(0);
	msg( thisFile, __LINE__, sprintf( "RETCODE__SUCCESS = %d.", RETCODE__SUCCESS ) );
	f0 = 0.2;
	f1 = 0.4;
	x0 = 0.6;
	p = 3.5;
	funchF_secret = @(x)( f0+f1*abs(x-x0).^p );
	%
	rhoArgs.xVals = linspace( -1.0, 2.0, 5 )';
	rhoArgs.fVals = funchF_secret(rhoArgs.xVals);
	rhoArgs.dVals = 0.4*ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
	%
	vecZ0 = [ 0.0, 2.0 ]';
	prm = [];
	prm.useCustomOmega = true;
	prm.funchOmega      = @(rho)( sum(log( 1.0 + rho.^2 )) );
	prm.funchOmegaDRho  = @(rho)( 2.0*rho ./ ( 1.0 + rho.^2 ) );
	prm.funchOmegaDRho2 = @(rho)(diag( 2.0*(1.0-rho.^2)./((1.0+rho.^2).^2) ));
	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm );
	echo__vecZ = vecZ
	echo__retCode = retCode;
	echo__datOut = datOut;
	assert( retCode == RETCODE__SUCCESS );
	msg( thisFile, __LINE__, "Finished test." );
	%
	%
	%
	xLo = min([ min(rhoArgs.xVals), vecZ(1) ]);
	xHi = max([ max(rhoArgs.xVals), vecZ(1) ]);
	numPts_viz = 1000;
	xVals_viz = linspace( xLo, xHi, numPts_viz );
	fVals_viz = funchF_secret(xVals_viz);
	[ errFlag, vecRho, bigF0, bigF1 ] = absPowSym_calcStuff( rhoArgs, vecZ )
	assert( ~errFlag );
	funchFModel = @(x)( bigF0 + bigF1*abs( x - vecZ(1) ).^vecZ(2) );
	fModelVals_viz_raw = funchFModel( xVals_viz );
	fLo = min(fVals_viz) - 0.5*(max(fVals_viz)-min(fVals_viz));
	fHi = max(fVals_viz) + 0.5*(max(fVals_viz)-min(fVals_viz));
	fModelVals_viz = cap( fModelVals_viz_raw, fLo, fHi );
	plot( ...
	  rhoArgs.xVals, rhoArgs.fVals, 'ko', 'markerSize', 25, 'linewidth', 5, ...
	  rhoArgs.xVals, rhoArgs.fVals, 'k+', 'markerSize', 25, 'linewidth', 5, ...
	  xVals_viz, fVals_viz, 's-', 'linewidth', 2, ...
	  xVals_viz, fModelVals_viz, 'x-', 'linewidth', 2 );
	grid on;
%
