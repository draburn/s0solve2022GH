	clear;
	commondefs;
	thisFile = "test6_asym";
	msg( thisFile, __LINE__, "Starting test." );
	%setprngstates();
	setprngstates(17315808);
	msg( thisFile, __LINE__, sprintf( "RETCODE__SUCCESS = %d.", RETCODE__SUCCESS ) );
	f0 = 0.0%randn();
	fL = 1.0%randn();
	fR = 1.0%randn();
	x0 = 0.0%randn();
	pL = 2.0%2.0+abs(randn());
	pR = 0.5%2.0+abs(randn());
	funchF_secret = @(x)( f0 ...
	  + (x<x0).*(fL*abs(x-x0).^pL) ...
	  + (x>x0).*(fR*abs(x-x0).^pR) );
	%
	rhoArgs.xVals = linspace( -1.0, 1.0, 201 )';
	rhoArgs.fVals = funchF_secret(rhoArgs.xVals);
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowAsym( ra, z );
	%
	[ errFlag, vecRho ] = funchRho( rhoArgs, [0.0; 2.0; 0.5] );
	%
	plot( ...
	  rhoArgs.xVals, rhoArgs.fVals, 'o-', ...
	  rhoArgs.xVals, rhoArgs.fVals+vecRho, 'x-' );
	grid on;	
	%
	return;
	%
	%%%vecZ0 = [ 0.0, 2.0 ]';
	%vecZ0 = [ 1.4, 2.8 ]';
	%vecZ0 = [ 1.4, 2.0 ]';
	%vecZ0 = [ 0.6, 2.0 ]';
	%vecZ0 = [ 0.3, 2.0 ]';
	vecZ0 = [ 0.249, 2.0 ]';
	%vecZ0 = [ 0.1, 2.0 ]';
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
