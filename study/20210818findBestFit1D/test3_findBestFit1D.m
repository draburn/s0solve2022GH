	clear;
	commondefs;
	thisFile = "test3_findBestFit1D";
	setprngstates(0);
	msg( thisFile, __LINE__, "Starting test." );
	%
	funchF_secret = @(x)( abs(x).^3 );
	rhoArgs.xVals = linspace( 0.0, 2.0, 3 )';
	rhoArgs.fVals = funchF_secret(rhoArgs.xVals);
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
	%
	prm = [];
	echo__prm = prm
	vecZ0 = [ 0.0, 2.0 ]'
	prm = [];
	prm.matZBounds = [ -Inf, 0.0; 1.0, +Inf ];
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