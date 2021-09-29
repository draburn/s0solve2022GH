	clear;
	commondefs;
	thisFile = "test10_studyAccel";
	numFigs = 0;
	%
	%
	%
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(73798480); % Fails due to bad initial guess.
	f0 = randn();
	f1 = randn();
	x0 = randn()
	p = abs(randn());
	funchF_secret = @(x)( f0 + f1*(abs(x-x0).^p) );
	%
	%
	%
	numPts = 5 + round(abs(4.0*randn()))
	rhoArgs.xVals = x0 + randn(numPts,1);
	rhoArgs.xVals(1:2) = x0 - abs(randn(2,1));
	rhoArgs.xVals(3:4) = x0 + abs(randn(2,1));
	rhoArgs.xVals = sort(rhoArgs.xVals);
	numPts = size(rhoArgs.xVals,1);
	rhoArgs.fVals = funchF_secret(rhoArgs.xVals);
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
	%
	prm = [];
	%vecZ0 = [ x0 + 0.5*randn(); 1.0 ];
	vecZ0 = [ x0 + 0.5*randn(); 2.0 ];
	%prm.useCustomOmega = true;
	%prm.funchOmega      = @(rho)( sum(log( 1.0 + rho.^2 )) );
	%prm.funchOmegaDRho  = @(rho)( 2.0*rho ./ ( 1.0 + rho.^2 ) );
	%prm.funchOmegaDRho2 = @(rho)(diag( 2.0*(1.0-rho.^2)./((1.0+rho.^2).^2) ));
	%prm.matZBounds = [ -Inf, +Inf; 0.0, +Inf; 0.0, +Inf ];
	%prm.iterLimit = 1;
	%
	%
	%
	xLo = min(rhoArgs.xVals);
	xHi = max(rhoArgs.xVals);
	numPts_viz = 1000;
	xVals_viz = linspace( xLo, xHi, numPts_viz );
	fVals_viz = funchF_secret(xVals_viz);
	%
	vecZ = vecZ0;
	[ errFlag, vecRho, bigF0, bigF1 ] = absPowSym_calcStuff( rhoArgs, vecZ );
	assert( ~errFlag );
	funchFModel = @(x)( bigF0 + bigF1 *( abs(x-vecZ(1)).^vecZ(2) ) );
	fModelVals_viz_raw = funchFModel( xVals_viz );
	fLo = min(fVals_viz) - 0.5*(max(fVals_viz)-min(fVals_viz));
	fHi = max(fVals_viz) + 0.5*(max(fVals_viz)-min(fVals_viz));
	fModelVals_viz = cap( fModelVals_viz_raw, fLo, fHi );
	bigF0_capped = cap( bigF0, fLo, fHi );
	%
	numFigs++; figure(numFigs);
	plot( ...
	  rhoArgs.xVals, rhoArgs.fVals, 'ko', 'markerSize', 25, 'linewidth', 2, ...
	  rhoArgs.xVals, rhoArgs.fVals, 'k+', 'markerSize', 25, 'linewidth', 2, ...
	  x0, f0, 'k+', 'markersize', 25, 'linewidth', 2, ...
	  xVals_viz, fVals_viz, '-', 'linewidth', 10, ...
	  vecZ(1), bigF0_capped, 'rx', 'markersize', 25, 'linewidth', 2, ...
	  xVals_viz, fModelVals_viz, 'r-', 'linewidth', 4 );
	grid on;
	clear vecZ;
	%
	%
	%
	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm );
	echo__vecZ = vecZ
	echo__retCode = retCode;
	echo__datOut = datOut;
	msg( thisFile, __LINE__, sprintf( "findBestFit1D() returned %s.", retcode2str(retCode) ) );
	msg( thisFile, __LINE__, "Finished test." );
	%
	%
	%
	[ errFlag, vecRho, bigF0, bigF1 ] = absPowSym_calcStuff( rhoArgs, vecZ );
	assert( ~errFlag );
	funchFModel = @(x)( bigF0 + bigF1 *( abs(x-vecZ(1)).^vecZ(2) ) );
	fModelVals_viz_raw = funchFModel( xVals_viz );
	fLo = min(fVals_viz) - 0.5*(max(fVals_viz)-min(fVals_viz));
	fHi = max(fVals_viz) + 0.5*(max(fVals_viz)-min(fVals_viz));
	fModelVals_viz = cap( fModelVals_viz_raw, fLo, fHi );
	bigF0_capped = cap( bigF0, fLo, fHi );
	numFigs++; figure(numFigs);
	plot( ...
	  rhoArgs.xVals, rhoArgs.fVals, 'ko', 'markerSize', 25, 'linewidth', 2, ...
	  rhoArgs.xVals, rhoArgs.fVals, 'k+', 'markerSize', 25, 'linewidth', 2, ...
	  x0, f0, 'k+', 'markersize', 25, 'linewidth', 2, ...
	  xVals_viz, fVals_viz, '-', 'linewidth', 10, ...
	  vecZ(1), bigF0_capped, 'rx', 'markersize', 25, 'linewidth', 2, ...
	  xVals_viz, fModelVals_viz, 'r-', 'linewidth', 4 );
	grid on;
