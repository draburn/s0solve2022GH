	clear;
	commondefs;
	thisFile = "test7a_asym";
	numFigs = 0;
	%
	% FOR THE TEST CASE THAT TRIGGERED THE "_TWEAK" STUFF, SEE test7_asym.
	use6723760Tweaks = true;
	%
	%
	%
	msg( thisFile, __LINE__, "Starting test." );
	%setprngstates();
	setprngstates(6723760);
	%setprngstates(17315808);
	msg( thisFile, __LINE__, sprintf( "RETCODE__SUCCESS = %d.", RETCODE__SUCCESS ) );
	f0 = randn();
	fL = abs(randn());
	fR = abs(randn());
	x0 = randn()
	pL = abs(randn())
	pR = abs(randn())
	funchF_secret = @(x)( f0 ...
	  + (x<x0).*(fL*abs(x-x0).^pL) ...
	  + (x>x0).*(fR*abs(x-x0).^pR) );
	vecZSecret = [ x0; pL; pR ];
	assert( fL*fR > 0.0 );
	%
	%
	%
	numPts = 6 + round(4*randn());
	rhoArgs.xVals = x0 + randn(numPts,1);
	rhoArgs.xVals(1:3) = x0 - abs(randn(3,1));
	rhoArgs.xVals(4:6) = x0 + abs(randn(3,1));
	%numPts = 10;
	%rhoArgs.xVals = x0+linspace( -2.0, 2.0, numPts )';
	rhoArgs.fVals = funchF_secret(rhoArgs.xVals);
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowAsym( ra, z );
	%
	vecZ0 = [ x0 + 0.1*randn(); 0.5; 0.3 ];
	prm = [];
	%prm.useCustomOmega = true;
	%prm.funchOmega      = @(rho)( sum(log( 1.0 + rho.^2 )) );
	%prm.funchOmegaDRho  = @(rho)( 2.0*rho ./ ( 1.0 + rho.^2 ) );
	%prm.funchOmegaDRho2 = @(rho)(diag( 2.0*(1.0-rho.^2)./((1.0+rho.^2).^2) ));
	%prm.matZBounds = [ -Inf, +Inf; 0.0, +Inf; 0.0, +Inf ];
	%
	if (use6723760Tweaks)
		vecZ0 = [ -1.2; 2.0; 2.0 ]
		%matZBounds = [ -Inf, +Inf; 0.0, +Inf; 0.0, +Inf ]
		%matZBounds = [ -1.321, -0.912; 0.0, +Inf; 0.0, +Inf ]
		matZBounds = [ -1.321, -0.912; -Inf, +Inf; -Inf, +Inf ]
		prm.matZBounds = matZBounds;
	end
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
	[ errFlag, vecRho, bigF0, bigFL, bigFR ] = absPowAsym_calcStuff( rhoArgs, vecZ );
	assert( ~errFlag );
	funchFModel = @(x)( bigF0 ...
	  + bigFL*(x<vecZ(1)).*(abs( x - vecZ(1) ).^vecZ(2)) ...
	  + bigFR*(x>vecZ(1)).*(abs( x - vecZ(1) ).^vecZ(3)) );
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
	[ errFlag, vecRho, bigF0, bigFL, bigFR ] = absPowAsym_calcStuff( rhoArgs, vecZ );
	assert( ~errFlag );
	funchFModel = @(x)( bigF0 ...
	  + bigFL*(x<vecZ(1)).*(abs( x - vecZ(1) ).^vecZ(2)) ...
	  + bigFR*(x>vecZ(1)).*(abs( x - vecZ(1) ).^vecZ(3)) );
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
