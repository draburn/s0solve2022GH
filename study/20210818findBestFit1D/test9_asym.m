	clear;
	commondefs;
	thisFile = "test9_asym";
	numFigs = 0;
	%
	%
	%
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(19597392); %N-I-C-E!
	%setprngstates(15638000); % Bad as-is, but easy to fix.
	use15638000Hacks = false;
	%setprngstates(4406128); % Screwy pts, fit is meh but wrong.
	%setprngstates(45081616); % Fails because X0 moves too far.
	use45081616Hacks = false;
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
	numPts = 6 + round(abs(4.0*randn()))
	rhoArgs.xVals = x0 + randn(numPts,1);
	rhoArgs.xVals(1:3) = x0 - abs(randn(3,1));
	rhoArgs.xVals(4:6) = x0 + abs(randn(3,1));
	%numPts = 10;
	%rhoArgs.xVals = x0+linspace( -2.0, 2.0, numPts )';
	rhoArgs.xVals = sort(rhoArgs.xVals);
	numPts = size(rhoArgs.xVals,1);
	rhoArgs.fVals = funchF_secret(rhoArgs.xVals);
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowAsym( ra, z );
	%
	prm = [];
	vecZ0 = [ x0 + 0.1*randn(); 2.0; 2.0 ];
	%prm.useCustomOmega = true;
	%prm.funchOmega      = @(rho)( sum(log( 1.0 + rho.^2 )) );
	%prm.funchOmegaDRho  = @(rho)( 2.0*rho ./ ( 1.0 + rho.^2 ) );
	%prm.funchOmegaDRho2 = @(rho)(diag( 2.0*(1.0-rho.^2)./((1.0+rho.^2).^2) ));
	prm.matZBounds = [ -Inf, +Inf; 0.0, +Inf; 0.0, +Inf ];
	if (use45081616Hacks)
		prm.matZBounds(1,:) = [ -0.9, -0.15 ]
	end
	if (use15638000Hacks)
		%matZBounds = [ -Inf, +Inf; 0.1, +Inf; 0.1, +Inf ]
		%prm.matZBounds = matZBounds
		%vecZ0(2:3) = 0.5;
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
