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
	rhoArgs.xVals = linspace( -1.0, 2.0, 5 );
	rhoArgs.fVals = f0+f1*abs(rhoArgs.xVals-x0).^p;
	rhoArgs.dVals = 0.4*ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
	%
	vecZ0 = [ 0.0, 2.0 ];
	prm = [];
	prm.useCustomOmega = true;
	prm.funchOmega = @(rho)( sum(log( 1.0 + rho.^2 )) );
	prm.funchOmegaP = @(rho)( 2.0*rho ./ ( 1.0 + rho.^2 ) );
	prm.funchOmegaPP = @(rho)(diag( 2.0*(1.0-rho.^2)./((1.0+rho.^2).^2) ));
	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm );
	echo__vecZ = vecZ
	echo__retCode = retCode;
	echo__datOut = datOut;
	assert( retCode == RETCODE__SUCCESS );
	msg( thisFile, __LINE__, "Finished test." );
%
