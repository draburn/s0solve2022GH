	clear;
	commondefs;
	thisFile = "test2_findBestFit1D";
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(0);
	f0 = 0.2;
	f1 = 0.4;
	x0 = 0.6;
	p = 3.5;
	rhoArgs.xVals = linspace( -1.0, 2.0, 5 );
	rhoArgs.fVals = f0+f1*abs(rhoArgs.xVals-x0).^p;
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
	%
	vecZ0 = [ 0.0, 2.0 ];
	prm = [];
	prm.useCustomOmega = true;
	prm.funchOmega = @(rho)( sum(rho.^2) );
	prm.funchOmegaP = @(rho)( 2.0*rho );
	%prm.funchOmegaPP = @(rho)( 2.0*eye(max(size(rho)),max(size(rho))) );
	prm.funchOmegaPP = @(rho)( diag(2.0*ones(size(rho))) );
	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm );
	echo__vecZ = vecZ
	echo__retCode = retCode;
	echo__datOut = datOut;
	assert( retCode == RETCODE__SUCCESS );
	msg( thisFile, __LINE__, "Finished test." );