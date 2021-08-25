	clear;
	commondefs;
	thisFile = "test3_findBestFit1D";
	setprngstates(0);
	msg( thisFile, __LINE__, "Starting test." );
	%
	rhoArgs.xVals = linspace( 0.0, 2.0, 3 );
	rhoArgs.fVals = abs(rhoArgs.xVals).^3;
	rhoArgs.dVals = ones(size(rhoArgs.xVals));
	funchRho = @(ra,z) calcRho_absPowSym( ra, z );
	%
	prm = [];
	echo__prm = prm
	vecZ0 = [ 0.0, 2.0 ]
	[ vecZ, retCode, datOut ] = findBestFit1D( funchRho, rhoArgs, vecZ0, prm=[] );
	echo__vecZ = vecZ
	echo__retCode = retCode;
	echo__datOut = datOut;
	assert( retCode == RETCODE__SUCCESS );
	msg( thisFile, __LINE__, "Finished test." );
%