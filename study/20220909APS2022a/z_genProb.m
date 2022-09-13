% Function...

function [ funchFOfX, datOut ] = z_genProb( strProbType, bigN0, probSeed=0, prm=[] )
	backup_prngStateDat = getprngstatedat();
	%
	setprngstates(probSeed,printState=false);
	gfaPrm = [];
	%
	switch ( tolower(strProbType) )
	case { "test0" }
		bigN = bigN0;
		bigM = bigN;
		funchSeed = probSeed;
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = 0.0;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0;
		gfaPrm.c2 = gfaPrm.c1 * 1.0E-2;
		gfaPrm.c3 = gfaPrm.c1 * 1.0E-2;
		gfaPrm.s1 = 0.0;
		gfaPrm.s2 = gfaPrm.s1 * 0.0;
		gfaPrm.s3 = gfaPrm.s1 * 0.0;
	case { "testr" }
		bigN = ceil( bigN0*exp(rand()) );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = 1.0 + 0.1*bigN*rand();
		gfaPrm.cx = exp(3.0*randn());
		gfaPrm.c0 = exp(3.0*randn());
		gfaPrm.c1 = exp(3.0*randn());
		gfaPrm.c2 = gfaPrm.c1 * exp(-abs(randn()));
		gfaPrm.c3 = gfaPrm.c1 * exp(-abs(randn()));
		gfaPrm.s1 = exp(3.0*randn());
		gfaPrm.s2 = gfaPrm.s1 * exp(-abs(randn()));
		gfaPrm.s3 = gfaPrm.s1 * exp(-abs(randn()));
	otherwise
		setprngstatedat(backup_prngStateDat);
		error([ "Invalid strProbType (\"" strProbType "\")." ]);
	endswitch
	%
	[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
	datOut.bigN = bigN;
	datOut.funchSeed = funchSeed;
	datOut.gfaPrm = gfaPrm;
	datOut.gfaDatOut = gfaDatOut;
	setprngstatedat(backup_prngStateDat);
return;
endfunction
