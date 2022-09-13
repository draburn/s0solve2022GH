% Function...

function [ funchFOfX, vecX0, datOut ] = zcompare_genProb( strProbType, bigN0, probSeed=0, prm=[] )
	backup_prngStateDat = getprngstatedat();
	%
	setprngstates(probSeed,printState=false);
	%
	switch ( tolower(strProbType) )
	case { "test0" }
		bigN = bigN0;
		bigM = bigN;
		funchSeed = probSeed;
		gfaPrm = [];
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
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "testr" }
		bigN = ceil( bigN0*exp(rand()) );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
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
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	otherwise
		setprngstatedat(backup_prngStateDat);
		error([ "Invalid strProbType (\"" strProbType "\")." ]);
	endswitch
	%
	datOut.bigN = bigN;
	datOut.funchSeed = funchSeed;
	datOut.genFunchPrm = gfaPrm;
	datOut.genFunchDatOut = gfaDatOut;
	setprngstatedat(backup_prngStateDat);
return;
endfunction
