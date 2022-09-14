function [ funchFOfX, vecX0, datOut ] = zrun_genProb( strProbType, bigN0, probSeed=0, prm=[] )
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
		gfaPrm.c0 = 1.0;
		gfaPrm.c1 = 1.0E0/bigN;
		gfaPrm.c2 = gfaPrm.c1 * 1.0E0;
		gfaPrm.c3 = gfaPrm.c1 * 1.0E0;
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
	case { "test1" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN);
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 1.0;
		gfaPrm.c1 = 1.0E-1*abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1 * abs(randn());
		gfaPrm.c3 = gfaPrm.c1 * abs(randn());
		gfaPrm.s1 = 1.0E-1 * abs(randn());
		gfaPrm.s2 = gfaPrm.s1 * abs(randn());
		gfaPrm.s3 = gfaPrm.s1 * abs(randn());
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	otherwise
		setprngstatedat(backup_prngStateDat);
		error([ "Invalid strProbType (\"" strProbType "\")." ]);
	endswitch
	%
	datOut.sizeX = bigN;
	datOut.sizeF = bigM;
	datOut.funchSeed = funchSeed;
	datOut.genFunchPrm = gfaPrm;
	datOut.genFunchDatOut = gfaDatOut;
	setprngstatedat(backup_prngStateDat);
return;
endfunction
