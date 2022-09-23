function [ funchFOfX, vecX0, datOut ] = genFunchAPS2022_fromType( strProbType, bigN0, probSeed=0, prm=[] )
	backup_prngStateDat = getprngstatedat();
	%
	setprngstates(probSeed,printState=false);
	%
	switch ( tolower(strProbType) )
	case { "lintest1" }
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
		gfaPrm.c2 = 0.0;
		gfaPrm.c3 = 0.0;
		gfaPrm.s1 = 0.0;
		gfaPrm.s2 = 0.0;
		gfaPrm.s3 = 0.0;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
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
		gfaPrm.c1 = 0.5 * abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1 * 1.0 * abs(randn());
		gfaPrm.c3 = gfaPrm.c1 * 0.5 * abs(randn());
		gfaPrm.s1 = 0.5 * abs(randn());
		gfaPrm.s2 = gfaPrm.s1 * 1.0 * abs(randn());
		gfaPrm.s3 = gfaPrm.s1 * 0.5 * abs(randn());
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "test_no_i" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN);
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0 * abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1 * 1.0 * abs(randn());
		gfaPrm.c3 = gfaPrm.c1 * 0.5 * abs(randn());
		gfaPrm.s1 = 1.0 * abs(randn());
		gfaPrm.s2 = gfaPrm.s1 * 1.0 * abs(randn());
		gfaPrm.s3 = gfaPrm.s1 * 0.5 * abs(randn());
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "test_sja0" }
		bigN = bigN0;
		bigM = bigN;
		funchSeed = probSeed;
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = 2.0;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 1.0;
		gfaPrm.c1 = 0.0;
		gfaPrm.c2 = 0.0;
		gfaPrm.c3 = 0.0;
		gfaPrm.s1 = 1.0;
		gfaPrm.s2 = 0.1;
		gfaPrm.s3 = 0.1;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "test_sja1" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN);
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 0.0;
		gfaPrm.c2 = 0.0;
		gfaPrm.c3 = 0.0;
		gfaPrm.s1 = 1.0;
		gfaPrm.s2 = 0.1;
		gfaPrm.s3 = 0.1;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja0" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN);
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0e-5*abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1*0.1;
		gfaPrm.c3 = gfaPrm.c1*0.1;
		gfaPrm.s1 = 1.0*abs(randn());
		gfaPrm.s2 = gfaPrm.s1*0.1*abs(randn());
		gfaPrm.s3 = gfaPrm.s1*0.1*abs(randn());
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja100" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN);
		gfaPrm.cx = abs(randn());
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0e-6*abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1*0.1*abs(randn());
		gfaPrm.c3 = gfaPrm.c1*0.1*abs(randn());
		gfaPrm.s1 = 1.0*abs(randn());
		gfaPrm.s2 = gfaPrm.s1*0.1*abs(randn());
		gfaPrm.s3 = gfaPrm.s1*0.1*abs(randn());
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja200" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN)/2.0;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0e-6*abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1*0.1;
		gfaPrm.c3 = gfaPrm.c1*0.1;
		gfaPrm.s1 = 1.0*abs(randn());
		gfaPrm.s2 = gfaPrm.s1*0.1;
		gfaPrm.s3 = gfaPrm.s1*0.1;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja200lin" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = sqrt(bigN)/2.0;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0e-6*abs(randn())/bigN;
		gfaPrm.c2 = 0.0;
		gfaPrm.c3 = 0.0;
		gfaPrm.s1 = 1.0*abs(randn());
		gfaPrm.s2 = 0.0;
		gfaPrm.s3 = 0.0;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja300" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = ceil( sqrt(bigN)/2.0 );
		gfaPrm.forceEqualNumberLambdaPerRow = true;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 0.0%1.0e-7/bigN;
		gfaPrm.c2 = gfaPrm.c1*0.1;
		gfaPrm.c3 = gfaPrm.c1*0.1;
		gfaPrm.s1 = 1.0;
		gfaPrm.s2 = gfaPrm.s1*0.1;
		gfaPrm.s3 = gfaPrm.s1*0.1;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja400" }
		bigN = ceil( bigN0*(2.0+rand())/2.0 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = ceil( sqrt(bigN) );
		gfaPrm.forceEqualNumberLambdaPerRow = true;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0e-8*abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1*0.1;
		gfaPrm.c3 = gfaPrm.c1*0.1;
		gfaPrm.s1 = 1.0*abs(randn());
		gfaPrm.s2 = gfaPrm.s1*0.1;
		gfaPrm.s3 = gfaPrm.s1*0.1;
		[ funchFOfX, gfaDatOut ] = genFunchAPS2022( bigN, funchSeed, gfaPrm );
		vecX0 = zeros(bigN,1);
	case { "sja500" }
		bigN = ceil( bigN0*(1.0+rand())/1.5 );
		bigM = bigN;
		funchSeed = round( 1.0E12*rand() );
		gfaPrm = [];
		gfaPrm.bigM = bigM;
		gfaPrm.bigP = min([ ceil(5.0*sqrt(bigN)), bigN, bigM ]);
		gfaPrm.lambda = 1 + sqrt(bigN)*abs(randn());
		gfaPrm.forceEqualNumberLambdaPerRow = false;
		gfaPrm.cx = 1.0;
		gfaPrm.c0 = 0.0;
		gfaPrm.c1 = 1.0e-10*abs(randn())/bigN;
		gfaPrm.c2 = gfaPrm.c1*0.1*abs(randn());
		gfaPrm.c3 = gfaPrm.c1*0.1*abs(randn());
		gfaPrm.s1 = 1.0*abs(randn());
		gfaPrm.s2 = gfaPrm.s1*0.1*abs(randn());
		gfaPrm.s3 = gfaPrm.s1*0.1*abs(randn());
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
