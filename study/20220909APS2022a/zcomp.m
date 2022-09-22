function zcompDatOut = zcomp( probSetPrm=[], algoSetPrm=[], default_solverPrm=[], zcompPrm=[] )
	backup_prngStateDat = getprngstatedat();
	mydefs;
	startTime = time();
	if ( isempty(probSetPrm) )
		probSetPrm.numProbs = 10;
		probSetPrm.probType = "test1";
		probSetPrm.numUnknowns = 20;
		probSetPrm.setSeed = 0;
		zcompPrm.verbLev = mygetfield( zcompPrm, "verbLev", VERBLEV__INFO );
	endif
	if ( isempty(algoSetPrm) )
		algoSetPrm.n = 3;
		algoSetPrm.s(1).f = @groot_jfnk_basic;
		algoSetPrm.s(1).p.btCoeff = 0.0;
		algoSetPrm.s(2).f = @groot_jfnk_basic;
		algoSetPrm.s(3).f = @groot_fsolve;
	endif
	zcompPrm.verbLev = mygetfield( zcompPrm, "verbLev", VERBLEV__INFO );
	zcompPrm.valdLev = mygetfield( zcompPrm, "valdLev", VALDLEV__HIGH );
	%
	dateStr = datestr(now,31);
	dateStr(" "==dateStr) = "_";
	dateStr("-"==dateStr) = "";
	dateStr(":"==dateStr) = "";
	runName = [ ...
	  num2str(probSetPrm.numProbs) "X_" probSetPrm.probType ...
	  "_SZ" num2str(probSetPrm.numUnknowns) ...
	  "_SD" num2str(probSetPrm.setSeed) ...
	  "__" dateStr ];
	if ( zcompPrm.verbLev >= VERBLEV__INFO )
		msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
		msg( __FILE__, __LINE__, sprintf( "Starting runName \"%s\".", runName ) );
		msg( __FILE__, __LINE__, "" );
	endif
	%
	%
	setprngstates( probSetPrm.setSeed, false );
	probSeeds = floor( 1E8*rand(1,probSetPrm.numProbs) );
	setprngstatedat(backup_prngStateDat); % Do this ASAP, in case of abrupt return.
	zcompDatOut.probSeeds = probSeeds;
	%
	assert( isposintscalar(probSetPrm.numProbs) );
	assert( isposintscalar(probSetPrm.numUnknowns) );
	zcompDatOut.probSetPrm = probSetPrm;
	zcompDatOut.algoSetPrm = algoSetPrm;
	zcompDatOut.default_solverPrm = default_solverPrm;
	zcompDatOut.zcompPrm = zcompPrm;
	zcompDatOut.runName = runName;
	%
	doCharLog = (  ( zcompPrm.verbLev >= VERBLEV__INFO ) && ( zcompPrm.verbLev < VERBLEV__PROGRESS )  );
	for probIndex = 1 : probSetPrm.numProbs
		if ( stopsignalpresent() )
			msg( __FILE__, __LINE__, "Found stopsignal." );
			break;
		endif
		this_probType = probSetPrm.probType;
		this_sizeXIsh = probSetPrm.numUnknowns;
		this_probSeed = probSeeds(probIndex);
		this_probGenPrm = [];
		%
		zcompDatOut.prob(probIndex).probType = this_probType;
		zcompDatOut.prob(probIndex).sizeXIsh = this_sizeXIsh;
		zcompDatOut.prob(probIndex).probSeed = this_probSeed;
		zcompDatOut.prob(probIndex).probGenPrm = this_probGenPrm;
		[ this_funchF, this_vecX0, this_probDat ] = genFunchAPS2022_fromType( this_probType, this_sizeXIsh, this_probSeed, this_probGenPrm );
		assert( is_function_handle(this_funchF) );
		assert( isrealarray(this_vecX0,[this_probDat.sizeX,1]) );
		this_probStr = sprintf( "Prob %d", probIndex);
		this_probSizeStr = sprintf( "(%dx%d)", this_probDat.sizeF, this_probDat.sizeX );
		if ( doCharLog )
			msgnnl( __FILE__, __LINE__, sprintf( "  %15s  %15s:  ", this_probStr, this_probSizeStr ) );
		elseif ( zcompPrm.verbLev >= VERBLEV__PROGRESS )
			msg( __FILE__, __LINE__, sprintf( " %15s  %15s....", this_probStr, this_probSizeStr ) );
		endif
		zcompDatOut.prob(probIndex).sizeX = this_probDat.sizeX;
		zcompDatOut.prob(probIndex).sizeF = this_probDat.sizeF;
		%
		%
		%if ( 1 )
		%	foo = this_probDat.genFunchDatOut;
		%	matJSecret = foo.matA11L*foo.matA11R + foo.matA12L*foo.matA12R + foo.matB1;
		%	figure(200);
		%	imagesc(matJSecret);
		%endif
		%
		%
		this_xPrm = [];
		if ( doCharLog )
			this_xPrm.verbLev = VERBLEV__INFO;
		endif
		zcompDatOut.prob(probIndex).xPrm = this_xPrm;
		this_grootXDatOut = groot_x( this_funchF, this_vecX0, algoSetPrm, default_solverPrm, this_xPrm );
		zcompDatOut.prob(probIndex).grootXDatOut = this_grootXDatOut;
		%
		if ( doCharLog )	
			printf( "\n" );
		endif
		clear this_*;
	endfor
	if ( zcompPrm.verbLev >= VERBLEV__INFO )
		numAlgos = mygetfield( algoSetPrm, "n", -1 );
		msg( __FILE__, __LINE__, sprintf( "Completed %dx%d algo(solves) in %gs.", probSetPrm.numProbs, numAlgos, time()-startTime )  );
	endif
	zcompDatOut.elapsedTime = time()-startTime;
	%
	if ( zcompPrm.verbLev >= VERBLEV__INFO )
		msg( __FILE__, __LINE__, "" );
		msg( __FILE__, __LINE__, sprintf( "Finished runName \"%s\".", runName ) );
		msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
		msg( __FILE__, __LINE__, sprintf( "Elapsed time is %gs.", time()-startTime ) );
	endif
	%
	setprngstatedat(backup_prngStateDat); % May be redundant.
return;
endfunction
